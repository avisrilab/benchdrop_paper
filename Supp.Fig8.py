"""Scrublet doublet calling on the PBMC long-read matrix (Methods 2.1, Supplementary Fig. 8).
threshold_minimum, ported from scikit-image, supplies the histogram-minimum threshold where
scikit-image is absent.

Stages: load the paper's 19,031 cells, Scrublet at two prior doublet rates, cluster before and
after removal, doublets per cell type, Supplementary Fig. 8, clustering seed baseline.
"""
from __future__ import annotations
import argparse
import gzip
import json
import numpy as np
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from scipy import ndimage as ndi
from time import perf_counter


FEED = Path(os.path.expandvars("$BENCHDROP_FEED"))
SILO = Path(os.path.expandvars("$BENCHDROP_RESULTS/pbmc_doublets"))

MATRIX_DIR = FEED / "pbmc/count_lr"          # the PBMC EM transcript matrix
# The paper's 19,031 PBMC cells with their Azimuth labels (predicted.celltype.l1 and .l2),
# exported from the LR.PBMC.S3.rds Seurat object; shipped input.
PAPER_CELLS_META = FEED / "pbmc/paper_cells_meta.tsv"
TX_MAP = FEED / "annotation_maps/gencode_v32.tx.tsv"      # tx -> gene/name/biotype
GENE_MAP = FEED / "annotation_maps/gencode_v32.gene.tsv"  # gene -> name/biotype

# Expected matrix header. Rows are BARCODES, columns are TRANSCRIPTS (the reverse of the
# cellranger orientation), so the matrix loads via scipy.io.mmread, not sc.read_10x_mtx.
EXPECT_MTX = (2_109_262, 199_138, 41_738_415)

BARS = {
    # ---- s01: every paper cell must be in the matrix.
    "cells_expected": 19_031,
    "cells_min_matched": 19_031,
    "min_txmap_frac": 0.99,             # tx->gene coverage floor

    # ---- s02: Scrublet at the headline prior (scanpy's default) and a sensitivity prior.
    "doublet_rate_headline": 0.05,
    "doublet_rate_sensitivity": 0.10,
    "random_state": 0,
    # scanpy's own floor inside sc.pp.scrublet (`pp.filter_cells(min_genes=3)`): cells below it
    # get no score; used to account for every unscored cell, never to relax anything.
    "scrublet_internal_min_genes": 3,
    "threshold_nbins": 256,             # skimage.filters.threshold_minimum's defaults
    "threshold_max_iter": 10_000,

    # ---- s03: before/after clustering, identical recipe on both sides and both arms.
    "cluster_n_top_genes": 2000,
    "cluster_n_pcs": 30,
    "cluster_n_neighbors": 15,
    "cluster_leiden_resolution": 0.5,
    "cluster_random_state": 0,
    # s06: Leiden re-run on the untouched before-graph at these seeds; the ARI against seed 0
    # is the algorithm's own reproducibility, the baseline the s03 before/after ARI is read against.
    "seed_baseline_seeds": [1, 2, 3],
}


# ---------- stage machinery ----------

@dataclass
class StageResult:
    metrics: dict = field(default_factory=dict)
    passed: bool | None = None      # None for a report-only stage
    note: str = ""


class Stage:
    name: str = ""
    bar: str = ""

    def run(self) -> StageResult:
        raise NotImplementedError


# ---------- histogram-minimum threshold ----------


class ThresholdMinimumError(RuntimeError):
    """Raised on the same two conditions the real skimage function raises on."""


def _find_local_maxima_idx(hist: np.ndarray) -> list[int]:
    """Verbatim port of skimage's inner `find_local_maxima_idx`."""
    maximum_idxs: list[int] = []
    direction = 1
    for i in range(hist.shape[0] - 1):
        if direction > 0:
            if hist[i + 1] < hist[i]:
                direction = -1
                maximum_idxs.append(i)
        else:
            if hist[i + 1] > hist[i]:
                direction = 1
    return maximum_idxs


def threshold_minimum(scores, nbins: int = 256, max_num_iter: int = 10_000) -> tuple[float, dict]:
    """Port of `skimage.filters.threshold_minimum(image, nbins, max_num_iter)`.

    Returns (threshold, diagnostics) where diagnostics carries nbins, the number of smoothing
    iterations actually run, and the number of maxima found at convergence.
    """
    scores = np.asarray(scores, dtype=float)
    counts, edges = np.histogram(scores, bins=nbins, range=(float(scores.min()), float(scores.max())))
    bin_centers = (edges[:-1] + edges[1:]) / 2.0

    smooth_hist = counts.astype("float32", copy=False)
    maximum_idxs: list[int] = []
    counter = 0
    for counter in range(max_num_iter):
        smooth_hist = ndi.uniform_filter1d(smooth_hist, 3)
        maximum_idxs = _find_local_maxima_idx(smooth_hist)
        if len(maximum_idxs) < 3:
            break

    diagnostics = {"nbins": nbins, "max_num_iter": max_num_iter,
                    "iterations": counter + 1, "n_maxima": len(maximum_idxs)}

    if len(maximum_idxs) != 2:
        raise ThresholdMinimumError(
            f"unable to find two maxima in histogram (found {len(maximum_idxs)} after "
            f"{counter + 1} smoothing passes)")
    if counter == max_num_iter - 1:
        raise ThresholdMinimumError("maximum iteration reached for histogram smoothing")

    lo, hi = maximum_idxs[0], maximum_idxs[1] + 1
    threshold_idx = int(np.argmin(smooth_hist[lo:hi]))
    threshold = float(bin_centers[maximum_idxs[0] + threshold_idx])
    return threshold, diagnostics



# ---------- shared helpers ----------

def _fresh(*paths: Path):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _made(*paths: Path):
    for p in paths:
        if not Path(p).exists() or Path(p).stat().st_size == 0:
            raise RuntimeError(f"stage exited but produced no {p}")


def _tx2gene() -> dict:
    tx2g = {}
    with open(TX_MAP) as fh:
        header = fh.readline()
        assert header.startswith("tx_id")
        for ln in fh:
            tx, g = ln.rstrip("\n").split("\t")[:2]
            tx2g[tx] = g
    return tx2g


def _gene_names() -> dict:
    g2name = {}
    with open(GENE_MAP) as fh:
        header = fh.readline()
        assert header.startswith("gene_id")
        for ln in fh:
            g, name = ln.rstrip("\n").split("\t")[:2]
            g2name[g] = name
    return g2name


def _load_cell_labels():
    """The paper's cell list in file order, with barcode -> predicted.celltype.l1 (primary) and
    barcode -> predicted.celltype.l2 (secondary)."""
    order: list = []
    primary: dict = {}
    secondary: dict = {}
    with open(PAPER_CELLS_META) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bi = header.index("barcode")
        l1i = header.index("predicted.celltype.l1")
        l2i = header.index("predicted.celltype.l2")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            b = f[bi]
            order.append(b)
            primary[b] = f[l1i]
            secondary[b] = f[l2i]
    return order, primary, secondary


# Fixed colors for the predicted.celltype.l1 labels; an unseen label draws from the spare palette.
_KNOWN_CELLTYPE_COLORS = {
    "CD4 T": "#2ca02c", "CD8 T": "#98df8a", "B": "#1f77b4", "Mono": "#d62728",
    "DC": "#ff7f0e", "NK": "#9467bd", "other T": "#8c564b", "other": "#7f7f7f",
}
_SPARE_PALETTE = ["#bcbd22", "#17becf", "#aec7e8", "#e377c2", "#c49c94", "#f7b6d2"]


def _celltype_colors(labels) -> dict:
    """label -> hex color, stable for every known label, spare palette for unseen ones."""
    spare = list(_SPARE_PALETTE)
    out = {}
    for ct in sorted(set(labels)):
        out[ct] = _KNOWN_CELLTYPE_COLORS.get(ct) or (spare.pop(0) if spare else "#000000")
    return out


def _cluster(adata):
    """normalize -> log1p -> HVG -> PCA -> neighbours -> leiden -> UMAP, always re-derived from
    raw counts (`layers['counts']`), never by sub-setting a fitted embedding."""
    import scanpy as sc

    a = adata.copy()
    a.X = a.layers["counts"].copy()
    a.uns.pop("log1p", None)   # X was just restored from raw counts; clear any stale flag
    b = BARS
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=b["cluster_n_top_genes"])
    sc.pp.pca(a, n_comps=b["cluster_n_pcs"], mask_var="highly_variable",
              random_state=b["cluster_random_state"])
    sc.pp.neighbors(a, n_neighbors=b["cluster_n_neighbors"], random_state=b["cluster_random_state"])
    sc.tl.leiden(a, resolution=b["cluster_leiden_resolution"], flavor="igraph", n_iterations=2,
                 directed=False, random_state=b["cluster_random_state"])
    sc.tl.umap(a, random_state=b["cluster_random_state"])
    return a


# ---------- s01: load the paper's cells, aggregate transcripts to genes ----------

class S01Load(Stage):
    name, bar = "s01_load", ("mtx header == expected dims; paper cells matched >= "
                             "cells_min_matched; tx->gene coverage >= min_txmap_frac")
    asks = "Do the paper's 19,031 cells load cleanly off this matrix, at gene resolution?"

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        import scipy.io as sio
        import scipy.sparse as sp

        SILO.mkdir(parents=True, exist_ok=True)
        out = SILO / "paper_gene.h5ad"
        _fresh(out)

        with gzip.open(MATRIX_DIR / "matrix.mtx.gz", "rt") as fh:
            fh.readline()
            fh.readline()
            dims = tuple(int(x) for x in fh.readline().split())
        if dims != EXPECT_MTX:
            return StageResult(metrics={"mtx_header": dims}, passed=False,
                               note=f"header {dims} != expected {EXPECT_MTX}")

        with gzip.open(MATRIX_DIR / "barcodes.tsv.gz", "rt") as fh:
            bcs = [l.strip() for l in fh]
        with gzip.open(MATRIX_DIR / "features.tsv.gz", "rt") as fh:
            feats = [l.strip() for l in fh]

        order, primary, secondary = _load_cell_labels()
        n_cells_in_file = len(order)

        bc_set = set(bcs)
        missing = [b for b in order if b not in bc_set]
        bc_idx = {b: i for i, b in enumerate(bcs)}
        rows = [bc_idx[b] for b in order if b in bc_idx]
        matched = len(rows)

        m = sio.mmread(MATRIX_DIR / "matrix.mtx.gz").tocsr()   # barcodes x transcripts
        sub = m[rows, :]
        obs_bc = [bcs[i] for i in rows]
        obs_ct = [primary[b] for b in obs_bc]
        obs_ct2 = [secondary[b] for b in obs_bc]

        tx2g = _tx2gene()
        g2name = _gene_names()
        mapped = sum(1 for t in feats if t in tx2g)
        txmap_frac = mapped / len(feats)

        genes = sorted({tx2g[t] for t in feats if t in tx2g})
        gidx = {g: j for j, g in enumerate(genes)}
        col = np.array([gidx.get(tx2g.get(t), -1) for t in feats])
        keep = col >= 0
        ind = sp.csr_matrix((np.ones(int(keep.sum())), (np.nonzero(keep)[0], col[keep])),
                            shape=(len(feats), len(genes)))
        gsub = (sub @ ind).tocsr()

        obs = pd.DataFrame({"celltype": pd.Categorical(obs_ct),
                            "celltype_l2": pd.Categorical(obs_ct2)}, index=obs_bc)
        var = pd.DataFrame({"gene_name": [g2name.get(g, g) for g in genes]}, index=genes)
        adata = ad.AnnData(gsub.astype(np.float32), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()
        adata.write_h5ad(out)
        _made(out)

        n_genes = len(genes)
        ok = (dims == EXPECT_MTX and matched >= BARS["cells_min_matched"]
              and txmap_frac >= BARS["min_txmap_frac"])
        metrics = {"mtx_header": list(dims), "cells_in_file": n_cells_in_file, "cells_matched": matched,
                   "cells_missing_from_matrix": len(missing),
                   "missing_barcodes_sample": missing[:20],
                   "n_transcripts": len(feats), "n_genes": n_genes,
                   "txmap_frac": round(txmap_frac, 4)}
        note = (f"{matched:,}/{n_cells_in_file:,} cells matched ({len(missing)} of "
                f"the cell list absent from this matrix's barcodes), {n_genes:,} genes, "
                f"txmap {txmap_frac:.4f}")
        return StageResult(metrics=metrics, passed=ok, note=note)


# ---------- s02: Scrublet at two expected-doublet-rate arms ----------

class S02Scrublet(Stage):
    # scanpy's sc.pp.scrublet drops cells with fewer than 3 detected genes before scoring, so
    # those come back with no score. Every unscored cell must be one of those, and the unscored
    # set must be identical across both arms; they are then dropped from s02 onward.
    name, bar = "s02_scrublet", ("Scrublet returns a finite doublet_score for every loaded cell "
                                 "at both prior rates, except cells with <3 detected genes, and "
                                 "that excluded set is identical across both arms; the threshold "
                                 "and its source (skimage, or the ported threshold_minimum) are "
                                 "recorded per arm")
    asks = "What fraction of cells does Scrublet flag at a prior doublet rate of 0.05 and of 0.10?"

    # Placeholder threshold for the first sc.pp.scrublet call when scikit-image is absent (the
    # `threshold=None` branch imports skimage). predicted_doublet from that call is discarded;
    # only doublet_score and the simulated-score histogram are kept.
    _PLACEHOLDER_THRESHOLD = 10.0

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        import scanpy as sc

        try:
            import skimage  # noqa: F401
            skimage_available = True
        except ModuleNotFoundError:
            skimage_available = False

        d = SILO / "scrublet"
        d.mkdir(parents=True, exist_ok=True)
        base = ad.read_h5ad(SILO / "paper_gene.h5ad")
        n_cells = base.n_obs

        # Which cells fall below scanpy's internal min_genes floor, from the loaded counts.
        n_genes_detected = np.asarray((base.layers["counts"] > 0).sum(axis=1)).ravel()
        below_floor = set(base.obs_names[n_genes_detected < BARS["scrublet_internal_min_genes"]])

        arms = {"headline": BARS["doublet_rate_headline"],
                "sensitivity": BARS["doublet_rate_sensitivity"]}
        metrics, problems = {}, []
        excluded_by_arm = {}

        for arm, rate in arms.items():
            tsv = d / f"{arm}.tsv"
            _fresh(tsv)
            a = base.copy()

            if skimage_available:
                # scanpy's own automatic threshold.
                sc.pp.scrublet(a, expected_doublet_rate=rate,
                               random_state=BARS["random_state"], threshold=None, verbose=False)
                threshold = float(a.uns["scrublet"]["threshold"])
                diag = {"source": "skimage", "nbins": None, "iterations": None, "n_maxima": None}
            else:
                # Pass 1: explicit placeholder threshold, never touches the skimage import;
                # yields doublet_score (obs) and doublet_scores_sim (uns) for every cell.
                sc.pp.scrublet(a, expected_doublet_rate=rate,
                               random_state=BARS["random_state"],
                               threshold=self._PLACEHOLDER_THRESHOLD, verbose=False)
                sim_scores = np.asarray(a.uns["scrublet"]["doublet_scores_sim"])
                try:
                    threshold, tdiag = threshold_minimum(
                        sim_scores, nbins=BARS["threshold_nbins"],
                        max_num_iter=BARS["threshold_max_iter"])
                except ThresholdMinimumError as exc:
                    problems.append(f"{arm}: ported threshold_minimum failed ({exc})")
                    continue
                diag = {"source": "ported", **tdiag}
                # Apply the threshold as `Scrublet.call_doublets` does (`doublet_score > threshold`).
                a.obs["predicted_doublet"] = (a.obs["doublet_score"].to_numpy() > threshold)
                a.uns["scrublet"]["threshold"] = threshold

            score = a.obs["doublet_score"].to_numpy()
            pred = a.obs["predicted_doublet"].to_numpy()
            finite = np.isfinite(score)
            nan_bcs = set(np.asarray(a.obs_names)[~finite])
            excluded_by_arm[arm] = nan_bcs
            unexplained = nan_bcs - below_floor
            if len(score) != n_cells:
                problems.append(f"{arm}: {len(score)} scores vs {n_cells} loaded cells")
            if unexplained:
                problems.append(f"{arm}: {len(unexplained)} unscored cell(s) NOT explained by "
                                f"the <{BARS['scrublet_internal_min_genes']}-genes floor: "
                                f"{sorted(unexplained)[:10]}")

            keep = finite
            n_scored = int(keep.sum())
            n_doub = int(pred[keep].sum())
            rate_obs = n_doub / n_scored if n_scored else float("nan")

            pd.DataFrame({"barcode": np.asarray(a.obs_names)[keep],
                          "celltype": a.obs["celltype"].astype(str).to_numpy()[keep],
                          "celltype_l2": a.obs["celltype_l2"].astype(str).to_numpy()[keep],
                          "doublet_score": score[keep], "predicted_doublet": pred[keep],
                          }).to_csv(tsv, sep="\t", index=False)
            _made(tsv)

            metrics[arm] = {"expected_doublet_rate": rate, "n_cells_loaded": n_cells,
                            "n_cells_scored": n_scored,
                            "excluded_unscored": sorted(nan_bcs),
                            "n_doublets": n_doub, "doublet_rate": round(rate_obs, 4),
                            "threshold": round(float(threshold), 6), "threshold_diag": diag}

        if len(arms) == 2 and not problems:
            a1, a2 = excluded_by_arm.get("headline"), excluded_by_arm.get("sensitivity")
            if a1 is not None and a2 is not None and a1 != a2:
                problems.append(f"unscored-cell set differs between arms: headline={sorted(a1)} "
                                f"vs sensitivity={sorted(a2)}")

        out = SILO / "scrublet_summary.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        ok = not problems and len(metrics) == len(arms)
        note = ("; ".join(problems) if problems else
                "; ".join(f"{a}: {v['n_doublets']}/{v['n_cells_scored']} "
                          f"({v['doublet_rate']*100:.2f}%) of {v['n_cells_scored']} scored "
                          f"({v['n_cells_loaded'] - v['n_cells_scored']} excluded, <3 genes), "
                          f"threshold {v['threshold']:.4f}" for a, v in metrics.items()))
        return StageResult(metrics=metrics, passed=ok, note=note)


# ---------- s03: before/after clustering + ARI, both arms ----------

class S03Cluster(Stage):
    name, bar = "s03_cluster", ("report-only: HVG/PCA/neighbours/leiden/UMAP re-derived before "
                                "and after doublet removal, identical recipe both sides and both "
                                "arms; ARI, cluster counts, and per-cluster removal reported")
    asks = ("How much does removing Scrublet-called doublets change the retained-cell cluster "
            "structure, at each expected-doublet-rate arm?")

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        from sklearn.metrics import adjusted_rand_score

        d = SILO
        before_h5 = d / "before.h5ad"
        _fresh(before_h5)
        full = ad.read_h5ad(d / "paper_gene.h5ad")

        # A cell with no doublet call cannot enter a doublet-removal compare, so "before" is the
        # scored set (identical on both arms, which s02 checks).
        scored_bcs = pd.read_csv(d / "scrublet" / "headline.tsv", sep="\t")["barcode"].tolist()
        n_dropped_unscored = full.n_obs - len(scored_bcs)
        full = full[scored_bcs].copy()

        before = _cluster(full)
        before.write_h5ad(before_h5)
        _made(before_h5)

        n_clusters_before = int(before.obs["leiden"].cat.categories.size)
        metrics = {"cells_before": int(before.n_obs), "n_clusters_before": n_clusters_before,
                  "cells_dropped_unscored": int(n_dropped_unscored),
                  "recipe": {k: BARS[k] for k in
                             ("cluster_n_top_genes", "cluster_n_pcs", "cluster_n_neighbors",
                              "cluster_leiden_resolution", "cluster_random_state")}}

        arms = {}
        for arm in ("headline", "sensitivity"):
            calls = pd.read_csv(d / "scrublet" / f"{arm}.tsv", sep="\t").set_index("barcode")
            pred = calls.loc[before.obs_names, "predicted_doublet"].to_numpy().astype(bool)
            keep = ~pred
            n_removed = int(pred.sum())

            after_h5 = d / f"after_{arm}.h5ad"
            _fresh(after_h5)
            after = _cluster(before[keep].copy())
            after.write_h5ad(after_h5)
            _made(after_h5)

            ari = float(adjusted_rand_score(before.obs["leiden"][keep].astype(str),
                                            after.obs["leiden"].astype(str)))
            n_clusters_after = int(after.obs["leiden"].cat.categories.size)

            loss = {}
            for cl in before.obs["leiden"].cat.categories:
                m = (before.obs["leiden"] == cl).to_numpy()
                n_cl = int(m.sum())
                n_cl_removed = int((m & pred).sum())
                loss[str(cl)] = {"n": n_cl, "removed": n_cl_removed,
                                 "frac_removed": round(n_cl_removed / n_cl, 4) if n_cl else 0.0}
            top_loss = sorted(loss.items(), key=lambda kv: -kv[1]["frac_removed"])[:5]

            arms[arm] = {"n_removed": n_removed, "cells_after": int(after.n_obs),
                         "n_clusters_after": n_clusters_after, "ari": round(ari, 4),
                         "top_clusters_by_frac_removed": top_loss}

        metrics["arms"] = arms
        out = d / "cluster_stability.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        note = "; ".join(f"{k}: ARI {v['ari']:.3f}, clusters {n_clusters_before}->"
                         f"{v['n_clusters_after']}, removed {v['n_removed']}"
                         for k, v in arms.items())
        return StageResult(metrics=metrics, passed=None, note=note)


# ---------- s04: doublet fraction per cell type, both arms ----------

class S04CellType(Stage):
    name, bar = "s04_celltype", ("report-only: doublet fraction per cell type, both arms, on "
                                 "predicted.celltype.l1, with a predicted.celltype.l2 table beside it")
    asks = "Which cell types carry the most Scrublet-called doublets?"

    @staticmethod
    def _by_column(df, column):
        rows = {}
        for celltype, g in df.groupby(column):
            n = len(g)
            nd = int(g["predicted_doublet"].sum())
            rows[celltype] = {"n_cells": int(n), "n_doublets": nd,
                              "frac": round(nd / n, 4) if n else 0.0}
        return rows

    def run(self) -> StageResult:
        import pandas as pd

        d = SILO
        metrics, metrics_l2 = {}, {}
        for arm in ("headline", "sensitivity"):
            df = pd.read_csv(d / "scrublet" / f"{arm}.tsv", sep="\t")
            metrics[arm] = self._by_column(df, "celltype")
            metrics_l2[arm] = self._by_column(df, "celltype_l2")

        out = d / "per_celltype_doublets.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)

        note = "; ".join(f"{arm}: " + ", ".join(f"{ct} {v['frac']*100:.1f}%"
                          for ct, v in sorted(rows.items())) for arm, rows in metrics.items())

        out_l2 = d / "per_celltype_doublets_l2.json"
        _fresh(out_l2)
        out_l2.write_text(json.dumps(metrics_l2, indent=2))
        _made(out_l2)

        return StageResult(metrics={"l1": metrics, "l2": metrics_l2}, passed=None, note=note)


# ---------- s05: the before/after UMAP panel ----------

class S05Panel(Stage):
    name, bar = "s05_panel", "report-only: Supplementary Fig. 8, the before/after UMAP for both arms"
    asks = "Draw Supplementary Fig. 8."

    def run(self) -> StageResult:
        import anndata as ad
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        d = SILO
        before = ad.read_h5ad(d / "before.h5ad")
        colors = _celltype_colors(before.obs["celltype"].astype(str))
        arms = ("headline", "sensitivity")

        prior = {"headline": "prior 0.05", "sensitivity": "prior 0.10"}
        fig, axes = plt.subplots(len(arms), 3, figsize=(7.5, 2.7 * len(arms)))
        axes = np.atleast_2d(axes)
        for r, arm in enumerate(arms):
            calls = pd.read_csv(d / "scrublet" / f"{arm}.tsv", sep="\t").set_index("barcode")
            pred = calls.loc[before.obs_names, "predicted_doublet"].to_numpy().astype(bool)
            after = ad.read_h5ad(d / f"after_{arm}.h5ad")
            u = before.obsm["X_umap"]
            lab = prior[arm]

            ax = axes[r, 0]
            for ct, color in colors.items():
                sel = (before.obs["celltype"] == ct).to_numpy()
                ax.scatter(u[sel, 0], u[sel, 1], s=2, c=color, label=f"{ct} (n={int(sel.sum())})",
                           rasterized=True)
            ax.set_title(f"{lab}, all cells", fontsize=8)

            ax = axes[r, 1]
            ax.scatter(u[~pred, 0], u[~pred, 1], s=2, c="#d9d9d9", label="singlet", rasterized=True)
            ax.scatter(u[pred, 0], u[pred, 1], s=6, c="#000000",
                       label=f"Scrublet doublet (n={int(pred.sum())})", rasterized=True)
            ax.set_title(f"{lab}, Scrublet-called doublets", fontsize=8)

            ax = axes[r, 2]
            ua = after.obsm["X_umap"]
            for ct, color in colors.items():
                sel = (after.obs["celltype"] == ct).to_numpy()
                ax.scatter(ua[sel, 0], ua[sel, 1], s=2, c=color, label=f"{ct} (n={int(sel.sum())})",
                           rasterized=True)
            ax.set_title(f"{lab}, doublets removed", fontsize=8)

            for c in range(3):
                axes[r, c].set_xticks([])
                axes[r, c].set_yticks([])
                axes[r, c].legend(fontsize=4.5, markerscale=2, loc="best", frameon=False,
                                  handletextpad=0.3, labelspacing=0.25)
                axes[r, c].set_xlabel("UMAP 1", fontsize=8)
                axes[r, c].set_ylabel("UMAP 2", fontsize=8)

        fig.tight_layout()
        png, pdf = d / "Supp.Fig8.png", d / "Supp.Fig8.pdf"
        _fresh(png, pdf)
        fig.savefig(png, dpi=300)
        fig.savefig(pdf)
        plt.close(fig)
        _made(png, pdf)
        return StageResult(metrics={"panel": str(pdf)}, passed=None, note="Supplementary Fig. 8 written")


# ---------- s06: Leiden seed baseline for reading the before/after ARI ----------

class S06SeedBaseline(Stage):
    name, bar = "s06_seed_baseline", ("report-only: leiden re-run on the untouched before-graph at "
                                      "seed_baseline_seeds; ARI of each against seed 0")
    asks = ("How much does the clustering move when nothing but the random seed changes, on the "
            "same cells and the same neighbour graph?")

    def run(self) -> StageResult:
        import anndata as ad
        import scanpy as sc
        from sklearn.metrics import adjusted_rand_score

        d = SILO
        before = ad.read_h5ad(d / "before.h5ad")
        base = before.obs["leiden"].astype(str).to_numpy()
        b = BARS
        runs = {}
        for seed in b["seed_baseline_seeds"]:
            key = f"leiden_seed{seed}"
            sc.tl.leiden(before, resolution=b["cluster_leiden_resolution"], flavor="igraph",
                         n_iterations=2, directed=False, random_state=int(seed), key_added=key)
            lab = before.obs[key].astype(str).to_numpy()
            runs[str(seed)] = {"n_clusters": int(len(set(lab))),
                               "ari_vs_seed0": round(float(adjusted_rand_score(base, lab)), 4)}
        aris = [v["ari_vs_seed0"] for v in runs.values()]
        metrics = {"cells": int(before.n_obs), "n_clusters_seed0": int(len(set(base))),
                   "seeds": runs, "ari_min": min(aris), "ari_max": max(aris)}
        out = d / "seed_baseline.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        note = ", ".join(f"seed {k}: {v['n_clusters']} clusters, ARI {v['ari_vs_seed0']:.3f}"
                         for k, v in runs.items())
        return StageResult(metrics=metrics, passed=None, note=note)


STAGES = [S01Load(), S02Scrublet(), S03Cluster(), S04CellType(), S05Panel(), S06SeedBaseline()]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from-stage", metavar="NAME", help="resume at this stage")
    args = ap.parse_args()
    names = [s.name for s in STAGES]
    stages = STAGES[names.index(args.from_stage):] if args.from_stage else STAGES
    SILO.mkdir(parents=True, exist_ok=True)
    metrics_path = SILO / "stage_metrics.json"
    metrics = json.loads(metrics_path.read_text()) if metrics_path.exists() else {}
    for stage in stages:
        print(f"[{stage.name}] {stage.asks}\n[{stage.name}] bar: {stage.bar}")
        t0 = perf_counter()
        result = stage.run()
        status = "ok" if result.passed is None else ("pass" if result.passed else "FAIL")
        print(f"[{stage.name}] {status} ({perf_counter() - t0:.1f}s) {result.note}".rstrip())
        metrics[stage.name] = result.metrics
        metrics_path.write_text(json.dumps(metrics, indent=2, default=str))
        if result.passed is False:
            sys.exit(f"{stage.name} failed its check; stopping")


if __name__ == "__main__":
    main()
