"""K562 + LCL cell-hashing experiment (Methods 1.6, Supplementary Fig. 7):
doublet rates from the hashtags and cluster structure before and after multiplet removal.

Stages: load the hashing libraries, cluster, call multiplets with hashsolo, measure where the
multiplets sit, recluster without them, draw Supplementary Fig. 7, summarize the stability.
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
from time import perf_counter


FEED = Path(os.path.expandvars("$BENCHDROP_FEED"))
SILO = Path(os.path.expandvars("$BENCHDROP_RESULTS/hashing"))

HASHING = FEED / "hashing"      # <loading>/RNA/filtered_matrix/<sensitivity> and <loading>/RNA/SNT/snt_raw_matrix
CONDS = ("5k", "7k", "10k")
HASHED = ("7k", "10k")          # 5k is the K562-only design control and carries no SNT channel
SENSITIVITY = "sensitivity_3"   # PIPseeker cell-calling sensitivity level

# Cell-line programs used to call the line of each cluster.
K562_MK = ["HBB", "HBA1", "HBA2", "HBG1", "HBG2", "GYPA", "GATA1", "ALAS2", "KLF1", "CA1", "AHSP"]
LCL_MK = ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "IGKC", "IGHM", "CD19", "BANK1", "HLA-DRB1"]

# Clustering recipe.
N_TOP_GENES = 2000
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RES = 0.3
RANDOM_STATE = 0

BARS = {
    # ---- s01: each arm loads with at least this many called cells.
    "load_min_cells": 1_000,

    # ---- s05: a cluster counts as MAJOR at >= 1% of cells. A cell's LINE CALL is the identity
    # (K562 or LCL) of the major cluster it sits in; concordance is the fraction of retained cells
    # whose line call survives re-clustering without the multiplets.
    "umap_major_frac": 0.01,
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




# ---------- shared helpers ----------

def _fresh(*paths: Path):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _made(*paths: Path):
    for p in paths:
        if not Path(p).exists() or Path(p).stat().st_size == 0:
            raise RuntimeError(f"stage exited but produced no {p}")


def _matrix_dir(cond: str) -> Path:
    return HASHING / cond / "RNA" / "filtered_matrix" / SENSITIVITY


def _snt_dir(cond: str) -> Path:
    return HASHING / cond / "RNA" / "SNT" / "snt_raw_matrix"


def _cluster(cond: str, adata=None):
    """normalize -> log1p -> HVG -> PCA -> neighbours -> leiden -> UMAP.

    `adata` lets s05 re-run the identical recipe on a cell SUBSET from raw counts, which is the
    honest form of the stability test: re-derive the structure, never sub-set a fitted embedding.
    """
    import scanpy as sc

    if adata is None:
        adata = sc.read_10x_mtx(_matrix_dir(cond))
        adata.var_names_make_unique()
        adata.layers["counts"] = adata.X.copy()
    else:
        adata = adata.copy()
        adata.X = adata.layers["counts"].copy()
        # The h5ad carries uns['log1p'] from the BEFORE run; scanpy reads that as "already
        # log-transformed" and warns spuriously. X was just restored from raw counts (verified
        # integral, max 5685), so clear the stale key rather than leave a misleading warning.
        adata.uns.pop("log1p", None)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    sc.pp.pca(adata, n_comps=N_PCS, mask_var="highly_variable", random_state=RANDOM_STATE)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2,
                 directed=False, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    sc.tl.score_genes(adata, [g for g in K562_MK if g in adata.var_names], score_name="k562_sig")
    sc.tl.score_genes(adata, [g for g in LCL_MK if g in adata.var_names], score_name="lcl_sig")
    return adata


def _label_clusters(adata, major_frac: float):
    """Per-cluster: (line call, n cells, mean K562 sig, mean LCL sig, is_major)."""
    n = adata.n_obs
    out = {}
    for cl in adata.obs["leiden"].cat.categories:
        m = (adata.obs["leiden"] == cl).to_numpy()
        ks = float(adata.obs.loc[m, "k562_sig"].mean())
        ls = float(adata.obs.loc[m, "lcl_sig"].mean())
        out[str(cl)] = {"line": "K562" if ks > ls else "LCL", "n": int(m.sum()),
                        "k562_sig": round(ks, 4), "lcl_sig": round(ls, 4),
                        "major": bool(m.sum() >= major_frac * n)}
    return out


def _add_snt(adata, cond: str):
    """Join SNT tag counts to the RNA cells by barcode string."""
    from scipy.io import mmread

    d = _snt_dir(cond)
    feats = [l.split("\t")[0] for l in gzip.open(d / "features.tsv.gz", "rt")]   # ['lcl','k562']
    bcs = [l.strip() for l in gzip.open(d / "barcodes.tsv.gz", "rt")]
    M = mmread(str(d / "matrix.mtx.gz")).tocsr().toarray()                        # tags x barcodes
    bidx = {b: j for j, b in enumerate(bcs)}
    lr, kr = feats.index("lcl"), feats.index("k562")
    lcl = np.zeros(adata.n_obs)
    k = np.zeros(adata.n_obs)
    for i, b in enumerate(adata.obs_names):
        j = bidx.get(b)
        if j is not None:
            lcl[i], k[i] = M[lr, j], M[kr, j]
    adata.obs["lcl"], adata.obs["k562"] = lcl, k
    return adata, int((lcl + k > 0).sum())


# ---------- s01: the three loading arms are on disk ----------

class S01Load(Stage):
    name, bar = "s01_load", ("all three loading arms load with >= load_min_cells cells; both "
                             "hashed arms carry an SNT channel")
    asks = "Are the three loading arms of the hashing experiment on disk?"

    def run(self) -> StageResult:
        SILO.mkdir(parents=True, exist_ok=True)
        out = SILO / "arm_dims.json"
        _fresh(out)

        dims, problems = {}, []
        for cond in CONDS:
            md = _matrix_dir(cond)
            if not (md / "matrix.mtx.gz").exists():
                problems.append(f"{cond}: no matrix at {md}")
                continue
            with gzip.open(md / "matrix.mtx.gz", "rt") as fh:
                fh.readline()
                fh.readline()
                g, c, nnz = (int(x) for x in fh.readline().split())
            dims[cond] = {"genes": g, "cells": c, "nnz": nnz,
                          "snt": (_snt_dir(cond) / "matrix.mtx.gz").exists()}
            if c < BARS["load_min_cells"]:
                problems.append(f"{cond}: {c} cells < {BARS['load_min_cells']}")
            if cond in HASHED and not dims[cond]["snt"]:
                problems.append(f"{cond}: hashed arm with no SNT channel")

        out.write_text(json.dumps(dims, indent=2))
        _made(out)
        ok = not problems and len(dims) == len(CONDS)
        return StageResult(metrics=dims, passed=ok,
                           note="; ".join(problems) if problems
                                else " ".join(f"{k}={v['cells']}c" for k, v in dims.items()))


# ---------- s02: cluster each arm, call the line per cluster ----------

class S02Cluster(Stage):
    name, bar = "s02_cluster", "report-only: cluster each arm and call the line of every cluster"
    asks = "How many clusters does each arm form, and which line is each?"

    def run(self) -> StageResult:
        d = SILO / "clusters"
        d.mkdir(parents=True, exist_ok=True)

        metrics = {}
        for cond in CONDS:
            h5 = d / f"{cond}.h5ad"
            _fresh(h5)
            a = _cluster(cond)
            lab = _label_clusters(a, BARS["umap_major_frac"])
            a.obs["line"] = [lab[str(c)]["line"] for c in a.obs["leiden"]]
            a.write_h5ad(h5)
            _made(h5)

            k_cells = sum(v["n"] for v in lab.values() if v["line"] == "K562")
            l_cells = sum(v["n"] for v in lab.values() if v["line"] == "LCL")
            metrics[cond] = {"cells": int(a.n_obs), "n_clusters": len(lab),
                             "rna_k562": k_cells, "rna_lcl": l_cells, "clusters": lab}

        out = SILO / "cluster_summary.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        return StageResult(metrics={k: {x: y for x, y in v.items() if x != "clusters"}
                                    for k, v in metrics.items()},
                           passed=None,
                           note="; ".join(f"{k}: {v['n_clusters']} clusters, K562 {v['rna_k562']:,} / "
                                          f"LCL {v['rna_lcl']:,} cells" for k, v in metrics.items()))


# ---------- s03: hashsolo doublet rates ----------

class S03Hashsolo(Stage):
    name, bar = "s03_hashsolo", ("report-only: cross-line doublet rate from the hashtags and the "
                                 "collision-corrected total rate, both hashed arms")
    asks = "What doublet rate does the hashing experiment measure at each loading?"

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        from scanpy.external.pp import hashsolo

        d = SILO / "hashsolo"
        d.mkdir(parents=True, exist_ok=True)

        metrics = {}
        for cond in HASHED:
            tsv = d / f"{cond}.tsv"
            _fresh(tsv)
            a = ad.read_h5ad(SILO / "clusters" / f"{cond}.h5ad")
            a, n_tagged = _add_snt(a, cond)
            hashsolo(a, ["lcl", "k562"], number_of_noise_barcodes=1)   # noise=1 fixes the 2-tag case
            cl = a.obs["Classification"]
            vc = cl.value_counts()
            n = a.n_obs
            n_sl, n_sk = int(vc.get("lcl", 0)), int(vc.get("k562", 0))
            n_doub, n_neg = int(vc.get("Doublet", 0)), int(vc.get("Negative", 0))
            cross = n_doub / n
            p_k = n_sk / max(n_sk + n_sl, 1)
            total = cross / (2 * p_k * (1 - p_k)) if 0 < p_k < 1 else float("nan")

            # Validation: a hashtag doublet should carry BOTH RNA programs.
            def sig(mask, s):
                return round(float(a.obs.loc[mask, s].mean()), 4) if mask.any() else float("nan")
            dm, lm, km = (cl == "Doublet"), (cl == "lcl"), (cl == "k562")

            pd.DataFrame({"barcode": a.obs_names, "leiden": a.obs["leiden"].astype(str),
                          "line": a.obs["line"], "lcl": a.obs["lcl"], "k562": a.obs["k562"],
                          "classification": cl.astype(str),
                          "umap1": a.obsm["X_umap"][:, 0], "umap2": a.obsm["X_umap"][:, 1],
                          }).to_csv(tsv, sep="\t", index=False)
            _made(tsv)

            cp, tp = round(100 * cross, 1), round(100 * total, 1)
            metrics[cond] = {"cells": n, "tagged": n_tagged,
                             "sing_lcl": n_sl, "sing_k562": n_sk,
                             "doublet": n_doub, "negative": n_neg,
                             "cross_pct": cp, "total_pct": tp,
                             "doublet_k562_sig": sig(dm, "k562_sig"),
                             "doublet_lcl_sig": sig(dm, "lcl_sig"),
                             "k562singlet_k562_sig": sig(km, "k562_sig"),
                             "lclsinglet_lcl_sig": sig(lm, "lcl_sig")}

        out = SILO / "hashsolo_rates.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        return StageResult(metrics=metrics, passed=None,
                           note="; ".join(f"{k}: cross-line {v['cross_pct']}%, corrected total "
                                          f"{v['total_pct']}%" for k, v in metrics.items()))


# ---------- s04: where do the multiplet calls actually sit? (descriptive) ----------

class S04Enrichment(Stage):
    name, bar = "s04_enrichment", "report-only: multiplet load per cluster, both hashed arms"
    asks = "Do hashtag multiplet calls concentrate anywhere in the transcriptome, or scatter?"

    def run(self) -> StageResult:
        import pandas as pd

        out = SILO / "multiplet_enrichment.json"
        _fresh(out)
        metrics = {}
        for cond in HASHED:
            hs = pd.read_csv(SILO / "hashsolo" / f"{cond}.tsv", sep="\t")
            rows = {}
            for cl, g in hs.groupby("leiden"):
                rows[str(cl)] = {"n": int(len(g)), "line": str(g["line"].iloc[0]),
                                 "doublets": int((g["classification"] == "Doublet").sum()),
                                 "pct_doublet": round(100 * float(
                                     (g["classification"] == "Doublet").mean()), 1)}
            top = max(rows.values(), key=lambda r: r["pct_doublet"])
            rest = max((r["pct_doublet"] for r in rows.values() if r is not top), default=0.0)
            metrics[cond] = {"clusters": rows, "max_pct_doublet": top["pct_doublet"],
                             "next_pct_doublet": rest, "enriched_cluster_n": top["n"]}
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        return StageResult(
            metrics={k: {x: y for x, y in v.items() if x != "clusters"} for k, v in metrics.items()},
            passed=None,
            note="; ".join(f"{k}: top cluster {v['max_pct_doublet']}% multiplet vs "
                           f"{v['next_pct_doublet']}% next" for k, v in metrics.items()))


# ---------- s05: re-derive the embedding without multiplets (descriptive) ----------

class S05Recluster(Stage):
    name, bar = "s05_recluster", ("report-only: re-derive the embedding from raw counts on the "
                                  "multiplet-free cells and measure ARI, major-cluster count, "
                                  "and line-call concordance")
    asks = "What happens to the cluster structure when the hashtag-called multiplets are removed?"

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        from sklearn.metrics import adjusted_rand_score

        d = SILO / "stability"
        d.mkdir(parents=True, exist_ok=True)
        metrics = {}

        for cond in HASHED:
            tsv = d / f"{cond}_after.tsv"
            _fresh(tsv)
            before = ad.read_h5ad(SILO / "clusters" / f"{cond}.h5ad")
            hs = pd.read_csv(SILO / "hashsolo" / f"{cond}.tsv", sep="\t").set_index("barcode")
            before.obs["classification"] = hs.loc[before.obs_names, "classification"].to_numpy()

            # Remove the MULTIPLETS only. `Negative` is a cell the tag channel failed to call,
            # not a called multiplet, so dropping it would answer a question nobody asked.
            keep = (before.obs["classification"] != "Doublet").to_numpy()
            n_removed = int((~keep).sum())
            after = _cluster(cond, adata=before[keep].copy())

            b_lab = _label_clusters(before, BARS["umap_major_frac"])
            a_lab = _label_clusters(after, BARS["umap_major_frac"])
            b_major = [c for c, v in b_lab.items() if v["major"]]
            a_major = [c for c, v in a_lab.items() if v["major"]]

            ari = float(adjusted_rand_score(before.obs["leiden"][keep].astype(str),
                                            after.obs["leiden"].astype(str)))
            # Line call travels with the cell: its cluster's call before vs after.
            b_line = np.array([b_lab[str(c)]["line"] for c in before.obs["leiden"][keep]])
            a_line = np.array([a_lab[str(c)]["line"] for c in after.obs["leiden"]])
            conc = float((b_line == a_line).mean())
            d_major = abs(len(a_major) - len(b_major))

            pd.DataFrame({"barcode": after.obs_names, "leiden_after": after.obs["leiden"].astype(str),
                          "line_after": a_line, "leiden_before": before.obs["leiden"][keep].astype(str),
                          "line_before": b_line,
                          "umap1": after.obsm["X_umap"][:, 0], "umap2": after.obsm["X_umap"][:, 1],
                          }).to_csv(tsv, sep="\t", index=False)
            _made(tsv)

            metrics[cond] = {"cells_before": int(before.n_obs), "multiplets_removed": n_removed,
                             "cells_after": int(after.n_obs),
                             "clusters_before": len(b_lab), "clusters_after": len(a_lab),
                             "major_before": len(b_major), "major_after": len(a_major),
                             "ari": round(ari, 4), "line_concordance": round(conc, 4),
                             "major_lines_before": sorted({b_lab[c]["line"] for c in b_major}),
                             "major_lines_after": sorted({a_lab[c]["line"] for c in a_major})}

        out = SILO / "umap_stability.json"
        _fresh(out)
        out.write_text(json.dumps(metrics, indent=2))
        _made(out)
        return StageResult(metrics=metrics, passed=None,
                           note="; ".join(f"{k}: ARI {v['ari']:.3f}, conc {v['line_concordance']:.3f}, "
                                          f"major {v['major_before']}->{v['major_after']}"
                                          for k, v in metrics.items()))


# ---------- s06: Supplementary Fig. 7 ----------

class S06Panel(Stage):
    name, bar = "s06_panel", "report-only: the before/after panel renders for both hashed arms"
    asks = "Assemble the before/after UMAP panel (Supplementary Fig. 7)."

    def run(self) -> StageResult:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        colors = {"K562": "#1f77b4", "LCL": "#d62728"}
        loading = {"7k": "7,000", "10k": "10,000"}
        fig, axes = plt.subplots(len(HASHED), 3, figsize=(7.2, 2.4 * len(HASHED)))
        axes = np.atleast_2d(axes)
        for r, cond in enumerate(HASHED):
            bef = pd.read_csv(SILO / "hashsolo" / f"{cond}.tsv", sep="\t")
            aft = pd.read_csv(SILO / "stability" / f"{cond}_after.tsv", sep="\t")
            lab = loading.get(cond, cond)

            ax = axes[r, 0]
            for line, sub in bef.groupby("line"):
                ax.scatter(sub["umap1"], sub["umap2"], s=2, c=colors.get(line, "grey"),
                           label=f"{line} (n={len(sub)})", rasterized=True)
            ax.set_title(f"{lab} loading, all barcodes", fontsize=8)

            ax = axes[r, 1]
            dm = bef["classification"] == "Doublet"
            ax.scatter(bef.loc[~dm, "umap1"], bef.loc[~dm, "umap2"], s=2, c="#d9d9d9",
                       label="hashtag singlet / negative", rasterized=True)
            ax.scatter(bef.loc[dm, "umap1"], bef.loc[dm, "umap2"], s=6, c="#000000",
                       label=f"hashtag multiplet (n={int(dm.sum())})", rasterized=True)
            ax.set_title(f"{lab} loading, hashtag-called multiplets", fontsize=8)

            ax = axes[r, 2]
            for line, sub in aft.groupby("line_after"):
                ax.scatter(sub["umap1"], sub["umap2"], s=2, c=colors.get(line, "grey"),
                           label=f"{line} (n={len(sub)})", rasterized=True)
            ax.set_title(f"{lab} loading, multiplets removed", fontsize=8)

            for c in range(3):
                axes[r, c].set_xticks([])
                axes[r, c].set_yticks([])
                axes[r, c].legend(fontsize=6, markerscale=3, loc="best", frameon=False)
                axes[r, c].set_xlabel("UMAP 1", fontsize=8)
                axes[r, c].set_ylabel("UMAP 2", fontsize=8)

        fig.tight_layout()
        png, pdf = SILO / "Supp.Fig7.png", SILO / "Supp.Fig7.pdf"
        _fresh(png, pdf)
        fig.savefig(png, dpi=300)
        fig.savefig(pdf)
        plt.close(fig)
        _made(png, pdf)
        return StageResult(metrics={"panel": str(pdf)}, passed=None, note="Supplementary Fig. 7 written")


# ---------- s07: the stability summary the text quotes ----------

class S07Stability(Stage):
    name, bar = "s07_stability", ("report-only: after removing hashtag-called multiplets and "
                                  "re-deriving the embedding from raw counts, the major-cluster "
                                  "count, the ARI, and the line-call concordance on both hashed arms")
    asks = ("How much does removing hashtag-called multiplets change the hashing experiment's "
            "cluster structure?")

    def run(self) -> StageResult:
        metrics = json.loads((SILO / "umap_stability.json").read_text())
        return StageResult(metrics=metrics, passed=None,
                           note="; ".join(f"{k}: major clusters {v['major_before']}->{v['major_after']}, "
                                          f"ARI {v['ari']:.3f}, line concordance "
                                          f"{v['line_concordance']:.1%}" for k, v in metrics.items()))


STAGES = [S01Load(), S02Cluster(), S03Hashsolo(), S04Enrichment(), S05Recluster(),
          S06Panel(), S07Stability()]


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
