"""Isoform usage along a within-lineage pseudotime axis (Fig. 6, Supplementary Table 3, Methods 2.8).

Stages: align, read geometry, transcript matrix, subtypes, pseudotime, ordering check, isoform
usage along the axis, depth control, Fig. 6.
"""
from __future__ import annotations
import argparse
import gzip
import json
import numpy as np
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter


FEED = Path(os.path.expandvars("$BENCHDROP_FEED"))
SILO = Path(os.path.expandvars("$BENCHDROP_RESULTS/pseudotime"))

RAW_FASTQ = FEED / "pbmc/raw/pbmc.top30Mreads.fastq.gz"
MATRIX_DIR = FEED / "pbmc/count_lr"          # the PBMC EM transcript matrix
CELLTYPES = FEED / "pbmc/pbmc_celltypes.tsv"   # barcode -> lineage (b_cell, t_cell, nk, mono), shipped input
TRANSCRIPTOME = FEED / "genome/transcriptome.fa"          # GENCODE v32, 199,138 bare-ENST records
TX_MAP = FEED / "annotation_maps/gencode_v32.tx.tsv"      # tx -> gene/name/biotype
GENE_MAP = FEED / "annotation_maps/gencode_v32.gene.tsv"  # gene -> name/biotype

READS_BIN = Path(os.path.expandvars("$BENCHDROP_READS_BIN"))
PROFILER = Path(__file__).resolve().parent / "misc" / "phasing_coverage_profile.py"

# Expected matrix header; a different quantification under the same path stops the run here
# instead of silently shifting every number.
EXPECT_MTX = (2_109_262, 199_138, 41_738_415)

# Exemplar genes checked for discriminating sequence inside the covered 3' window.
EXEMPLARS = ("IGHM", "CD8A", "NKG7", "IL7R")

# Canonical within-lineage orderings:
#   T: naive -> memory;  B: naive -> memory -> plasmablast.
CANONICAL_PAIRS = (("t_naive", "t_memory"), ("b_naive", "b_memory"), ("b_memory", "plasmablast"))

# Marker sets naming the subtypes (the labels themselves come from the reference transfer).
MARKERS = {
    "t_naive": ["CCR7", "SELL", "TCF7", "LEF1"],
    "t_memory": ["S100A4", "GZMK", "CCL5", "ITGB1"],
    "b_naive": ["TCL1A", "IGHD", "FCER2", "IL4R"],
    "b_memory": ["CD27", "TNFRSF13B", "AIM2"],
    "plasmablast": ["MZB1", "XBP1", "PRDM1", "JCHAIN", "CD38"],
}

# Azimuth PBMC reference (Hao et al. 2021 CITE-seq, 161,764 cells; Human PBMC Reference 2.10).
AZIMUTH_REF = FEED / "azimuth_pbmc/pbmc_multimodal.h5seurat"
SHORT_READ = FEED / "pbmc/count_sr/Gene/filtered"   # matched SR, same barcodes

BARS = {
    # s04 label transfer: k nearest reference cells per query cell.
    "az_k": 50,

    # s01: enough primary alignments for stable percentiles.
    "align_min_primary": 100_000,
    # s02: the profile is valid only if reads are 3'-anchored (median distance from the 3' end
    # within this many bp) and there are enough of them.
    "geometry_max_median_tpdist": 200,
    "geometry_min_n": 50_000,
    # s03: labelled-cell recovery + annotation coverage.
    "matrix_min_typed_cells": 16_000,
    "matrix_min_txmap_frac": 0.99,
    # s04/s05: an ordering pair is scoreable only with both ends populated.
    "subtype_min_cells": 30,
    # s05: a degenerate axis (constant, or non-finite for stragglers) cannot be ordered.
    "pseudotime_min_finite_frac": 0.99,
    # s06: within-lineage label permutations for the ordering AUC.
    "ordering_n_perm": 10_000,
    # s07:
    #   expression-invariant  = gene-level fold change between top and bottom pseudotime
    #                           quintile <= 1.5x
    #   coverage criterion    = >=10% unique 31-mer fraction, computed on the isoform's
    #                           terminal window (last W bp, W = the s02 median span)
    #   usage shift           = per-isoform binned usage vs pseudotime, Spearman, BH FDR<0.05,
    #                           with a >=10-percentage-point usage swing
    "usage_min_gene_umi": 1_000,     # gene EM-UMIs across typed cells; 20 bins keep >=50 each
    "usage_min_isoform_frac": 0.10,  # an isoform participates at >=10% pseudobulk usage
    "usage_kmer": 31,
    "usage_min_unique_frac": 0.10,   # the coverage criterion, on the terminal window
    "usage_n_bins": 20,
    "usage_min_bin_umi": 25,         # a bin with fewer gene UMIs is masked
    "usage_min_used_bins": 15,
    "usage_max_fdr": 0.05,
    "usage_min_swing": 0.10,
    "usage_max_gene_fc": 1.5,        # the expression-invariance bound
    "threads": 8,
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


def _run(cmd, log_path: Path, **kw):
    """Subprocess with stderr kept (never muted) in a log beside the results."""
    with open(log_path, "ab") as log:
        log.write((" ".join(str(c) for c in cmd) + "\n").encode())
        log.flush()
        subprocess.run([str(c) for c in cmd], stderr=log, check=True, **kw)


def _read_fasta(path: Path) -> dict:
    """transcriptome.fa -> {bare ENST: sequence}. 404 MB, fits comfortably."""
    seqs, cur, buf = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = ln[1:].split()[0].strip(), []
            else:
                buf.append(ln.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def _tx_tables():
    """Annotation maps -> (tx->gene, gene->[tx], gene->name, name->gene, gene->biotype)."""
    tx2g, g2tx, g2name, g2bio = {}, defaultdict(list), {}, {}
    with open(TX_MAP) as fh:
        header = fh.readline()
        assert header.startswith("tx_id")
        for ln in fh:
            tx, g, name, gbio = ln.rstrip("\n").split("\t")[:4]
            tx2g[tx] = g
            g2tx[g].append(tx)
            g2name[g] = name
            g2bio[g] = gbio
    name2g = {}
    for g, name in g2name.items():
        name2g.setdefault(name, g)
    return tx2g, g2tx, g2name, name2g, g2bio


def _window_unique_frac(tx: str, siblings: list, seqs: dict, window: int, k: int) -> float:
    """Fraction of the last-`window`-bp k-mers of `tx` absent from every full-length sibling.

    The staircase identifiability rule (31-mer uniqueness vs the gene's other isoforms),
    restricted to the terminal window a typical 3'-anchored read actually spans: one number
    that enforces the >=10% floor ON the covered-3'-window sequence.
    """
    s = seqs.get(tx, "")
    if len(s) < k:
        return 0.0
    tail = s[-window:] if len(s) > window else s
    kmers = {tail[i:i + k] for i in range(len(tail) - k + 1)}
    if not kmers:
        return 0.0
    sib = set()
    for o in siblings:
        so = seqs.get(o, "")
        sib.update(so[i:i + k] for i in range(len(so) - k + 1))
    return len(kmers - sib) / len(kmers)


# ---------- s01: subsample + transcriptome alignment ----------

class S01Align(Stage):
    name, bar = "s01_align", "primary mapped alignments >= align_min_primary (of 500k reads)"
    asks = "Do enough PBMC reads transcriptome-align to profile the run's read geometry?"

    def run(self) -> StageResult:
        d = SILO / "geom"
        d.mkdir(parents=True, exist_ok=True)
        sub, bam = d / "pbmc_sub500k.fastq.gz", d / "pbmc_geom.bam"
        _fresh(sub, bam)
        # First 500k reads of the top-30M subset; the profiler caps at 200k mapped primaries.
        with open(d / "subsample.log", "wb") as log:
            p1 = subprocess.Popen([READS_BIN / "unpigz", "-c", RAW_FASTQ],
                                  stdout=subprocess.PIPE, stderr=log)
            p2 = subprocess.Popen(["head", "-n", "2000000"], stdin=p1.stdout,
                                  stdout=subprocess.PIPE, stderr=log)
            p1.stdout.close()
            with open(sub, "wb") as out:
                p3 = subprocess.Popen([READS_BIN / "pigz", "-p", "4"], stdin=p2.stdout,
                                      stdout=out, stderr=log)
                p2.stdout.close()
                p3.wait()
            p2.wait()
            p1.wait()   # SIGPIPE from head is expected; p3's status is the one that matters
            if p3.returncode:
                raise RuntimeError("subsample recompression failed")
        _made(sub)
        # Same minimap2 flags as the staircase.
        aln = subprocess.Popen(
            [str(READS_BIN / "minimap2"), "-ax", "map-ont", "--for-only", "-N", "200",
             "-p", "0.9", "-t", str(BARS["threads"]), str(TRANSCRIPTOME), str(sub)],
            stdout=subprocess.PIPE, stderr=open(d / "minimap2.log", "wb"))
        _run([READS_BIN / "samtools", "view", "-b", "-o", bam, "-"],
             log_path=d / "samtools.log", stdin=aln.stdout)
        aln.stdout.close()
        if aln.wait():
            raise RuntimeError("minimap2 failed; see geom/minimap2.log")
        _made(bam)
        flag = subprocess.run([str(READS_BIN / "samtools"), "flagstat", str(bam)],
                              capture_output=True, text=True, check=True).stdout
        (d / "flagstat.txt").write_text(flag)
        primary = int(re.search(r"(\d+) \+ \d+ primary mapped", flag).group(1))
        ok = primary >= BARS["align_min_primary"]
        return StageResult(metrics={"primary_mapped": primary},
                           passed=ok, note=f"{primary:,} primary mapped")


# ---------- s02: geometry profile + exemplar 3'-window check ----------

class S02Geometry(Stage):
    name, bar = "s02_geometry", ("median tp_dist <= geometry_max_median_tpdist and "
                                 "n >= geometry_min_n (3'-anchoring check)")
    asks = ("What 3' window do PBMC reads cover, and does each exemplar's discriminating "
            "sequence sit inside it?")

    def run(self) -> StageResult:
        d = SILO / "geom"
        bam, prof = d / "pbmc_geom.bam", d / "geometry_profile.txt"
        _fresh(prof)
        # misc/phasing_coverage_profile.py, run under the reads environment (pysam lives there);
        # its stdout is kept beside the parsed JSON.
        out = subprocess.run([str(READS_BIN / "python"), str(PROFILER), "--bam", str(bam)],
                             capture_output=True, text=True, check=True).stdout
        prof.write_text(out)
        _made(prof)
        # Parse the profiler's fixed-format lines: "  span_bp  n=..  p10=.. p25=..  MED=.." .
        stats = {}
        for ln in out.splitlines():
            m = re.match(r"\s+(\w+)\s+n=\s*([\d,]+)\s+p10=(\S+) p25=(\S+)\s+MED=(\S+)\s+p75=(\S+) p90=(\S+)", ln)
            if m:
                stats[m.group(1)] = dict(n=int(m.group(2).replace(",", "")),
                                         p10=float(m.group(3)), p25=float(m.group(4)),
                                         med=float(m.group(5)), p75=float(m.group(6)),
                                         p90=float(m.group(7)))
        span, tpd = stats["span_bp"], stats["tp_dist_3end"]
        window = int(span["med"])
        ok = tpd["med"] <= BARS["geometry_max_median_tpdist"] and span["n"] >= BARS["geometry_min_n"]

        # Exemplar check: windowed unique fraction per annotated isoform of each exemplar,
        # window = this run's median span.
        seqs = _read_fasta(TRANSCRIPTOME)
        _, g2tx, g2name, name2g, _ = _tx_tables()
        k = BARS["usage_kmer"]
        rows = []
        for name in EXEMPLARS:
            g = name2g.get(name)
            txs = g2tx.get(g, [])
            for t in txs:
                wu = _window_unique_frac(t, [o for o in txs if o != t], seqs, window, k)
                rows.append((name, g, t, len(seqs.get(t, "")), round(wu, 4)))
        chk = d / "exemplar_window_check.tsv"
        with open(chk, "w") as fh:
            fh.write("gene\tgene_id\ttx_id\ttx_len\twindow_unique_frac\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        _made(chk)
        best = {name: max((r[4] for r in rows if r[0] == name), default=0.0)
                for name in EXEMPLARS}
        (d / "geometry.json").write_text(json.dumps(
            {"window_bp": window, "span": span, "tp_dist": tpd,
             "exemplar_best_window_unique": best}, indent=1))
        return StageResult(
            metrics={"window_bp": window, "median_span": span["med"], "p25_span": span["p25"],
                     "median_tpdist": tpd["med"], "n": span["n"],
                     "exemplar_best_window_unique": best},
            passed=ok,
            note=f"window={window}bp (median span), tp_dist med={tpd['med']:.0f}, n={span['n']:,}")


# ---------- s03: matrix load, typed-cell subset, gene aggregation ----------

class S03Matrix(Stage):
    name, bar = "s03_matrix", ("mtx header == expected dims; labelled cells matched >= "
                               "matrix_min_typed_cells; tx->gene coverage >= matrix_min_txmap_frac")
    asks = "Is this the expected PBMC EM matrix, and do the labelled cells load cleanly?"

    def run(self) -> StageResult:
        import anndata as ad
        import scipy.io as sio
        import scipy.sparse as sp

        d = SILO
        d.mkdir(parents=True, exist_ok=True)
        out_tx, out_gene = d / "pbmc_typed_tx.h5ad", d / "pbmc_typed_gene.h5ad"
        _fresh(out_tx, out_gene)

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
        ct = {}
        with open(CELLTYPES) as fh:
            fh.readline()
            for ln in fh:
                b, t = ln.rstrip("\n").split("\t")
                ct[b] = t
        rows = [i for i, b in enumerate(bcs) if b in ct]
        matched = len(rows)

        m = sio.mmread(MATRIX_DIR / "matrix.mtx.gz").tocsr()      # barcodes x tx
        sub = m[rows, :]
        obs_bc = [bcs[i] for i in rows]

        tx2g, _, g2name, _, g2bio = _tx_tables()
        mapped = sum(1 for t in feats if t in tx2g)
        frac = mapped / len(feats)

        genes = sorted({tx2g[t] for t in feats if t in tx2g})
        gidx = {g: j for j, g in enumerate(genes)}
        col = np.array([gidx.get(tx2g.get(t), -1) for t in feats])
        keep = col >= 0
        # tx->gene aggregation as a sparse indicator product: (cells x tx) @ (tx x genes).
        ind = sp.csr_matrix((np.ones(keep.sum()), (np.nonzero(keep)[0], col[keep])),
                            shape=(len(feats), len(genes)))
        gsub = (sub @ ind).tocsr()

        import pandas as pd
        obs = pd.DataFrame({"celltype": [ct[b] for b in obs_bc]}, index=obs_bc)
        atx = ad.AnnData(sub, obs=obs.copy(), var=pd.DataFrame(index=feats))
        agene = ad.AnnData(gsub, obs=obs.copy(), var=pd.DataFrame(
            {"gene_name": [g2name[g] for g in genes],
             "biotype": [g2bio[g] for g in genes]}, index=genes))
        atx.write_h5ad(out_tx)
        agene.write_h5ad(out_gene)
        _made(out_tx, out_gene)
        ok = matched >= BARS["matrix_min_typed_cells"] and frac >= BARS["matrix_min_txmap_frac"]
        return StageResult(metrics={"mtx_header": dims, "typed_matched": matched,
                                    "txmap_frac": round(frac, 4), "n_genes": len(genes)},
                           passed=ok, note=f"{matched:,} typed cells, txmap {frac:.3f}")


# ---------- s04: subtypes by transfer from the Azimuth PBMC reference ----------

class S04Subtypes(Stage):
    name, bar = "s04_subtypes", ("report-only: reference transfer; each canonical subtype "
                                 "counted (an ordering pair needs >= subtype_min_cells at each end)")
    asks = "What are these cells, according to the Azimuth PBMC reference?"

    # Labels come from the Azimuth PBMC reference, transferred by SPCA projection + kNN vote onto
    # the matched short-read matrix, then carried to the long-read cells by barcode. This is not
    # the Azimuth package; misc/azimuth_labels.py states what differs.

    def run(self) -> StageResult:
        import pandas as pd

        sys.path.insert(0, str(Path(__file__).resolve().parent / "misc"))
        import azimuth_labels as AZ

        ref = AZ.read_reference(AZIMUTH_REF)
        X, genes, bcs = AZ.load_query()
        l1, l2, votes = AZ.transfer(ref, X, genes, k=int(BARS["az_k"]))
        az = pd.DataFrame({"barcode": bcs, "az_l1": l1, "az_l2": l2, "vote": votes}
                          ).set_index("barcode")
        lin = pd.read_csv(CELLTYPES, sep="\t").set_index("barcode")
        j = lin.join(az, how="inner")
        j["subtype"] = j["az_l2"].map(AZ.L2_TO_SUBTYPE)
        j["lineage"] = j["subtype"].map(lambda s: "T" if isinstance(s, str) and s.startswith("t_")
                                        else ("B" if isinstance(s, str) else None))
        med_vote = float(j["vote"].median())

        out = SILO / "pbmc_subtypes.tsv"
        _fresh(out)
        j.to_csv(out, sep="\t")
        _made(out)

        counts = {s: int((j["subtype"] == s).sum()) for s in
                  ("t_naive", "t_memory", "b_naive", "b_memory", "plasmablast")}
        return StageResult(
            metrics={"median_vote": round(med_vote, 3), "subtype_counts": counts,
                     "n_labelled": int(len(j))},
            passed=None,
            note=(f"median vote {med_vote:.2f}; "
                  + ", ".join(f"{s}={counts[s]}" for s in counts)))


# ---------- s05: the pseudotime axis ----------

class S05Pseudotime(Stage):
    name, bar = "s05_pseudotime", ("per lineage: finite pseudotime for >= "
                                   "pseudotime_min_finite_frac of cells, non-degenerate "
                                   "(sd > 0), and the root cell inside that lineage")
    asks = ("Does a within-lineage diffusion-pseudotime axis order the T and B compartments "
            "from their naive poles?")

    # One axis per lineage: PBMCs are co-existing lineages, not one continuum, so a single
    # global axis would order composition rather than cell state. Lineage membership and the
    # naive pole come from the transferred labels.
    LINEAGES = {"T": ("T", "t_naive"), "B": ("B", "b_naive")}

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        import scanpy as sc

        base = ad.read_h5ad(SILO / "pbmc_typed_gene.h5ad")
        sub = pd.read_csv(SILO / "pbmc_subtypes.tsv", sep="\t", index_col=0)
        base.obs = base.obs.join(sub[["subtype", "lineage", "az_l2", "vote"]])
        frames, metrics = [], {}
        for lin, (lineage_key, naive_subtype) in self.LINEAGES.items():
            a = base[base.obs["lineage"] == lineage_key].copy()
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
            sc.pp.highly_variable_genes(a, n_top_genes=2000)
            a = a[:, a.var.highly_variable].copy()
            sc.pp.scale(a, max_value=10)
            sc.tl.pca(a, n_comps=50, svd_solver="arpack")
            sc.pp.neighbors(a, n_neighbors=15)
            sc.tl.umap(a)
            sc.tl.diffmap(a)
            # Root: the medoid, in diffusion space, of the cells the reference calls naive (a
            # medoid rather than an argmax so one outlier cell cannot set the root).
            cand = np.nonzero((a.obs["subtype"] == naive_subtype).values)[0]
            if len(cand) > 500:      # cap the pairwise block; 500 medoid candidates is ample
                cand = cand[np.linspace(0, len(cand) - 1, 500).astype(int)]
            dm = a.obsm["X_diffmap"][cand, 1:6]
            d = np.linalg.norm(dm[:, None, :] - dm[None, :, :], axis=-1).sum(1)
            a.uns["iroot"] = int(cand[int(np.argmin(d))])
            root_ct = str(a.obs["lineage"].iloc[a.uns["iroot"]])
            root_sub = str(a.obs["subtype"].iloc[a.uns["iroot"]])
            sc.tl.dpt(a)
            pt = a.obs["dpt_pseudotime"].values
            frames.append(pd.DataFrame(
                {"barcode": a.obs_names, "celltype": a.obs["celltype"],
                 "subtype": a.obs["subtype"], "lineage": lin, "pseudotime": pt,
                 "umap1": a.obsm["X_umap"][:, 0], "umap2": a.obsm["X_umap"][:, 1]}))
            metrics[lin] = {"n": int(a.n_obs), "finite_frac": round(float(np.isfinite(pt).mean()), 4),
                            "sd": round(float(np.nanstd(pt)), 4), "root_celltype": root_ct,
                            "root_subtype": root_sub}
        out = SILO / "pseudotime.tsv"
        _fresh(out)
        pd.concat(frames).to_csv(out, sep="\t", index=False)
        _made(out)
        ok = all(m["finite_frac"] >= BARS["pseudotime_min_finite_frac"] and m["sd"] > 0
                 and m["root_celltype"] == self.LINEAGES[lin][0]
                 for lin, m in metrics.items())
        return StageResult(metrics=metrics, passed=ok,
                           note="; ".join(f"{l}: n={m['n']}, sd={m['sd']}, root={m['root_subtype']}"
                                          for l, m in metrics.items()))


# ---------- s06: ordering validation (inversion rate + permutation null) ----------

class S06Ordering(Stage):
    name, bar = "s06_ordering", ("report-only: mean canonical-pair AUC against a within-lineage "
                                 "label-permutation null")
    asks = ("Does the axis respect the canonical immune orderings (naive->memory T; "
            "naive B->memory B->plasmablast)?")

    def run(self) -> StageResult:
        import pandas as pd

        df = pd.read_csv(SILO / "pseudotime.tsv", sep="\t")
        rng = np.random.default_rng(20260813)

        def auc(early, late):
            """P(pt_early < pt_late) via the Mann-Whitney U statistic."""
            from scipy.stats import rankdata
            x = np.concatenate([early, late])
            r = rankdata(x)
            n1, n2 = len(early), len(late)
            u = r[n1:].sum() - n2 * (n2 + 1) / 2          # rank-sum of the LATE group
            return u / (n1 * n2)

        # Each canonical pair is scored inside its own lineage axis.
        lineages = {"T": ["t_naive", "t_memory"], "B": ["b_naive", "b_memory", "plasmablast"]}
        lin_of = {s: lin for lin, subs in lineages.items() for s in subs}
        pairs = []
        for a_, b_ in CANONICAL_PAIRS:
            lin = lin_of[a_]
            d = df[df.lineage == lin]
            e = d.loc[d.subtype == a_, "pseudotime"].values
            l = d.loc[d.subtype == b_, "pseudotime"].values
            if min(len(e), len(l)) >= BARS["subtype_min_cells"]:
                pairs.append((a_, b_, e, l, lin))
        obs = {f"{a_}->{b_}": auc(e, l) for a_, b_, e, l, _ in pairs}
        mean_auc = float(np.mean(list(obs.values())))
        inversions = sum(1 for v in obs.values() if v < 0.5)

        # Null: shuffle subtype labels within each lineage's own axis and rescore.
        groups = {}
        for lin in lineages:
            d = df[(df.lineage == lin) & (df.subtype.isin(lineages[lin]))]
            groups[lin] = (d["pseudotime"].values, d["subtype"].values)
        nperm = BARS["ordering_n_perm"]
        null = np.empty(nperm)
        for i in range(nperm):
            vals = []
            for lin, (pt, labs) in groups.items():
                perm = rng.permutation(labs)
                for a_, b_, _e, _l, plin in pairs:
                    if plin == lin:
                        vals.append(auc(pt[perm == a_], pt[perm == b_]))
            null[i] = np.mean(vals)
        p = float((1 + (null >= mean_auc).sum()) / (1 + nperm))
        out = SILO / "ordering_validation.json"
        _fresh(out)
        out.write_text(json.dumps(
            {"per_pair_auc": {k: round(v, 4) for k, v in obs.items()},
             "mean_auc": round(mean_auc, 4), "inversions": inversions,
             "n_pairs": len(pairs),
             "inversion_rate": round(inversions / len(pairs), 4) if pairs else None,
             "perm_p": p, "n_perm": nperm}, indent=1))
        _made(out)
        return StageResult(metrics={"per_pair_auc": obs, "mean_auc": mean_auc,
                                    "inversions": inversions, "n_pairs": len(pairs), "p": p},
                           passed=None,
                           note=f"mean AUC {mean_auc:.3f}, {inversions}/{len(pairs)} inverted, p={p:.4g}")


# ---------- s07: expression-invariant isoform-usage shifts ----------

class S07Usage(Stage):
    name, bar = "s07_usage", ("report-only: genes with isoform-usage shift at BH FDR < "
                              "usage_max_fdr, swing >= usage_min_swing, gene-level quintile FC "
                              "<= usage_max_gene_fc, above the windowed unique-sequence criterion; "
                              "the qualifying rows are Supplementary Table 3")
    asks = ("Which genes switch isoform usage along the axis while their total expression "
            "stays flat, above the floor and inside the covered 3' window?")

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        from scipy.stats import spearmanr
        from statsmodels.stats.multitest import multipletests

        atx_all = ad.read_h5ad(SILO / "pbmc_typed_tx.h5ad")
        agene_all = ad.read_h5ad(SILO / "pbmc_typed_gene.h5ad")
        pt_all = pd.read_csv(SILO / "pseudotime.tsv", sep="\t", index_col=0)

        tx2g, g2tx, g2name, _, g2bio = _tx_tables()
        window = json.loads((SILO / "geom" / "geometry.json").read_text())["window_bp"]
        k = BARS["usage_kmer"]
        seqs = _read_fasta(TRANSCRIPTOME)
        n_bins = BARS["usage_n_bins"]

        # Each lineage is tested on its own axis, so a shift inside one lineage cannot be
        # cross-lineage composition.
        results, curves, screened = [], [], {}
        for lin in sorted(pt_all["lineage"].unique()):
            pt = pt_all[pt_all.lineage == lin]
            atx = atx_all[pt.index]
            order = np.argsort(pt["pseudotime"].values)
            X = atx.X.tocsc()[order, :]
            bins = np.array_split(np.arange(X.shape[0]), n_bins)
            feats = list(atx.var_names)
            fidx = {t: j for j, t in enumerate(feats)}

            tx_tot = np.asarray(X.sum(0)).ravel()
            gene_tot = defaultdict(float)
            for t, j in fidx.items():
                g = tx2g.get(t)
                if g:
                    gene_tot[g] += tx_tot[j]
            candidates = [g for g, tot in gene_tot.items()
                          if tot >= BARS["usage_min_gene_umi"]
                          and g2bio.get(g) == "protein_coding" and len(g2tx[g]) >= 2]

            # Expression-invariance is judged on the SAME cells the usage test uses.
            agene = agene_all[pt.index]
            G = agene.X.tocsc()[order, :]
            lib = np.asarray(G.sum(1)).ravel()
            lib[lib == 0] = 1
            quint = np.array_split(np.arange(G.shape[0]), 5)
            gidx = {g: j for j, g in enumerate(agene.var_names)}

            n_pass_iso, n_pass_floor = 0, 0
            for g in candidates:
                txs = [t for t in g2tx[g] if t in fidx]
                if len(txs) < 2:
                    continue
                cols = X[:, [fidx[t] for t in txs]].toarray()
                bin_tx = np.vstack([cols[b].sum(0) for b in bins])      # bins x isoforms
                bin_gene = bin_tx.sum(1)
                usable = bin_gene >= BARS["usage_min_bin_umi"]
                if usable.sum() < BARS["usage_min_used_bins"]:
                    continue
                usage = bin_tx[usable] / bin_gene[usable, None]
                pseudo = bin_tx.sum(0) / bin_tx.sum()
                part = np.nonzero(pseudo >= BARS["usage_min_isoform_frac"])[0]
                if len(part) < 2:
                    continue
                n_pass_iso += 1
                wu = {txs[i]: _window_unique_frac(txs[i], [o for o in g2tx[g] if o != txs[i]],
                                                  seqs, window, k) for i in part}
                if not all(v >= BARS["usage_min_unique_frac"] for v in wu.values()):
                    continue
                n_pass_floor += 1
                v = np.asarray(G[:, gidx[g]].todense()).ravel() / lib * 1e4
                q1, q5 = v[quint[0]].mean(), v[quint[-1]].mean()
                gfc = float(max(q1, q5) / max(min(q1, q5), 1e-9))
                xb = np.arange(len(bins))[usable]
                for i in part:
                    rho, pval = spearmanr(xb, usage[:, i])
                    results.append(dict(lineage=lin, gene_id=g, gene=g2name.get(g, g),
                                        tx_id=txs[i], rho=float(rho), p=float(pval),
                                        swing=float(usage[:, i].max() - usage[:, i].min()),
                                        window_unique=round(wu[txs[i]], 4),
                                        pseudobulk_usage=round(float(pseudo[i]), 4),
                                        gene_quintile_fc=gfc))
                    curves.append(pd.DataFrame(
                        {"lineage": lin, "gene": g2name.get(g, g), "tx_id": txs[i],
                         "bin": xb, "usage": usage[:, i]}))
            screened[lin] = {"candidates_expr": len(candidates), "pass_isoform_screen": n_pass_iso,
                             "pass_floor_window": n_pass_floor}

        res = pd.DataFrame(results)
        # One BH correction across every (lineage, gene, isoform) test in the study.
        if len(res):
            res["fdr"] = multipletests(res["p"], method="fdr_bh")[1]
            res["qualifies"] = ((res["fdr"] < BARS["usage_max_fdr"]) &
                                (res["swing"] >= BARS["usage_min_swing"]) &
                                (res["gene_quintile_fc"] <= BARS["usage_max_gene_fc"]))
        else:
            res["fdr"], res["qualifies"] = [], []
        out = SILO / "usage_shifts.tsv"
        _fresh(out)
        res.sort_values("fdr").to_csv(out, sep="\t", index=False)
        _made(out)
        cv = SILO / "usage_curves.tsv"
        _fresh(cv)
        (pd.concat(curves) if curves else pd.DataFrame(
            columns=["gene", "tx_id", "bin", "usage"])).to_csv(cv, sep="\t", index=False)
        qual = sorted(set(res.loc[res["qualifies"], "gene"])) if len(res) else []
        per_lin = ({lin: sorted(set(d.loc[d["qualifies"], "gene"]))
                    for lin, d in res.groupby("lineage")} if len(res) else {})
        (SILO / "qualifying_genes.tsv").write_text("\n".join(qual) + "\n")
        st3 = SILO / "Supp.Table3.tsv"
        _fresh(st3)
        (res[res["qualifies"]].sort_values(["lineage", "gene", "fdr"]) if len(res) else res
         ).to_csv(st3, sep="\t", index=False)
        return StageResult(
            metrics={"screened": screened, "tested_isoforms": len(res),
                     "qualifying_genes": qual, "n_qualifying": len(qual),
                     "qualifying_by_lineage": per_lin},
            passed=None,
            note=f"{len(qual)} qualifying genes ("
                 + ", ".join(f"{l}:{len(g)}" for l, g in per_lin.items()) + ")")


# ---------- s08: depth control ----------

class S08DepthControl(Stage):
    name, bar = "s08_depth_control", ("report-only: partial Spearman of usage vs axis "
                                      "controlling bin mean library size")
    asks = "Do the qualifying usage shifts survive control for sequencing depth?"

    # A usage shift along an axis that tracks library size can be a depth artifact; each
    # qualifying isoform is re-scored as a partial correlation controlling bin library size.

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd
        from scipy.stats import rankdata, spearmanr

        res = pd.read_csv(SILO / "usage_shifts.tsv", sep="\t")
        qual = res[res["qualifies"]] if len(res) else res
        pt_all = pd.read_csv(SILO / "pseudotime.tsv", sep="\t", index_col=0)
        atx_all = ad.read_h5ad(SILO / "pbmc_typed_tx.h5ad")
        agene_all = ad.read_h5ad(SILO / "pbmc_typed_gene.h5ad")
        _, g2tx, _, _, _ = _tx_tables()

        def partial_spearman(x, y, z):
            rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
            res_x = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
            res_y = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
            return float(np.corrcoef(res_x, res_y)[0, 1])

        rows, lineage_depth = [], {}
        for lin in sorted(pt_all["lineage"].unique()):
            d = pt_all[pt_all.lineage == lin]
            atx, agene = atx_all[d.index], agene_all[d.index]
            lib = np.asarray(agene.X.sum(1)).ravel()
            o = np.argsort(d["pseudotime"].values)
            X, libo = atx.X.tocsc()[o, :], lib[o]
            lineage_depth[lin] = round(float(spearmanr(d["pseudotime"].values, lib).statistic), 3)
            fidx = {t: j for j, t in enumerate(atx.var_names)}
            bins = np.array_split(np.arange(X.shape[0]), BARS["usage_n_bins"])
            bin_lib = np.array([libo[b].mean() for b in bins])
            for _, r in qual[qual.lineage == lin].iterrows():
                txs = [t for t in g2tx[r.gene_id] if t in fidx]
                cols = X[:, [fidx[t] for t in txs]].toarray()
                bt = np.vstack([cols[b].sum(0) for b in bins])
                bg = bt.sum(1)
                ok = bg >= BARS["usage_min_bin_umi"]
                usage = bt[ok] / bg[ok, None]
                i = txs.index(r.tx_id)
                xb = np.arange(len(bins))[ok]
                rho_p = partial_spearman(xb, usage[:, i], bin_lib[ok])
                rows.append(dict(lineage=lin, gene=r.gene, tx_id=r.tx_id, rho_raw=round(r.rho, 3),
                                 rho_partial_depth=round(rho_p, 3),
                                 attenuation=round(1 - abs(rho_p) / max(abs(r.rho), 1e-9), 2),
                                 depth_robust=abs(rho_p) >= 0.5))
        df = pd.DataFrame(rows)
        out = SILO / "depth_control.tsv"
        _fresh(out)
        df.to_csv(out, sep="\t", index=False)
        _made(out)
        robust = sorted(set(df.loc[df.depth_robust, "gene"])) if len(df) else []
        (SILO / "qualifying_genes_depth_robust.tsv").write_text("\n".join(robust) + "\n")
        return StageResult(
            metrics={"axis_depth_spearman": lineage_depth, "n_tested": len(df),
                     "n_depth_robust_isoforms": int(df.depth_robust.sum()) if len(df) else 0,
                     "depth_robust_genes": robust},
            passed=None,
            note=f"{len(robust)} genes depth-robust; axis-depth rho {lineage_depth}")


# ---------- s09: figure-ready panel ----------

class S09Panel(Stage):
    name, bar = "s09_panel", "report-only: Fig. 6 renders (axis, cell-state ordering, exemplar usage curves)"
    asks = "Draw Fig. 6."

    def run(self) -> StageResult:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        pt = pd.read_csv(SILO / "pseudotime.tsv", sep="\t")
        res = pd.read_csv(SILO / "usage_shifts.tsv", sep="\t")
        qual = res[res["qualifies"]].sort_values("fdr") if len(res) else res

        fig = plt.figure(figsize=(7.2, 5.6))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])
        lins = sorted(pt["lineage"].unique())
        letters = iter("ABCD")

        def letter(ax):
            ax.text(-0.12, 1.04, next(letters), transform=ax.transAxes, fontsize=10,
                    fontweight="bold", ha="left", va="bottom")
        # Display only: raw DPT values are right-skewed, so the colour is the within-lineage
        # percentile, the rank quantity every test in s06/s07 uses. Raw values stay in
        # pseudotime.tsv.
        pt["pct"] = pt.groupby("lineage")["pseudotime"].rank(pct=True)

        for col, lin in enumerate(lins):
            d = pt[pt.lineage == lin]
            ax = fig.add_subplot(gs[0, col])
            s_ = ax.scatter(d["umap1"], d["umap2"], c=d["pct"], s=1, cmap="viridis", rasterized=True)
            cb = fig.colorbar(s_, ax=ax)
            cb.ax.tick_params(labelsize=7); cb.set_label("pseudotime percentile", fontsize=7)
            ax.set_title(f"{lin} lineage", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlabel("UMAP 1", fontsize=7); ax.set_ylabel("UMAP 2", fontsize=7)
            letter(ax)

        ax = fig.add_subplot(gs[0, len(lins)])
        letter(ax)
        order = [s for s in ["t_naive", "t_memory", "b_naive", "b_memory", "plasmablast"]
                 if (pt["subtype"] == s).sum() >= BARS["subtype_min_cells"]]
        ax.violinplot([pt.loc[pt.subtype == s, "pct"] for s in order], showmedians=True)
        pretty = {"t_naive": "naive T", "t_memory": "memory T", "b_naive": "naive B",
                  "b_memory": "memory B", "plasmablast": "plasmablast"}
        ax.set_xticks(range(1, len(order) + 1), [pretty.get(s, s) for s in order],
                      rotation=30, fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("pseudotime percentile", fontsize=8)
        ax.set_title("cell state", fontsize=9)

        # Bottom row: three depth-robust qualifying genes, isoform usage vs pseudotime bin.
        curves = pd.read_csv(SILO / "usage_curves.tsv", sep="\t")
        robust = set((SILO / "qualifying_genes_depth_robust.tsv").read_text().split())
        pick = qual[qual["gene"].isin(robust)] if len(qual) and robust else qual
        top = pick.drop_duplicates(["lineage", "gene"]).head(3) if len(pick) else pick
        for j, (_, row) in enumerate(top.iterrows()):
            axu = fig.add_subplot(gs[1, j])
            if j == 0:
                letter(axu)
            cur = curves[(curves["gene"] == row["gene"]) & (curves["lineage"] == row["lineage"])]
            keep = set(qual.loc[(qual["gene"] == row["gene"]) &
                                (qual["lineage"] == row["lineage"]), "tx_id"])
            cur = cur[cur["tx_id"].isin(keep)]
            for tx, cg in cur.groupby("tx_id"):
                axu.plot(cg["bin"], cg["usage"], marker="o", ms=2.5, label=tx)
            axu.set_ylim(0, 1)
            axu.legend(fontsize=6, frameon=False)
            axu.set_title(row["gene"], fontsize=9, fontstyle="italic")
            axu.set_xlabel("pseudotime bin", fontsize=8)
            axu.set_ylabel("isoform usage" if j == 0 else "", fontsize=8)
            axu.tick_params(labelsize=7)
        fig.tight_layout()
        out_png, out_pdf = SILO / "Fig6.png", SILO / "Fig6.pdf"
        _fresh(out_png, out_pdf)
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        plt.close(fig)
        _made(out_png, out_pdf)
        return StageResult(metrics={"panel": str(out_pdf)}, passed=None, note="Fig. 6 written")


STAGES = [S01Align(), S02Geometry(), S03Matrix(), S04Subtypes(), S05Pseudotime(),
          S06Ordering(), S07Usage(), S08DepthControl(), S09Panel()]


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
