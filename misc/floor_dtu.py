"""What fraction of PBMC multi-isoform genes clear the coverage criterion (Methods 2.5: 60.6% of
the 4,136 genes with two or more expressed isoforms, within the 382 bp terminal window).

Stages: check the matrix, list the multi-isoform gene population and the isoforms expressed,
apply the unique-sequence criterion under four counting rules, summarize on one denominator.
"""
from __future__ import annotations
import argparse
import gzip
import json
import numpy as np
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter


FEED = Path(os.path.expandvars("$BENCHDROP_FEED"))
REV = Path(os.path.expandvars("$BENCHDROP_RESULTS"))
SILO = REV / "floor_dtu"

# The PBMC EM matrix and its typed cells, the same artifacts Fig6.py uses.
MATRIX_DIR = FEED / "pbmc/count_lr"
TYPED_TX = REV / "pseudotime/pbmc_typed_tx.h5ad"     # 16,240 typed cells x 199,138 tx, built by Fig6.py
TRANSCRIPTOME = FEED / "genome/transcriptome.fa"     # GENCODE v32, bare-ENST records
TX_MAP = FEED / "annotation_maps/gencode_v32.tx.tsv"

EXPECT_MTX = (2_109_262, 199_138, 41_738_415)        # expected header; stop if it differs

# The floor rule, taken verbatim from staircase/identifiability.py rather than re-invented:
#   K = 31; a transcript's unique fraction = the share of its 31-mers absent from every
#   comparison sibling; the GENE-level measure is the MIN unique fraction across the gene's
#   isoforms, computed BOTH over all annotated siblings and over expressed-only siblings.
K = 31
FLOOR = 0.10

# An isoform "participates" at >= 10% pseudobulk usage within its gene, the same threshold as
# Fig6.py's `usage_min_isoform_frac`.
MIN_ISOFORM_FRAC = 0.10

# The terminal window a typical PBMC read spans, measured by Fig6.py (30M-read
# subset, 195,397 primary alignments): median span 382 bp. The windowed variant is the floor for a
# per-cell claim built on these reads, and it is the harshest of the three.
WINDOW_BP = 382

# Gene population for the primary track: protein-coding, multi-isoform, detected here.
MIN_GENE_UMI = 1_000        # Fig6.py's `usage_min_gene_umi`
BIOTYPE = "protein_coding"  # the staircase's 75% figure is over protein-coding genes

BARS = {
    # ---- data-identity checks.
    "mtx_header_exact": True,
    "min_typed_cells": 16_000,
    "min_txmap_frac": 0.99,

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




def _fresh(*paths: Path):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _made(*paths: Path):
    for p in paths:
        if not Path(p).exists() or Path(p).stat().st_size == 0:
            raise RuntimeError(f"stage exited but produced no {p}")


def _tx_tables():
    """gencode_v32.tx.tsv -> (tx->gene, gene->[tx], gene->name, gene->biotype)."""
    tx2g, g2tx, g2name, g2bio = {}, defaultdict(list), {}, {}
    with open(TX_MAP) as fh:
        assert fh.readline().startswith("tx_id")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tx, g, name, gbio = f[0], f[1], f[2], f[3]
            tx2g[tx] = g
            g2tx[g].append(tx)
            g2name[g] = name
            g2bio[g] = gbio
    return tx2g, g2tx, g2name, g2bio


def _read_fasta(path: Path, wanted: set) -> dict:
    """transcriptome.fa -> {bare ENST: sequence}, restricted to `wanted` to bound memory."""
    seqs, cur, buf = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur and cur in wanted:
                    seqs[cur] = "".join(buf)
                cur, buf = ln[1:].split()[0].strip(), []
            elif cur in wanted:
                buf.append(ln.strip())
    if cur and cur in wanted:
        seqs[cur] = "".join(buf)
    return seqs


def _kmers(seq: str, k: int, window: int | None = None) -> set:
    """The transcript's k-mer set; `window` restricts to the terminal (3') `window` bases."""
    s = seq[-window:] if window else seq
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def _gene_min_unique(txs: list, seqs: dict, k: int, window: int | None) -> float:
    """MIN over the gene's isoforms of (fraction of that isoform's k-mers absent from siblings).

    identifiability.py's GENE rule. A gene is callable only if EVERY participating isoform can be
    told apart, which is what a differential-USAGE claim between isoforms actually requires.
    """
    present = [t for t in txs if seqs.get(t)]
    if len(present) < 2:
        return 1.0 if len(present) == 1 else 0.0
    km = {t: _kmers(seqs[t], k, window) for t in present}
    km = {t: s for t, s in km.items() if s}
    if len(km) < 2:
        return 0.0
    out = []
    for t, s in km.items():
        others = set().union(*[v for u, v in km.items() if u != t])
        out.append(len(s - others) / len(s))
    return float(min(out))


# ---------- s01: is this the PBMC matrix, with its typed cells? ----------

class S01Identity(Stage):
    name, bar = "s01_identity", ("mtx header == expected dims; typed cells >= min_typed_cells; "
                                 "tx->gene coverage >= min_txmap_frac")
    asks = "Is this the PBMC matrix and its typed cells?"

    def run(self) -> StageResult:
        import anndata as ad

        SILO.mkdir(parents=True, exist_ok=True)
        with gzip.open(MATRIX_DIR / "matrix.mtx.gz", "rt") as fh:
            fh.readline(); fh.readline()
            dims = tuple(int(x) for x in fh.readline().split())
        if BARS["mtx_header_exact"] and dims != EXPECT_MTX:
            return StageResult(metrics={"mtx_header": dims}, passed=False,
                               note=f"header {dims} != expected {EXPECT_MTX}")

        a = ad.read_h5ad(TYPED_TX)
        tx2g, _, _, _ = _tx_tables()
        feats = list(a.var_names)
        frac = sum(1 for t in feats if t in tx2g) / len(feats)
        ok = a.n_obs >= BARS["min_typed_cells"] and frac >= BARS["min_txmap_frac"]
        return StageResult(metrics={"mtx_header": dims, "typed_cells": int(a.n_obs),
                                    "n_tx": len(feats), "txmap_frac": round(frac, 4)},
                           passed=ok, note=f"{a.n_obs:,} typed cells, txmap {frac:.3f}")


# ---------- s02: the multi-isoform gene population ----------

class S02Population(Stage):
    name, bar = "s02_population", ("report-only: protein-coding multi-isoform genes detected at "
                                   ">= MIN_GENE_UMI, with per-isoform pseudobulk usage")
    asks = "Which protein-coding genes have multiple isoforms here, and which isoforms are expressed?"

    def run(self) -> StageResult:
        import anndata as ad
        import pandas as pd

        out = SILO / "gene_population.tsv"
        _fresh(out)
        a = ad.read_h5ad(TYPED_TX)
        tx2g, g2tx, g2name, g2bio = _tx_tables()

        # Pseudobulk EM-UMIs per transcript across all typed cells.
        tot = np.asarray(a.X.sum(axis=0)).ravel()
        feats = list(a.var_names)
        by_gene_tx = defaultdict(dict)
        for i, t in enumerate(feats):
            g = tx2g.get(t)
            if g is not None and tot[i] > 0:
                by_gene_tx[g][t] = float(tot[i])

        rows = []
        for g, txcounts in by_gene_tx.items():
            if g2bio.get(g) != BIOTYPE:
                continue
            gsum = sum(txcounts.values())
            if gsum < MIN_GENE_UMI:
                continue
            n_annot = len(g2tx[g])
            if n_annot < 2:
                continue                      # single-isoform genes cannot show isoform usage
            expressed = sorted([t for t, v in txcounts.items()
                                if v / gsum >= MIN_ISOFORM_FRAC])
            rows.append({"gene_id": g, "gene_name": g2name.get(g, ""), "gene_umi": gsum,
                         "n_annotated": n_annot, "n_expressed": len(expressed),
                         "expressed": ",".join(expressed)})
        df = pd.DataFrame(rows).sort_values("gene_umi", ascending=False)
        df.to_csv(out, sep="\t", index=False)
        _made(out)
        multi_expr = int((df["n_expressed"] >= 2).sum())
        return StageResult(metrics={"genes_pc_multiisoform_detected": int(len(df)),
                                    "genes_with_2plus_expressed_isoforms": multi_expr,
                                    "min_gene_umi": MIN_GENE_UMI,
                                    "min_isoform_frac": MIN_ISOFORM_FRAC},
                           passed=None,
                           note=f"{len(df):,} protein-coding multi-isoform genes detected; "
                                f"{multi_expr:,} with >=2 expressed isoforms")


# ---------- s03: the floor, three variants ----------

class S03Floor(Stage):
    name, bar = "s03_floor", ("report-only: gene-level MIN unique 31-mer fraction, over all "
                              "annotated siblings, over expressed-only siblings, and over the "
                              "measured 382 bp terminal window")
    asks = "What fraction of that gene population clears the 10% identifiability floor?"

    def run(self) -> StageResult:
        import pandas as pd

        out = SILO / "floor_per_gene.tsv"
        _fresh(out)
        pop = pd.read_csv(SILO / "gene_population.tsv", sep="\t")
        _, g2tx, _, _ = _tx_tables()

        wanted = set()
        for g in pop["gene_id"]:
            wanted.update(g2tx[g])
        seqs = _read_fasta(TRANSCRIPTOME, wanted)

        rows = []
        for _, r in pop.iterrows():
            g = r["gene_id"]
            annot = g2tx[g]
            expressed = [t for t in str(r["expressed"]).split(",") if t]
            rows.append({
                "gene_id": g, "gene_name": r["gene_name"], "gene_umi": r["gene_umi"],
                "n_annotated": r["n_annotated"], "n_expressed": r["n_expressed"],
                "uniq_all": _gene_min_unique(annot, seqs, K, None),
                "uniq_expressed": _gene_min_unique(expressed, seqs, K, None),
                "uniq_window": _gene_min_unique(annot, seqs, K, WINDOW_BP),
                # The fourth rule, expressed siblings scored within the terminal window, is the
                # per-cell rule Methods 2.5 names.
                "uniq_window_expressed": _gene_min_unique(expressed, seqs, K, WINDOW_BP),
            })
        df = pd.DataFrame(rows)
        for col in ("uniq_all", "uniq_expressed", "uniq_window", "uniq_window_expressed"):
            df[f"pass_{col[5:]}"] = df[col] >= FLOOR
        df.to_csv(out, sep="\t", index=False)
        _made(out)

        n = len(df)
        # Genes with >=2 EXPRESSED isoforms are the ones a usage claim can even be made about.
        sub = df[df["n_expressed"] >= 2]
        m = {"n_genes": int(n),
             "pass_frac_all": round(float(df["pass_all"].mean()), 4),
             "pass_frac_expressed": round(float(df["pass_expressed"].mean()), 4),
             "pass_frac_window": round(float(df["pass_window"].mean()), 4),
             "pass_frac_window_expressed": round(float(df["pass_window_expressed"].mean()), 4),
             "n_multiexpressed": int(len(sub)),
             "pass_frac_expressed_multiexpressed": round(float(sub["pass_expressed"].mean()), 4)
             if len(sub) else None,
             # Every rule on the one denominator the paper quotes (genes with >= 2 expressed isoforms)
             "multiexpressed_pass_frac": {k: round(float(sub[f"pass_{k}"].mean()), 4)
                                         for k in ("all", "expressed", "window", "window_expressed")}
             if len(sub) else None}
        (SILO / "floor_summary.json").write_text(json.dumps(m, indent=2))
        return StageResult(metrics=m, passed=None,
                           note=f"pass rate: all-siblings {m['pass_frac_all']:.3f}, "
                                f"expressed-only {m['pass_frac_expressed']:.3f}, "
                                f"382bp window {m['pass_frac_window']:.3f}")


# ---------- s04: the four rates on the one denominator ----------

class S04Summary(Stage):
    name, bar = "s04_summary", ("report-only: the fraction of genes with >= 2 expressed isoforms "
                                "clearing the criterion under each of the four counting rules")
    asks = "What fraction of the multi-isoform gene population clears the criterion?"

    def run(self) -> StageResult:
        floor = json.loads((SILO / "floor_summary.json").read_text())
        # Rates are taken over genes with >= 2 EXPRESSED isoforms. A gene with one expressed
        # isoform has no sibling to be confused with and scores 1.0 by construction; counting it
        # would inflate the rate with genes no usage comparison applies to.
        rates = floor["multiexpressed_pass_frac"]
        m = {"n_multiexpressed_genes": floor["n_multiexpressed"],
             "rate_all_siblings_full_length": rates["all"],
             "rate_all_siblings_382bp_window": rates["window"],
             "rate_expressed_siblings_full_length": rates["expressed"],
             "rate_expressed_siblings_382bp_window": rates["window_expressed"]}
        return StageResult(metrics=m, passed=None,
                           note=(f"of {floor['n_multiexpressed']:,} genes with >= 2 expressed "
                                 f"isoforms, {rates['window_expressed']:.1%} clear the criterion "
                                 f"within the 382 bp window (expressed siblings); "
                                 f"{rates['window']:.1%} against all annotated siblings"))


STAGES = [S01Identity(), S02Population(), S03Floor(), S04Summary()]


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
