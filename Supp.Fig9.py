"""Matched-bulk isoform-switch detection (Methods 2.6, Supplementary Fig. 9).
Per-lineage isoform usage from the PBMC long reads is scored against ENCODE sorted-population
bulk RNA-seq by the area under the precision-recall curve (AUPRC). Bulk fixes which isoforms
vary between populations, so this scores switch detection rather than per-cell fraction accuracy.

Stages: pseudobulk the labelled cells and load the bulk quantifications, score at the default
thresholds, sweep the thresholds, Supplementary Fig. 9.
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
SILO = REV / "auprc"

MATRIX_DIR = FEED / "pbmc/count_lr"                      # the PBMC EM transcript matrix
CELLTYPES = FEED / "pbmc/pbmc_celltypes.tsv"             # barcode -> lineage, shipped input
BULK = FEED / "public/encode_immune/quant"               # salmon quant.sf per sorted population
TX_MAP = FEED / "annotation_maps/gencode_v32.tx.tsv"     # tx -> gene

POPS = ["b_cell", "t_cell", "nk", "mono"]                # ENCODE populations matched to the lineages
TPM_COL = 3                                              # quant.sf: Name Length EffectiveLength TPM NumReads

# Default thresholds (Methods 2.6) and the one-at-a-time sweep around them.
DEFAULT = {"pos": 0.20,     # bulk usage range at or above this: switching (positive)
           "neg": 0.05,     # bulk usage range at or below this: not switching (negative)
           "tpm": 10.0,     # bulk: mean gene TPM per population at least this
           "sc": 20}        # single-cell: gene reads per lineage to trust its usage
SWEEP = {"pos": [0.15, 0.20, 0.25, 0.30],
         "neg": [0.02, 0.05, 0.10],
         "tpm": [5.0, 10.0, 20.0],
         "sc": [10, 20, 50, 100]}
SWEEP_TITLE = {"pos": "switching:\nbulk range at least",
               "neg": "not switching:\nbulk range at most",
               "tpm": "bulk gene TPM\nat least",
               "sc": "lineage gene reads\nat least"}

BARS = {
    # ---- s02/s03: the figure must reproduce the numbers Methods 2.6 quotes (rounded).
    "quote_auprc": 0.57,
    "quote_baseline": 0.39,
    "quote_sweep_min": 0.46,
    "quote_sweep_max": 0.68,
    "quote_lift_min": 1.3,
    "quote_lift_max": 1.7,
    "min_scored": 20,                   # fewer scored transcripts than this is not an AUPRC
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


def _tx_tables():
    """gencode_v32.tx.tsv -> (tx->gene, gene->[tx])."""
    tx2g, g2tx = {}, defaultdict(list)
    with open(TX_MAP) as fh:
        assert fh.readline().startswith("tx_id")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tx2g[f[0]] = f[1]
            g2tx[f[1]].append(f[0])
    return tx2g, g2tx


def _load_tables():
    """The s01 tables as {transcript: array over POPS}: reads per lineage, TPM per population."""
    import pandas as pd
    tables = []
    for name in ("sc_reads.tsv", "bulk_tpm.tsv"):
        df = pd.read_csv(SILO / name, sep="\t", index_col=0)[POPS]
        tables.append(dict(zip(df.index, df.to_numpy(dtype=float))))
    return tables


def _score(reads, bulk, g2tx, pos, neg, tpm, sc):
    """Bulk labels and single-cell usage range at one threshold setting.

    Returns (scored transcripts, labels, scores) over the transcripts that carry both a bulk
    label and a single-cell usage range.
    """
    zero = np.zeros(len(POPS))
    label, srange = {}, {}
    for g, txs in g2tx.items():
        if len(txs) < 2:
            continue
        # bulk: gene expressed in every population, usage range across the populations
        m = np.array([bulk.get(t, zero) for t in txs])
        gt = m.sum(axis=0)
        if gt.sum() >= tpm * len(POPS) and gt.min() > 0:
            u = m / gt
            for t, r in zip(txs, u.max(axis=1) - u.min(axis=1)):
                if r >= pos:
                    label[t] = 1
                elif r <= neg:
                    label[t] = 0
        # single cell: usage range across the lineages with enough gene reads
        stx = [t for t in txs if t in reads]
        if stx:
            m = np.array([reads[t] for t in stx])
            gr = m.sum(axis=0)
            usable = gr >= sc
            if usable.sum() >= 2:
                u = m[:, usable] / gr[usable]
                srange.update(zip(stx, u.max(axis=1) - u.min(axis=1)))
    scored = [t for t in label if t in srange]
    y = np.array([label[t] for t in scored], dtype=int)
    s = np.array([srange[t] for t in scored], dtype=float)
    return scored, y, s


def _auprc(y, s):
    from sklearn.metrics import average_precision_score
    if len(y) < BARS["min_scored"] or y.sum() == 0 or y.sum() == len(y):
        return float("nan"), float("nan"), int(len(y))
    return float(average_precision_score(y, s)), float(y.mean()), int(len(y))


def _rounds_to(x: float, quote: float, ndigits: int) -> bool:
    return round(float(x), ndigits) == round(quote, ndigits)


# ---------- s01: pseudobulk the labelled cells, load the bulk quantifications ----------

class S01Tables(Stage):
    name, bar = "s01_tables", ("report-only: reads per transcript per lineage over the labelled "
                               "cells, and TPM per transcript per sorted population")
    asks = "What do the four lineages and the four sorted populations express, per isoform?"

    def run(self) -> StageResult:
        import pandas as pd

        sc_out, bulk_out = SILO / "sc_reads.tsv", SILO / "bulk_tpm.tsv"
        _fresh(sc_out, bulk_out)

        lineage = {}
        with open(CELLTYPES) as fh:
            fh.readline()
            for ln in fh:
                b, c = ln.rstrip("\n").split("\t")[:2]
                lineage[b] = c
        feats = [l.strip() for l in gzip.open(MATRIX_DIR / "features.tsv.gz", "rt")]
        bcs = [l.strip() for l in gzip.open(MATRIX_DIR / "barcodes.tsv.gz", "rt")]
        keeprow = {i: lineage[b] for i, b in enumerate(bcs) if b in lineage}

        # Rows are barcodes, columns are transcripts; stream the triplets, summing per lineage.
        reads = defaultdict(lambda: dict.fromkeys(POPS, 0.0))
        with gzip.open(MATRIX_DIR / "matrix.mtx.gz", "rt") as fh:
            for ln in fh:
                if not ln.startswith("%"):
                    break                                # the dimensions header
            for ln in fh:
                r, c, v = ln.split()
                lin = keeprow.get(int(r) - 1)
                if lin is not None:
                    reads[feats[int(c) - 1]][lin] += float(v)
        sc = pd.DataFrame.from_dict(reads, orient="index")[POPS]
        sc.index.name = "transcript"
        sc.to_csv(sc_out, sep="\t")

        tpm = {}
        for p in POPS:
            col = {}
            with open(BULK / p / "quant.sf") as fh:
                fh.readline()
                for ln in fh:
                    f = ln.split("\t")
                    col[f[0]] = float(f[TPM_COL])
            tpm[p] = col
        bulk = pd.DataFrame(tpm).fillna(0.0)[POPS]
        bulk.index.name = "transcript"
        bulk.to_csv(bulk_out, sep="\t")
        _made(sc_out, bulk_out)

        cells = pd.Series(list(keeprow.values())).value_counts().to_dict()
        return StageResult(metrics={"cells_labelled": len(lineage), "cells_in_matrix": len(keeprow),
                                    "cells_per_lineage": cells, "sc_transcripts": int(len(sc)),
                                    "bulk_transcripts": int(len(bulk))},
                           passed=None,
                           note=f"{len(keeprow):,} labelled cells in the matrix, "
                                f"{len(sc):,} transcripts with reads")


# ---------- s02: score at the default thresholds ----------

class S02Score(Stage):
    name, bar = "s02_score", ("AUPRC rounds to the quoted 0.57 and the prevalence baseline to the "
                              "quoted 0.39 (Methods 2.6)")
    asks = "Does the single-cell usage range rank the bulk-variable isoforms first?"

    def run(self) -> StageResult:
        import pandas as pd

        out = SILO / "scores.tsv"
        _fresh(out)
        tx2g, g2tx = _tx_tables()
        reads, bulk = _load_tables()
        scored, y, s = _score(reads, bulk, g2tx, **DEFAULT)
        ap, prev, n = _auprc(y, s)

        df = pd.DataFrame({"transcript": scored, "gene_id": [tx2g[t] for t in scored],
                           "bulk_switching": y, "sc_usage_range": s})
        df.to_csv(out, sep="\t", index=False)
        _made(out)

        ok = (_rounds_to(ap, BARS["quote_auprc"], 2)
              and _rounds_to(prev, BARS["quote_baseline"], 2))
        return StageResult(metrics={"n_scored": n, "positives": int(y.sum()),
                                    "auprc": round(ap, 4), "baseline": round(prev, 4),
                                    "lift": round(ap / prev, 3), "thresholds": DEFAULT},
                           passed=ok,
                           note=f"AUPRC {ap:.3f} (baseline {prev:.3f}, n={n:,})")


# ---------- s03: the threshold sweep ----------

class S03Sweep(Stage):
    name, bar = "s03_sweep", ("across the sweep, AUPRC min/max round to the quoted 0.46/0.68 and "
                              "the lift over baseline to 1.3/1.7 (Methods 2.6)")
    asks = "Is the AUPRC a property of the data or of the four thresholds?"

    def run(self) -> StageResult:
        import pandas as pd

        out = SILO / "sweep.tsv"
        _fresh(out)
        _, g2tx = _tx_tables()
        reads, bulk = _load_tables()
        rows = []
        for key, values in SWEEP.items():
            for v in values:
                setting = {**DEFAULT, key: v}
                _, y, s = _score(reads, bulk, g2tx, **setting)
                ap, prev, n = _auprc(y, s)
                rows.append({"varied": key, "value": v, **setting, "n_scored": n,
                             "positives": int(y.sum()), "baseline": prev, "auprc": ap,
                             "lift": ap / prev})
        df = pd.DataFrame(rows)
        df.to_csv(out, sep="\t", index=False)
        _made(out)

        lo, hi = df["auprc"].min(), df["auprc"].max()
        llo, lhi = df["lift"].min(), df["lift"].max()
        ok = (_rounds_to(lo, BARS["quote_sweep_min"], 2) and _rounds_to(hi, BARS["quote_sweep_max"], 2)
              and _rounds_to(llo, BARS["quote_lift_min"], 1) and _rounds_to(lhi, BARS["quote_lift_max"], 1))
        return StageResult(metrics={"settings": int(len(df)), "auprc_min": round(lo, 4),
                                    "auprc_max": round(hi, 4), "lift_min": round(llo, 3),
                                    "lift_max": round(lhi, 3)},
                           passed=ok,
                           note=f"AUPRC {lo:.3f} to {hi:.3f}, lift {llo:.2f} to {lhi:.2f}x "
                                f"over {len(df)} settings")


# ---------- s04: the panel ----------

class S04Panel(Stage):
    name, bar = "s04_panel", "report-only: Supplementary Fig. 9, the PR curve and the sweep"
    asks = "Draw Supplementary Fig. 9."

    def run(self) -> StageResult:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from sklearn.metrics import precision_recall_curve

        scores = pd.read_csv(SILO / "scores.tsv", sep="\t")
        sweep = pd.read_csv(SILO / "sweep.tsv", sep="\t")

        fig = plt.figure(figsize=(7.5, 2.9))
        gs = fig.add_gridspec(1, 5, width_ratios=[1.9, 1, 1, 1, 1], wspace=0.45,
                              left=0.07, right=0.98, top=0.8, bottom=0.2)

        # A: the precision-recall curve at the default thresholds
        ax = fig.add_subplot(gs[0, 0])
        y, s = scores["bulk_switching"].to_numpy(), scores["sc_usage_range"].to_numpy()
        ap, prev, n = _auprc(y, s)
        prec, rec, _ = precision_recall_curve(y, s)
        ax.plot(rec, prec, color="#000000", lw=1.2, label=f"AUPRC {ap:.2f} (n = {n:,} isoforms)")
        ax.axhline(prev, color="#000000", lw=0.8, ls=":", label=f"prevalence baseline {prev:.2f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("recall", fontsize=8)
        ax.set_ylabel("precision", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6.5, loc="upper right", frameon=False)
        ax.text(-0.18, 1.06, "A", transform=ax.transAxes, fontsize=11, fontweight="bold")

        # B: one-at-a-time threshold sweep
        axes = []
        for i, key in enumerate(SWEEP):
            ax = fig.add_subplot(gs[0, i + 1], sharey=axes[0] if axes else None)
            axes.append(ax)
            sub = sweep[sweep["varied"] == key].sort_values("value")
            x = np.arange(len(sub))
            ax.plot(x, sub["auprc"], "o-", color="#000000", ms=4, lw=1, label="AUPRC")
            ax.plot(x, sub["baseline"], "o--", color="#000000", ms=4, lw=1, mfc="white",
                    label="prevalence baseline")
            d = list(sub["value"]).index(DEFAULT[key])
            ax.axvspan(d - 0.3, d + 0.3, color="#dddddd", lw=0, zorder=0)
            ax.set_xticks(x)
            ax.set_xticklabels([f"{v:g}" for v in sub["value"]], fontsize=7)
            ax.set_xlim(-0.6, len(sub) - 0.4)
            ax.set_title(SWEEP_TITLE[key], fontsize=7)
            ax.tick_params(labelsize=7)
            if i == 0:
                ax.set_ylabel("AUPRC", fontsize=8)
                ax.text(-0.42, 1.12, "B", transform=ax.transAxes, fontsize=11, fontweight="bold")
                ax.legend(fontsize=5.5, loc="lower left", frameon=False)
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
        axes[0].set_ylim(0, 0.8)

        png, pdf = SILO / "Supp.Fig9.png", SILO / "Supp.Fig9.pdf"
        _fresh(png, pdf)
        fig.savefig(png, dpi=300)
        fig.savefig(pdf)
        plt.close(fig)
        _made(png, pdf)
        return StageResult(metrics={"panel": str(pdf)}, passed=None, note="Supplementary Fig. 9 written")


STAGES = [S01Tables(), S02Score(), S03Sweep(), S04Panel()]


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
