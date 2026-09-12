#!/usr/bin/env python
"""Matched-bulk AUPRC (Methods 2.6): does Bagpiper detect which isoforms are cell-type-variable?
Per-cell-type isoform fractions from the PBMC long reads are scored against ENCODE
sorted-population bulk RNA-seq as the reference.
"""
import gzip
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score

FEED = Path(os.path.expandvars("$BENCHDROP_FEED"))
BENCH = FEED / "pbmc/count_lr"
BULK = FEED / "public/encode_immune/quant"
OUT = Path(os.path.expandvars("$BENCHDROP_RESULTS"))
CELLTYPES = FEED / "pbmc/pbmc_celltypes.tsv"   # barcode -> lineage (b_cell, t_cell, nk, mono), shipped input
POPS = ["b_cell", "t_cell", "nk", "mono"]        # ENCODE pops matched to the 4 SC lineages
MIN_GENE_TPM = 10.0          # bulk: gene must be expressed
MIN_GENE_SC = 20            # single-cell: gene reads per lineage to trust its usage
POS_THR, NEG_THR = 0.20, 0.05   # bulk usage range: variable (positive) vs stable (negative)
K = 31


def load_annotation():
    tx2gene = {}
    for ln in open(FEED / "annotation_maps/gencode_v32.tx.tsv"):
        if ln.startswith("tx_id"):
            continue
        f = ln.split("\t")
        tx2gene[f[0]] = f[1]
    gene_tx = defaultdict(list)
    for t, g in tx2gene.items():
        gene_tx[g].append(t)
    return tx2gene, gene_tx


def bulk_usage(tx2gene, gene_tx):
    """per-population isoform-usage fraction; per gene the range across POPS."""
    tpm = {p: {} for p in POPS}
    for p in POPS:
        for ln in open(BULK / p / "quant.sf"):
            if ln.startswith("Name"):
                continue
            f = ln.split("\t")
            tpm[p][f[0]] = float(f[3])
    label, brange = {}, {}
    for g, txs in gene_tx.items():
        if len(txs) < 2:
            continue
        gtpm = {p: sum(tpm[p].get(t, 0.0) for t in txs) for p in POPS}
        if sum(gtpm.values()) < MIN_GENE_TPM * len(POPS) or min(gtpm.values()) <= 0:
            continue
        for t in txs:
            usage = [tpm[p].get(t, 0.0) / gtpm[p] for p in POPS]
            r = max(usage) - min(usage)
            brange[t] = r
            if r >= POS_THR:
                label[t] = 1
            elif r <= NEG_THR:
                label[t] = 0
            # between -> ambiguous, left unlabeled
    return label, brange


def sc_usage(tx2gene, gene_tx):
    """Stream the long-read isoform matrix (barcodes x transcripts), pseudobulk per
    lineage, return per-transcript usage-range across the 4 lineages."""
    ct = {}
    for ln in open(CELLTYPES):
        if ln.startswith("barcode"):
            continue
        b, c = ln.rstrip("\n").split("\t")
        ct[b] = c
    feats = [l.strip() for l in gzip.open(BENCH / "features.tsv.gz", "rt")]
    bcs = [l.strip() for l in gzip.open(BENCH / "barcodes.tsv.gz", "rt")]
    keeprow = {}
    for i, b in enumerate(bcs):
        if b in ct:
            keeprow[i] = ct[b]
    # (lineage, transcript_col) -> summed reads
    s = defaultdict(float)
    with gzip.open(BENCH / "matrix.mtx.gz", "rt") as fh:
        for ln in fh:
            if ln.startswith("%"):
                continue
            break                                    # dims header
        for ln in fh:
            r, c, v = ln.split()
            r0 = int(r) - 1
            lin = keeprow.get(r0)
            if lin is not None:
                s[(lin, int(c) - 1)] += float(v)
    # per gene per lineage -> usage fractions -> range
    # gather per (gene, lineage) total and per (tx, lineage) reads
    tx_reads = defaultdict(lambda: {p: 0.0 for p in POPS})
    for (lin, col), val in s.items():
        tx_reads[feats[col]][lin] += val
    srange = {}
    for g, txs in gene_tx.items():
        if len(txs) < 2:
            continue
        gtot = {p: sum(tx_reads[t][p] for t in txs if t in tx_reads) for p in POPS}
        usable = [p for p in POPS if gtot[p] >= MIN_GENE_SC]
        if len(usable) < 2:
            continue
        for t in txs:
            if t not in tx_reads:
                continue
            usage = [tx_reads[t][p] / gtot[p] for p in usable]
            srange[t] = max(usage) - min(usage)
    return srange


def uniq_fraction(scored_tx, tx2gene, gene_tx):
    """Unique-31mer fraction of each scored isoform vs its gene siblings (from the
    transcriptome). Bounded to scored transcripts."""
    need = set(scored_tx)
    sib = set()
    for t in need:
        for s in gene_tx[tx2gene[t]]:
            sib.add(s)
    seq, cur, buf = {}, None, []
    for ln in open(FEED / "genome/transcriptome.fa"):
        if ln.startswith(">"):
            if cur in sib:
                seq[cur] = "".join(buf)
            cur, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if cur in sib:
        seq[cur] = "".join(buf)
    km = {t: {s[i:i + K] for i in range(len(s) - K + 1)} for t, s in seq.items()}
    uf = {}
    for t in need:
        others = [s for s in gene_tx[tx2gene[t]] if s != t and s in km]
        if t not in km or not km[t]:
            continue
        if not others:
            uf[t] = 1.0
        else:
            uf[t] = len(km[t] - set().union(*[km[o] for o in others])) / len(km[t])
    return uf


def auprc(y, score):
    y, score = np.array(y), np.array(score)
    ap = average_precision_score(y, score)
    return ap, y.mean(), len(y)


def main():
    tx2gene, gene_tx = load_annotation()
    blabel, brange = bulk_usage(tx2gene, gene_tx)
    srange = sc_usage(tx2gene, gene_tx)

    scored = [t for t in blabel if t in srange]
    y = [blabel[t] for t in scored]
    sc = [srange[t] for t in scored]
    ap, prev, n = auprc(y, sc)
    print(f"=== matched-bulk AUPRC (bagpiper vs ENCODE silver standard) ===")
    print(f"scored isoforms: {n} | positives (bulk-variable): {sum(y)} ({100*prev:.1f}%) = the random baseline")
    print(f"AUPRC = {ap:.3f}   (baseline {prev:.3f}; lift {ap/prev:.2f}x)")

    # stratify by the coverage criterion
    uf = uniq_fraction(scored, tx2gene, gene_tx)
    for lo, hi, name in [(0.10, 1.01, "ABOVE ~10% unique (reliable regime)"),
                         (0.0, 0.10, "BELOW 10% unique (the floor)")]:
        idx = [i for i, t in enumerate(scored) if lo <= uf.get(t, -1) < hi]
        if len(idx) > 20 and 0 < sum(y[i] for i in idx) < len(idx):
            a, p, m = auprc([y[i] for i in idx], [sc[i] for i in idx])
            print(f"  {name}: AUPRC {a:.3f} (baseline {p:.3f}, n={m})")

    print("\nCeiling: bulk fixes WHICH isoforms vary, not the per-cell CHOICE -> this validates "
          "isoform-switch DETECTION, not fraction accuracy (the abundance-not-choice limit).")
    import json
    (OUT / "auprc_bulk.json").write_text(json.dumps(
        {"auprc": ap, "baseline": prev, "n": n, "positives": int(sum(y))}, indent=2))


if __name__ == "__main__":
    main()
