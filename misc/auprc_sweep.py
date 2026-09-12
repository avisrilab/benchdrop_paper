#!/usr/bin/env python
"""Threshold-robustness sweep for the matched-bulk AUPRC (auprc_bulk.py, Methods 2.6).
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
POPS = ["b_cell", "t_cell", "nk", "mono"]
K = 31
DEF = dict(pos=0.20, neg=0.05, tpm=10.0, sc=20)   # defaults


def load_once():
    tx2gene, gene_tx = {}, defaultdict(list)
    for ln in open(FEED / "annotation_maps/gencode_v32.tx.tsv"):
        if ln.startswith("tx_id"):
            continue
        f = ln.split("\t")
        tx2gene[f[0]] = f[1]
        gene_tx[f[1]].append(f[0])
    tpm = {p: {} for p in POPS}
    for p in POPS:
        for ln in open(BULK / p / "quant.sf"):
            if not ln.startswith("Name"):
                f = ln.split("\t")
                tpm[p][f[0]] = float(f[3])
    ct = {}
    for ln in open(CELLTYPES):
        if not ln.startswith("barcode"):
            b, c = ln.rstrip("\n").split("\t")
            ct[b] = c
    feats = [l.strip() for l in gzip.open(BENCH / "features.tsv.gz", "rt")]
    bcs = [l.strip() for l in gzip.open(BENCH / "barcodes.tsv.gz", "rt")]
    keeprow = {i: ct[b] for i, b in enumerate(bcs) if b in ct}
    tx_reads = defaultdict(lambda: {p: 0.0 for p in POPS})
    with gzip.open(BENCH / "matrix.mtx.gz", "rt") as fh:
        for ln in fh:
            if ln.startswith("%"):
                continue
            break
        for ln in fh:
            r, c, v = ln.split()
            lin = keeprow.get(int(r) - 1)
            if lin is not None:
                tx_reads[feats[int(c) - 1]][lin] += float(v)
    return tx2gene, gene_tx, tpm, tx_reads


def uniq_once(tx2gene, gene_tx, tpm, tx_reads):
    # candidate superset: isoforms in >=2-iso genes, present in sc, with any bulk tpm
    cand = set()
    for g, txs in gene_tx.items():
        if len(txs) < 2:
            continue
        if any(t in tx_reads for t in txs) and any(tpm["b_cell"].get(t, 0) or tpm["t_cell"].get(t, 0)
                                                   or tpm["nk"].get(t, 0) or tpm["mono"].get(t, 0) for t in txs):
            cand.update(txs)
    seq, cur, buf = {}, None, []
    for ln in open(FEED / "genome/transcriptome.fa"):
        if ln.startswith(">"):
            if cur in cand:
                seq[cur] = "".join(buf)
            cur, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if cur in cand:
        seq[cur] = "".join(buf)
    km = {t: {s[i:i + K] for i in range(len(s) - K + 1)} for t, s in seq.items()}
    uf = {}
    for t in cand:
        others = [s for s in gene_tx[tx2gene[t]] if s != t and s in km]
        if t in km and km[t]:
            uf[t] = 1.0 if not others else \
                len(km[t] - set().union(*[km[o] for o in others])) / len(km[t])
    return uf


def score(gene_tx, tpm, tx_reads, uf, pos, neg, min_tpm, min_sc):
    blabel = {}
    for g, txs in gene_tx.items():
        if len(txs) < 2:
            continue
        gtpm = {p: sum(tpm[p].get(t, 0.0) for t in txs) for p in POPS}
        if sum(gtpm.values()) < min_tpm * len(POPS) or min(gtpm.values()) <= 0:
            continue
        for t in txs:
            u = [tpm[p].get(t, 0.0) / gtpm[p] for p in POPS]
            r = max(u) - min(u)
            if r >= pos:
                blabel[t] = 1
            elif r <= neg:
                blabel[t] = 0
    srange = {}
    for g, txs in gene_tx.items():
        if len(txs) < 2:
            continue
        gtot = {p: sum(tx_reads[t][p] for t in txs if t in tx_reads) for p in POPS}
        usable = [p for p in POPS if gtot[p] >= min_sc]
        if len(usable) < 2:
            continue
        for t in txs:
            if t in tx_reads:
                u = [tx_reads[t][p] / gtot[p] for p in usable]
                srange[t] = max(u) - min(u)
    scored = [t for t in blabel if t in srange]
    y = np.array([blabel[t] for t in scored])
    s = np.array([srange[t] for t in scored])
    if len(y) < 20 or y.sum() == 0 or y.sum() == len(y):
        return None
    ap = average_precision_score(y, s)
    above = [i for i, t in enumerate(scored) if uf.get(t, -1) >= 0.10]
    ya, sa = y[above], s[above]
    apa = average_precision_score(ya, sa) if 0 < ya.sum() < len(ya) else float("nan")
    return len(y), y.mean(), ap, apa, ya.mean()


def main():
    tx2gene, gene_tx, tpm, tx_reads = load_once()
    uf = uniq_once(tx2gene, gene_tx, tpm, tx_reads)
    print(f"cached: {len(tx_reads)} sc transcripts, {len(uf)} with uniqueness\n")
    print(f"{'param':<18}{'n':>7}{'prev':>7}{'AUPRC':>8}{'lift':>6} | {'AUPRC>floor':>11}{'base':>7}")
    print("-" * 66)

    def row(tag, **kw):
        p = {**DEF, **kw}
        r = score(gene_tx, tpm, tx_reads, uf, p["pos"], p["neg"], p["tpm"], p["sc"])
        if r:
            n, prev, ap, apa, ba = r
            print(f"{tag:<18}{n:>7}{prev:>7.3f}{ap:>8.3f}{ap/prev:>6.2f} | {apa:>11.3f}{ba:>7.3f}")

    row("DEFAULT")
    print("-- vary POS_THR --")
    for v in (0.15, 0.20, 0.25, 0.30):
        row(f"  pos={v}", pos=v)
    print("-- vary NEG_THR --")
    for v in (0.02, 0.05, 0.10):
        row(f"  neg={v}", neg=v)
    print("-- vary MIN_GENE_TPM --")
    for v in (5.0, 10.0, 20.0):
        row(f"  tpm={v}", tpm=v)
    print("-- vary MIN_GENE_SC --")
    for v in (10, 20, 50, 100):
        row(f"  sc={v}", sc=v)


if __name__ == "__main__":
    main()
