#!/usr/bin/env python
"""Does recovery track unique sequence? Reported at the transcript and gene level.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import os
# STAIRCASE_DIR overrides the results dir (misc/error_control.py scores each arm in its own).
STA = Path(os.environ.get("STAIRCASE_DIR",
                          os.path.expandvars("$BENCHDROP_RESULTS/staircase")))
K = 31

# ---- template ----
txseq, cur, buf = {}, None, []
for ln in open(STA / "template.fa"):
    if ln.startswith(">"):
        if cur:
            txseq[cur] = "".join(buf)
        cur, buf = ln[1:].strip(), []
    else:
        buf.append(ln.strip())
if cur:
    txseq[cur] = "".join(buf)
gene_tx = defaultdict(list)
for ln in open(STA / "template_meta.tsv"):
    if ln.startswith("transcript_id"):
        continue
    t, g, i, l = ln.rstrip("\n").split("\t")
    gene_tx[g].append(t)

KM = {t: {s[i:i + K] for i in range(len(s) - K + 1)} for t, s in txseq.items()}


def uniq(t, others):
    """Fraction of t's k-mers absent from every transcript in `others`."""
    if not KM.get(t):
        return 0.0
    if not others:
        return 1.0
    return len(KM[t] - set().union(*[KM[o] for o in others if o in KM])) / len(KM[t])


# global (all-annotated) per-transcript unique fraction
tx_uniq = {}
for g, txs in gene_tx.items():
    for t in txs:
        tx_uniq[t] = uniq(t, [x for x in txs if x != t])

# ---- per-cell-transcript accuracy ----
# DESIGNED isoform class = a depth-INVARIANT denominator. The simulator gives the
# state-determined dominant isoform 80% of a gene's reads and splits 20% across the rest, so
# "is this the designed dominant?" is fixed by construction. Conditioning on "did it happen to
# receive a read" instead makes the denominator grow and get rarer with depth, which is what
# made detection appear to FALL as depth rose.
state = {}
for ln in open(STA / "cells.rung1.tsv"):
    if ln.startswith("cid"):
        continue
    f = ln.rstrip("\n").split("\t")
    state[f[0]] = int(f[5])

rows_tx = []                       # (uniq, abs frac error, detected, is_designed_dominant)
cellgene = defaultdict(list)       # (cid,g) -> [(tx, true_count, est_count)]
for ln in open(STA / "percelltx.rung1.tsv"):
    if ln.startswith("cid"):
        continue
    cid, g, t, tc, ec, tf, ef = ln.rstrip("\n").split("\t")
    tc, ec, tf, ef = int(tc), float(ec), float(tf), float(ef)
    cellgene[(cid, g)].append((t, tc, ec))
    txs = gene_tx.get(g, [])
    s = state.get(cid, 0)
    dom_idx = s if s < len(txs) else 0
    is_dom = int(bool(txs) and txs[dom_idx] == t)
    # Detection must be SCALE-INVARIANT: bagpiper's count output is per-cell normalized (rows
    # sum to 1), not read counts, so an absolute threshold like est>=1 read is meaningless, and
    # est>0 is trivially true for any isoform in an ambiguous equivalence class. Score detection
    # on the within-gene estimated FRACTION instead: called at >=5% of the gene.
    if t in tx_uniq:
        rows_tx.append((tx_uniq[t], abs(tf - ef), int(ef >= 0.05), is_dom))

# ---- gene-level identifiability, all-annotated vs expressed-only ----
dom = {}
for ln in open(STA / "percellgene.rung1.tsv"):
    if ln.startswith("cid"):
        continue
    cid, g, n_iso, reads, mae, ok = ln.rstrip("\n").split("\t")
    dom[(cid, g)] = int(ok)

# DESIGNED cell-gene denominator, matching the transcript level above. Iterating over OBSERVED
# cell-genes had two depth-dependent faults, both of which flattered low depth:
#   1. the denominator shrank. A gene that drew no reads simply vanished, so depth 1 scored 2,860
#      cell-genes and depth 32 scored 3,000, dropping exactly the hardest (lowest-count) cases.
#      That made depth 1 score HIGHER than depth 2 (0.860 vs 0.840), which is backwards.
#   2. the BIN ITSELF moved. u_all was taken over the transcripts that happened to appear, so at
#      low depth a gene was scored against a subset of its isoforms and landed in a
#      higher-unique-fraction bin than its annotation warrants.
# Both are fixed by taking the denominator and the unique fraction from the DESIGN: every gene in
# every cell, and u_all over the gene's ANNOTATED isoforms. A cell-gene with no reads is a
# dominant-isoform call that was not made, so it scores 0. u_all is then depth-invariant, which is
# what makes it the annotation-computable screening rule the prose claims it is.
cells_all = [ln.split("\t")[0]
             for ln in open(STA / "cells.rung1.tsv").read().splitlines()[1:] if ln.strip()]
genes_all = [g for g, txs in gene_tx.items() if len(txs) >= 2]

rows_gene = []                     # (min_uniq_all, min_uniq_expressed, dom_correct)
rows_gene_id = []                  # the (cid, gene) key of each rows_gene entry
for cid in cells_all:
    for g in genes_all:
        items = cellgene.get((cid, g), [])
        u_all = min(tx_uniq[t] for t in gene_tx[g] if t in tx_uniq)
        expressed = [t for t, tc, _ in items if tc > 0]
        if len(expressed) < 2:
            u_exp = 1.0            # <2 isoforms actually present: nothing to confuse it with
        else:
            u_exp = min(uniq(t, [x for x in expressed if x != t]) for t in expressed)
        rows_gene.append((u_all, u_exp, dom.get((cid, g), 0)))
        rows_gene_id.append((cid, g))
# Persist the per-cell-gene table so a downstream comparison can re-bin on a gene subset.
with open(STA / "cellgene_uniq.tsv", "w") as fo:
    fo.write("cid\tgene_id\tu_all\tu_exp\tdom_correct\n")
    for (cid, g), (u_all, u_exp, ok) in zip(rows_gene_id, rows_gene):
        fo.write(f"{cid}\t{g}\t{u_all:.5f}\t{u_exp:.5f}\t{ok}\n")

BINS = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, 0.25), (0.25, 0.50), (0.50, 1.01)]


LABEL = sys.argv[1] if len(sys.argv) > 1 else None
emit = []


def table(title, data, idx, val, fmt="{:.4f}", val_name="metric", level=None, metric=None):
    print(f"\n{title}")
    print(f"{'unique-seq fraction':>22} | {'n':>8} | {val_name:>18}")
    print("-" * 55)
    for lo, hi in BINS:
        sel = [d for d in data if lo <= d[idx] < hi]
        if not sel:
            continue
        v = val(sel)
        print(f"{lo:>10.2f}-{hi:<10.2f} | {len(sel):>8} | {fmt.format(v):>18}")
        if LABEL and level:
            emit.append((level, metric, f"{lo:.2f}-{hi:.2f}", len(sel), v))
    if LABEL and level and data:
        emit.append((level, metric, "ALL", len(data), val(data)))


print(f"transcripts scored: {len(tx_uniq)} | per-cell-transcript obs: {len(rows_tx)} | "
      f"cell-genes: {len(rows_gene)}")

# TRANSCRIPT level, on a designed (depth-invariant) denominator
dom_tx = [r for r in rows_tx if r[3] == 1]
minor_tx = [r for r in rows_tx if r[3] == 0]
table("TRANSCRIPT level, by the transcript's OWN unique fraction: isoform-fraction error",
      rows_tx, 0, lambda s: np.median([r[1] for r in s]), val_name="median |frac err|",
      level="tx", metric="frac_err")
table("TRANSCRIPT level: detection of DESIGNED-DOMINANT isoforms (80% of the gene)",
      dom_tx, 0, lambda s: np.mean([r[2] for r in s]), "{:.3f}", "detected",
      level="tx", metric="detect_dom")
table("TRANSCRIPT level: detection of DESIGNED-MINOR isoforms (fixed denominator)",
      minor_tx, 0, lambda s: np.mean([r[2] for r in s]), "{:.3f}", "detected",
      level="tx", metric="detect_minor")

# GENE level: all-annotated vs expressed-only
table("GENE level, min unique across ALL annotated isoforms: dominant-isoform recovery",
      rows_gene, 0, lambda s: np.mean([r[2] for r in s]), "{:.3f}", "dominant OK",
      level="gene", metric="dom_choice")
table("GENE level, min unique across EXPRESSED isoforms only: dominant-isoform recovery",
      rows_gene, 1, lambda s: np.mean([r[2] for r in s]), "{:.3f}", "dominant OK")

# which gene summary predicts better (Spearman-free: correlation of summary with correctness)
a = np.array([(r[0], r[2]) for r in rows_gene], dtype=float)
b = np.array([(r[1], r[2]) for r in rows_gene], dtype=float)
print(f"\ncorrelation with dominant-correct: all-annotated r={np.corrcoef(a[:,0],a[:,1])[0,1]:.3f} | "
      f"expressed-only r={np.corrcoef(b[:,0],b[:,1])[0,1]:.3f}")

if LABEL and emit:
    with open(STA / "sweep_results.tsv", "a") as fo:
        for lev, met, bin_, n, v in emit:
            fo.write(f"{LABEL}\t{lev}\t{met}\t{bin_}\t{n}\t{v:.4f}\n")
