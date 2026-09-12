#!/usr/bin/env python
"""Evaluation: Bagpiper-recovered counts against the known per-cell isoform truth.
"""
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import os
# STAIRCASE_DIR overrides the results dir (misc/error_control.py scores each arm in its own).
STA = Path(os.environ.get("STAIRCASE_DIR",
                          os.path.expandvars("$BENCHDROP_RESULTS/staircase")))


def _open_maybe_gz(path):
    """Open the count-matrix file whether it is plain or gzipped."""
    path = Path(path)
    if path.exists():
        return open(path)
    gz = path.with_name(path.name + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    raise FileNotFoundError(f"neither {path} nor {gz}")

# CELL-CALL the recovered barcodes, then map CB -> cid (read name = cid.gene.tx_CB_UMI).
# ONT errors miscorrect some barcodes into OTHER valid whitelist combinations, producing many
# low-count phantom cells (672 recovered barcodes from 10 real ones at the empirical error
# rate). The real pipeline discards those by cell calling, so the benchmark must too. Summing
# phantom barcodes back into their true cell inflates the estimates and was what corrupted the
# earlier rung-2/rung-3 numbers. Separation is stark (thousands of reads vs a few hundred), so
# a 10%-of-max knee is sufficient.
cb_reads = Counter()
cb_cid = defaultdict(Counter)
for ln in gzip.open(STA / "out/barcode/passed.bcd.nanopore.fa.gz", "rt"):
    if not ln.startswith(">"):
        continue
    body, cb, _umi = ln[1:].split()[0].rsplit("_", 2)
    cb_reads[cb] += 1
    cb_cid[cb][body.split(".")[0]] += 1
mx = max(cb_reads.values()) if cb_reads else 0
called = {cb for cb, n in cb_reads.items() if n >= 0.10 * mx}
cb2cid = {cb: cb_cid[cb].most_common(1)[0][0] for cb in called}
print(f"barcodes recovered {len(cb_reads)} -> called cells {len(called)} "
      f"({sum(cb_reads[c] for c in called)}/{sum(cb_reads.values())} reads kept)",
      file=sys.stderr)

# bagpiper count matrix (barcodes x transcripts)
feats = [l.strip() for l in _open_maybe_gz(STA / "counts/features.tsv")]
bcs = [l.strip() for l in _open_maybe_gz(STA / "counts/barcodes.tsv")]
est = defaultdict(float)
with _open_maybe_gz(STA / "counts/matrix.mtx") as fh:
    for l in fh:
        if l.startswith("%"):
            continue
        break
    for l in fh:
        r, c, v = l.split()
        cid = cb2cid.get(bcs[int(r) - 1])
        if cid:
            est[(cid, feats[int(c) - 1])] += float(v)

# template gene map + truth + cell states
tx_gene = {}
for ln in open(STA / "template_meta.tsv"):
    if ln.startswith("transcript_id"):
        continue
    t, g, i, l = ln.rstrip("\n").split("\t")
    tx_gene[t] = g
tmpl_tx = set(tx_gene)
truth = defaultdict(float)
for ln in open(STA / "truth.rung1.tsv"):
    if ln.startswith("cid"):
        continue
    cid, t, n = ln.rstrip("\n").split("\t")
    truth[(cid, t)] = float(n)
cids = sorted({cid for cid, _ in truth})
gene_tx = defaultdict(list)
for t, g in tx_gene.items():
    gene_tx[g].append(t)

# detection + count correlation over template transcripts
TP = FP = FN = 0
et, tt = [], []
for cid in cids:
    for t in tmpl_tx:
        e, tr = est.get((cid, t), 0.0), truth.get((cid, t), 0.0)
        et.append(e)
        tt.append(tr)
        if tr > 0 and e > 0:
            TP += 1
        elif tr == 0 and e > 0:
            FP += 1
        elif tr > 0 and e == 0:
            FN += 1
prec = TP / (TP + FP) if TP + FP else 0
rec = TP / (TP + FN) if TP + FN else 0
r = np.corrcoef(et, tt)[0, 1]

# isoform-fraction MAE + dominant-isoform recovery, per (cell, gene) with truth
maes, dom_ok, dom_tot = [], 0, 0
pcg = []   # per-cell-gene rows, joined against isoform identifiability downstream
pct = []   # per-cell-TRANSCRIPT rows: accuracy at the transcript level
for cid in cids:
    for g, txs in gene_tx.items():
        tsum = sum(truth.get((cid, t), 0) for t in txs)
        if tsum == 0:
            continue
        esum = sum(est.get((cid, t), 0) for t in txs)
        tf = {t: truth.get((cid, t), 0) / tsum for t in txs}
        ef = {t: (est.get((cid, t), 0) / esum if esum else 0) for t in txs}
        mae = float(np.mean([abs(tf[t] - ef[t]) for t in txs]))
        ok = int(max(txs, key=lambda t: truth.get((cid, t), 0))
                 == max(txs, key=lambda t: est.get((cid, t), 0)))
        maes.append(mae)
        dom_tot += 1
        dom_ok += ok
        pcg.append((cid, g, len(txs), int(tsum), round(mae, 5), ok))
        for t in txs:
            pct.append((cid, g, t, int(truth.get((cid, t), 0)),
                        round(est.get((cid, t), 0.0), 4), round(tf[t], 5), round(ef[t], 5)))
with open(STA / "percellgene.rung1.tsv", "w") as fo:
    fo.write("cid\tgene_id\tn_iso\ttrue_reads\tmae\tdom_correct\n")
    for row in pcg:
        fo.write("\t".join(map(str, row)) + "\n")
with open(STA / "percelltx.rung1.tsv", "w") as fo:
    fo.write("cid\tgene_id\ttranscript_id\ttrue_count\test_count\ttrue_frac\test_frac\n")
    for row in pct:
        fo.write("\t".join(map(str, row)) + "\n")

print(f"cells {len(cids)} | template transcripts {len(tmpl_tx)}")
print(f"detection: precision {prec:.3f}  recall {rec:.3f}  (TP {TP} FP {FP} FN {FN})")
print(f"count correlation (est vs true): r {r:.3f}")
print(f"isoform-fraction MAE: mean {np.mean(maes):.4f}  median {np.median(maes):.4f}"
      f"  over {len(maes)} cell-genes")
print(f"dominant isoform recovered: {dom_ok}/{dom_tot} ({dom_ok / dom_tot:.3f})")
