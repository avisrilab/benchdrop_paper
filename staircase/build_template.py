#!/usr/bin/env python
"""Endogenous multi-isoform template for the NanoSim simulation (Methods 2.5).
"""
import os
import random
import re
from collections import defaultdict
from pathlib import Path

GEN = Path(os.path.expandvars("$BENCHDROP_FEED/genome"))
OUT = Path(os.path.expandvars("$BENCHDROP_RESULTS/staircase"))
OUT.mkdir(parents=True, exist_ok=True)

MIN_ISO, MAX_ISO = 2, 5   # multi-isoform, capped so quantification stays interpretable
N_GENES = 300
SEED = 0

# 1. group transcripts by gene from the GTF; note gene biotype
gene_txs = defaultdict(set)
gene_type = {}
with open(GEN / "genes.gtf") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        a = f[8]
        g = re.search(r'gene_id "([^"]+)"', a).group(1)
        t = re.search(r'transcript_id "([^"]+)"', a).group(1)
        gt = re.search(r'gene_type "([^"]+)"', a)
        gene_txs[g].add(t)
        gene_type[g] = gt.group(1) if gt else ""

# 2. candidate multi-isoform protein-coding genes
cand = [g for g, ts in gene_txs.items()
        if gene_type.get(g) == "protein_coding" and MIN_ISO <= len(ts) <= MAX_ISO]
random.seed(SEED)
random.shuffle(cand)
want_tx = {t for g in cand for t in gene_txs[g]}
print(f"{len(cand)} candidate multi-isoform protein-coding genes")

# 3. pull sequences for candidate transcripts (stream the transcriptome once)
seq = {}
cur, buf = None, []
with open(GEN / "transcriptome.fa") as fh:
    for ln in fh:
        if ln.startswith(">"):
            if cur in want_tx:
                seq[cur] = "".join(buf)
            cur = ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if cur in want_tx:
        seq[cur] = "".join(buf)

# 4. keep genes that still have >= MIN_ISO transcripts WITH sequences; take N_GENES
sel = []
for g in cand:
    ts = [t for t in sorted(gene_txs[g]) if t in seq]
    if len(ts) >= MIN_ISO:
        sel.append((g, ts))
    if len(sel) >= N_GENES:
        break

with open(OUT / "template.fa", "w") as fo, open(OUT / "template_meta.tsv", "w") as fm:
    fm.write("transcript_id\tgene_id\tisoform_idx\tlength\n")
    n_tx = 0
    for g, ts in sel:
        for i, t in enumerate(ts):
            fo.write(f">{t}\n{seq[t]}\n")
            fm.write(f"{t}\t{g}\t{i}\t{len(seq[t])}\n")
            n_tx += 1
print(f"template: {len(sel)} genes, {n_tx} transcripts -> {OUT/'template.fa'}")
