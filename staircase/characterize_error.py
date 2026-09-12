#!/usr/bin/env python
"""ONT error profile (substitution, insertion, deletion rates) from the real K562 BenchDrop-seq
alignments; consumes `samtools view` on stdin.
"""
import re
import sys
from collections import Counter

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
tot = Counter()
nreads = 0

for ln in sys.stdin:
    f = ln.rstrip("\n").split("\t")
    if len(f) < 11 or f[5] == "*":
        continue
    nm = None
    for t in f[11:]:
        if t.startswith("NM:i:"):
            nm = int(t[5:])
            break
    if nm is None:
        continue
    M = I = D = 0
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op in "M=X":
            M += n
        elif op == "I":
            I += n
        elif op == "D":
            D += n
    subs = nm - I - D
    if M <= 0 or subs < 0:
        continue
    tot["M"] += M
    tot["I"] += I
    tot["D"] += D
    tot["S"] += subs
    nreads += 1

M = tot["M"]
if not M:
    sys.exit("no usable alignments")
err = tot["S"] + tot["I"] + tot["D"]
print(f"reads profiled: {nreads:,} | aligned bases: {M:,}")
print(f"substitution rate: {tot['S'] / M:.4f}")
print(f"insertion rate:    {tot['I'] / M:.4f}")
print(f"deletion rate:     {tot['D'] / M:.4f}")
print(f"TOTAL error rate:  {err / M:.4f}")
print(f"indel share of all errors: {(tot['I'] + tot['D']) / err:.3f}")
