#!/usr/bin/env python
"""Assemble the coverage-criterion surfaces, one per (level, metric).
"""
import os
from collections import defaultdict
from pathlib import Path

STA = Path(os.path.expandvars("$BENCHDROP_RESULTS/staircase"))

TITLES = {
    ("gene", "dom_choice"): "GENE level: dominant-isoform CHOICE recovery (higher better)",
    ("tx", "detect_minor"): "TRANSCRIPT level: detection of DESIGNED-MINOR isoforms, "
                            "fixed denominator (higher better)",
    ("tx", "detect_dom"): "TRANSCRIPT level: detection of DESIGNED-DOMINANT isoforms (higher better)",
    ("tx", "frac_err"): "TRANSCRIPT level: median |isoform-fraction error| (LOWER better)",
}
ORDER = [("gene", "dom_choice"), ("tx", "detect_minor"), ("tx", "detect_dom"),
         ("tx", "frac_err")]

data = defaultdict(lambda: defaultdict(dict))   # (level,metric) -> depth -> bin -> (v,n)
bins_seen = defaultdict(list)
for ln in open(STA / "sweep_results.tsv"):
    lab, lev, met, b, n, v = ln.rstrip("\n").split("\t")
    key = (lev, met)
    if b not in bins_seen[key]:
        bins_seen[key].append(b)
    data[key][lab][b] = (float(v), int(n))


def depth_of(lab):
    try:
        return float(lab.replace("depth", ""))
    except ValueError:
        return 1e9


for key in ORDER:
    if key not in data:
        continue
    order = [b for b in bins_seen[key] if b != "ALL"] + (
        ["ALL"] if "ALL" in bins_seen[key] else [])
    print(f"\n{TITLES.get(key, key)}")
    print("rows = mean reads/gene/cell, cols = unique-sequence fraction\n")
    print(f"{'depth':<8}" + "".join(f"{b:>12}" for b in order))
    print("-" * (8 + 12 * len(order)))
    for lab in sorted(data[key], key=depth_of):
        row = "".join(
            f"{data[key][lab][b][0]:>12.3f}" if b in data[key][lab] else f"{'-':>12}"
            for b in order)
        print(f"{lab:<8}" + row)
