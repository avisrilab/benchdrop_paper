#!/usr/bin/env python3
"""Split a CB-tagged, coordinate-sorted genome BAM into per-cell BAMs (one per CB), for bambu's
per-sample single-cell quantification.

Usage: split_cells.py <in.cb.bam> <out_dir> [top_n=300] [min_reads=10]
"""
import collections
import glob
import os
import sys
import pysam


def keep(r):   # primary, mapped, CB-tagged -- exactly what bambu quantifies
    return not (r.is_secondary or r.is_supplementary or r.is_unmapped) and r.has_tag("CB")


def main():
    inb = sys.argv[1]
    out = sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    min_reads = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    os.makedirs(out, exist_ok=True)
    for old in glob.glob(os.path.join(out, "*.bam*")):   # start clean; stale cells would re-break bambu
        os.remove(old)

    # pass 1: count primary-mapped reads per CB
    counts = collections.Counter()
    with pysam.AlignmentFile(inb, "rb") as f:
        for r in f:
            if keep(r):
                counts[r.get_tag("CB")] += 1
    eligible = collections.Counter({cb: c for cb, c in counts.items() if c >= min_reads})
    top = {cb for cb, _ in eligible.most_common(top_n)}
    print(f"[split] {len(counts)} cells; {len(eligible)} clear >={min_reads} reads; keeping top {len(top)}",
          file=sys.stderr)

    # pass 2: stream (coord order) into one writer per kept CB; primary-mapped reads only
    writers = {}
    with pysam.AlignmentFile(inb, "rb") as f:
        for r in f:
            if not keep(r):
                continue
            cb = r.get_tag("CB")
            if cb not in top:
                continue
            w = writers.get(cb)
            if w is None:
                w = writers[cb] = pysam.AlignmentFile(os.path.join(out, f"{cb}.bam"), "wb", template=f)
            w.write(r)
    for w in writers.values():
        w.close()
    for cb in top:
        pysam.index(os.path.join(out, f"{cb}.bam"))
    print(f"[split] wrote {len(top)} per-cell BAMs -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
