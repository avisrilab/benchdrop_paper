#!/usr/bin/env python3
"""Inject Bagpiper's cell barcode into a genome BAM as the CB tag, so IsoQuant and bambu consume the
same per-cell assignment (Methods 2.9).

Usage: inject_cb.py <in.bam> <out.bam>
"""
import sys
import pysam


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: inject_cb.py <in.bam> <out.bam>")
    inb, outb = sys.argv[1], sys.argv[2]
    kept = dropped = 0
    with pysam.AlignmentFile(inb, "rb") as fin, \
            pysam.AlignmentFile(outb, "wb", template=fin) as fout:
        for r in fin:
            parts = r.query_name.rsplit("_", 2)          # origid, CB(16nt), UMI
            if len(parts) != 3 or len(parts[1]) != 16:
                dropped += 1
                continue
            r.set_tag("CB", parts[1], "Z")
            r.set_tag("UB", parts[2], "Z")
            fout.write(r)
            kept += 1
    print(f"[inject_cb] tagged {kept} reads, dropped {dropped} (no 16nt CB in name)", file=sys.stderr)
    if kept == 0:
        sys.exit("[inject_cb] FATAL: 0 reads carried a 16nt CB in the name. The aligned reads were "
                 "almost certainly the FAILED (barcodeless) set, not passed. Aborting under set -e so "
                 "an empty CB BAM does not silently void every downstream arm (IsoQuant + bambu).")


if __name__ == "__main__":
    main()
