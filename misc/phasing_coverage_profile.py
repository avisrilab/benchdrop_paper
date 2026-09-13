"""Read length and 3'-anchored coverage profile of the BenchDrop-seq long reads, behind the 382 bp
window (Methods 2.5).
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pysam


def pctl(xs, ps):
    """Pure-python percentiles."""
    if not xs:
        return {p: float("nan") for p in ps}
    s = sorted(xs)
    n = len(s)
    out = {}
    for p in ps:
        if n == 1:
            out[p] = float(s[0])
            continue
        k = (n - 1) * (p / 100.0)
        lo = int(k)
        hi = min(lo + 1, n - 1)
        out[p] = s[lo] + (s[hi] - s[lo]) * (k - lo)
    return out


def frac_ge(xs, t):
    return sum(1 for x in xs if x >= t) / len(xs) if xs else float("nan")


def fastq_lengths(fq: Path, max_reads: int):
    lens = []
    op = gzip.open if str(fq).endswith(".gz") else open
    with op(fq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lens.append(len(line.strip()))
                if len(lens) >= max_reads:
                    break
    return lens


def bam_profile(bam: Path, max_reads: int):
    rlen, span, covf, tpd, fpos = [], [], [], [], []
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as bf:
        for r in bf:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            tlen = bf.get_reference_length(r.reference_name)
            if not tlen or r.reference_length is None:
                continue
            rlen.append(r.infer_read_length() or r.query_length or 0)
            span.append(r.reference_length)
            covf.append(r.reference_length / tlen)
            tpd.append(tlen - r.reference_end)
            fpos.append(r.reference_start)
            if len(span) >= max_reads:
                break
    return dict(rlen=rlen, span=span, covf=covf, tpd=tpd, fpos=fpos)


def show(name, xs, unit="bp", frac_thresholds=()):
    q = pctl(xs, [10, 25, 50, 75, 90])
    print(f"  {name:16s} n={len(xs):>7,}  "
          f"p10={q[10]:.0f} p25={q[25]:.0f}  MED={q[50]:.0f}  p75={q[75]:.0f} p90={q[90]:.0f} {unit}")
    for t in frac_thresholds:
        print(f"      frac >= {t:>5}{unit}: {100*frac_ge(xs, t):5.1f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="transcriptome-aligned reads")
    ap.add_argument("--fastq", default=None, help="optional raw fastq for a length cross-check")
    ap.add_argument("--max-reads", type=int, default=200000)
    args = ap.parse_args()

    if args.fastq:
        fl = fastq_lengths(Path(args.fastq), args.max_reads)
        print(f"\n== RAW read length ({Path(args.fastq).name}, unbiased of alignment) ==")
        show("raw_read_len", fl, frac_thresholds=(500, 1000, 2000))

    p = bam_profile(Path(args.bam), args.max_reads)
    print(f"\n== TRANSCRIPTOME-ALIGNED profile ({Path(args.bam).name}) ==")
    show("read_len_full", p["rlen"], frac_thresholds=())
    show("span_bp", p["span"], frac_thresholds=(500, 1000, 2000, 3000))
    show("cov_frac", [c * 100 for c in p["covf"]], unit="%", frac_thresholds=())
    show("tp_dist_3end", p["tpd"], frac_thresholds=())
    show("five_start", p["fpos"], frac_thresholds=())


if __name__ == "__main__":
    main()
