#!/usr/bin/env python3
"""Render Table 1 and Supplementary Table 2 from the scored benchmark outputs (metrics.py,
matched_basis.py).

Usage:  python3 table.py <OUT_dir> [--check]
"""
import sys
from pathlib import Path

# The numbers the manuscript quotes from this benchmark, as (label, expected, tolerance);
# `--check` passes only if every one reproduces from the scored outputs.
CLAIMS = [
    # Matched-annotation concordance, against the ambiguous-aware IsoQuant arm.
    ("matched bagpiper|bambu", 0.59, 0.005),
    ("matched bagpiper|isoquant_wa", 0.58, 0.005),
    ("matched isoquant_wa|bambu (reference)", 0.62, 0.005),
    ("transcripts_detected bagpiper (k)", 34.4, 0.05),
    ("transcripts_detected bambu (k)", 34.6, 0.05),
    # Efficiency: quoted from the CLEAN unique-only arm (see NO_WALL).
    ("wall bagpiper (min)", 13.0, 0.5),
    ("wall isoquant (h)", 2.0, 0.05),
    ("wall bambu (h)", 4.05, 0.05),
    ("peak_rss bagpiper (GB)", 1.1, 0.05),
    ("peak_rss isoquant (GB)", 4.3, 0.05),
    ("peak_rss isoquant_wa (GB)", 4.3, 0.05),      # prose: memory does not move with mode
    ("peak_rss bambu (GB)", 27.0, 0.5),
    ("speedup vs isoquant (x)", 9.1, 0.5),
    ("speedup vs bambu (x)", 18.1, 0.5),
    ("memory ratio vs isoquant (x)", 4.0, 0.2),
    ("memory ratio vs bambu (x)", 25.6, 0.5),
]

# Which aligner each arm has to run. Bagpiper is transcriptome-native; the others are
# genome-based and consume the same genome BAM.
ALIGN_OF = {"bagpiper": "txome_align", "isoquant": "genome_align",
            "isoquant_wa": "genome_align", "bambu": "genome_align"}

# `isoquant` is IsoQuant under its `unique_only` default; `isoquant_wa` is the same tool under
# ambiguous-aware quantification. The second arm is optional: everything below renders with or
# without it, so this runs unchanged before and after that rerun. When both are present the
# unique-only arm is relabelled, because "IsoQuant" unqualified would then be ambiguous to a reader.
BASE_TOOLS = ["bagpiper", "isoquant", "bambu"]
LABEL = {"bagpiper": "Bagpiper", "isoquant": "IsoQuant", "bambu": "bambu",
         "isoquant_wa": "IsoQuant (ambiguous-aware)"}


def read_metrics(path):
    """tidy tool/metric/value TSV -> {(tool, metric): float}."""
    d = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            tool, metric, value = ln.rstrip("\n").split("\t")
            d[(tool, metric)] = float(value)
    return d


def hhmm(seconds):
    """Wall clock in the unit a reader compares at a glance."""
    return f"{seconds / 60:.1f} min" if seconds < 3600 else f"{seconds / 3600:.2f} h"


def main():
    out = Path(sys.argv[1])
    check = "--check" in sys.argv[2:]
    com = read_metrics(out / "metrics_common.tsv")
    raw = read_metrics(out / "metrics.tsv")

    # The ambiguous-aware IsoQuant arm counts as present only once it has been scored.
    tools = list(BASE_TOOLS)
    if ("isoquant_wa|common", "cells") in com:
        tools.insert(2, "isoquant_wa")
        LABEL["isoquant"] = "IsoQuant (unique-only)"

    versions = (out / "versions.txt").read_text().strip() if (out / "versions.txt").exists() else ""
    rec = {}
    if (out / "bc" / "recovery.tsv").exists():
        for ln in (out / "bc" / "recovery.tsv").read_text().splitlines():
            k, v = ln.split("\t")
            rec[k] = float(v)

    # align+quant per arm, and the derived ratios the manuscript quotes
    wall = {t: raw[(t, "wall_s")] + raw[(ALIGN_OF[t], "wall_s")] for t in tools}
    mem = {t: raw[(t, "peak_rss_gb")] for t in tools}
    ncells = int(com[("common", "cells")])

    L = []
    L.append(f"**Table 1, per tool. Head-to-head on public PIP-seq plus Nanopore data, identical inputs.** "
             f"All arms consume one Bagpiper barcode assignment and the same GENCODE v32 "
             f"annotation, scored on the {ncells} cells called by all {len(tools)} arms.")
    L.append("")
    L.append("| Metric | " + " | ".join(LABEL[t] for t in tools) + " |")
    L.append("|---" * (len(tools) + 1) + "|")

    def row(name, fn):
        L.append(f"| {name} | " + " | ".join(fn(t) for t in tools) + " |")

    row("Cells scored", lambda t: f"{int(com[(t + '|common', 'cells')]):,}")
    row("Transcripts detected", lambda t: f"{int(com[(t + '|common', 'transcripts_detected')]):,}")
    row("Median transcripts per cell", lambda t: f"{com[(t + '|common', 'median_tx_per_cell')]:,.0f}")
    row("Median counts per cell", lambda t: f"{com[(t + '|common', 'median_counts_per_cell')]:,.0f}")
    # The ambiguous-aware arm's wall clock is not quoted: it includes IsoQuant's one-time genedb
    # build, which no other arm paid. Efficiency is taken from the unique-only run; peak memory
    # is reported for both.
    NO_WALL = {"isoquant_wa"}
    row("Wall clock, align plus quantification",
        lambda t: "not comparable" if t in NO_WALL else hhmm(wall[t]))
    row("Peak memory, quantification step (GB)", lambda t: f"{mem[t]:.2f}")
    L.append("")
    # Matched-annotation rescoring (matched_basis.py). IsoQuant emits only 25,346 of the
    # 199,138 annotated transcripts, so scoring each pair on its OWN shared annotation puts
    # Bagpiper|bambu on a 2.9x larger basis than the reference pair. The matched column
    # rescores every pair on IsoQuant's annotation, which is the like-for-like comparison.
    mb = {}
    if (out / "matched_basis.tsv").exists():
        with open(out / "matched_basis.tsv") as fh:
            next(fh)
            for ln in fh:
                p, basis, n, rho = ln.rstrip("\n").split("\t")
                mb[(p, basis)] = (int(n), float(rho))

    ref = LABEL["isoquant_wa" if "isoquant_wa" in tools else "isoquant"]
    L.append(f"**Table 1, pairwise. Pairwise agreement on the same {ncells} cells.** The {ref} / bambu "
             f"row is the reference pair: two established tools scored against each other on these "
             f"same inputs.")
    L.append("")
    head = ("| Tool pair | Co-detected transcripts | Pseudobulk Spearman (co-detected) "
            "| Per-cell Spearman | Detection Jaccard |")
    rule = "|---|---|---|---|---|"
    if mb:
        head = ("| Tool pair | Co-detected transcripts | Pseudobulk Spearman (co-detected) "
                "| Same, matched annotation | Per-cell Spearman | Detection Jaccard |")
        rule = "|---|---|---|---|---|---|"
    L.append(head)
    L.append(rule)
    # The reference pair is whichever two ESTABLISHED tools are compared like-for-like. With the
    # ambiguous-aware arm run, that is IsoQuant-wa vs bambu (EM against EM); the unique-only rows
    # stay in as disclosure, not as the yardstick.
    if "isoquant_wa" in tools:
        pairs = [("bagpiper", "bambu", ""), ("bagpiper", "isoquant_wa", ""),
                 ("isoquant_wa", "bambu", " (reference pair)"),
                 ("bagpiper", "isoquant", ""), ("isoquant", "bambu", "")]
    else:
        pairs = [("bagpiper", "bambu", ""), ("bagpiper", "isoquant", ""),
                 ("isoquant", "bambu", " (reference pair)")]
    for a, b, note in pairs:
        k = f"{a}|{b}" if (f"{a}|{b}", "codetected_tx") in com else f"{b}|{a}"
        cols = [f"{LABEL[a]} vs {LABEL[b]}{note}",
                f"{int(com[(k, 'codetected_tx')]):,}",
                f"{com[(k, 'pb_spearman_codet')]:.2f}"]
        if mb:
            hit = mb.get((k, "isoquant_annotation"))
            cols.append(f"{hit[1]:.2f} ({hit[0]:,} tx)" if hit else "n/a")
        cols += [f"{com[(k, 'cell_spearman')]:.2f}", f"{com[(k, 'jaccard')]:.2f}"]
        L.append("| " + " | ".join(cols) + " |")
    L.append("")
    foot = [f"Bagpiper 0.1.0 (tag), {versions}, bambu via bambu_renv."]
    if rec:
        foot.append(f"Barcode recovery is a shared constant across arms "
                    f"({rec['passed']:,.0f} of {rec['total']:,.0f} reads, {rec['recovery_pct']:.1f}%), "
                    f"set by the ONT barcode error rate in the source data and not by the quantifier.")
    foot.append(f"Peak memory in decimal GB; exact peak resident set "
                + ", ".join(f"{LABEL[t]} {mem[t] * 1e9:,.0f} B" for t in tools) + ".")
    if mb:
        n, rho = mb[("bagpiper|bambu", "isoquant_bambu_codetected")]
        foot.append("The matched-annotation column rescores every pair on IsoQuant's 25,346-transcript "
                    "output annotation, because IsoQuant emits only that subset and so caps the basis "
                    "of any pair it appears in. On the strictest basis, the exact transcripts IsoQuant "
                    f"and bambu co-detect, Bagpiper and bambu agree at {rho:.2f} over {n:,} transcripts; "
                    "that set is conditioned on IsoQuant having made a call, so it is reported as a "
                    "bound and not as the headline.")
    foot.append("Speed and memory ratios versus Bagpiper: "
                + ", ".join((f"{LABEL[t]} {wall[t] / wall['bagpiper']:.1f}x wall, " if t not in NO_WALL
                             else f"{LABEL[t]} ")
                            + f"{mem[t] / mem['bagpiper']:.1f}x memory"
                            for t in tools if t != "bagpiper") + ".")
    if NO_WALL & set(tools):
        foot.append("The ambiguous-aware arm is scored for concordance only, not efficiency: its "
                    "run included a one-time ~100 min gene-database build that no other arm paid, "
                    "on a machine under heavy concurrent load. IsoQuant's runtime figure is taken "
                    "from the clean unique-only run, whose quantification mode does not change how "
                    "long alignment and traversal take.")
    L.append("Notes. " + " ".join(foot))
    table = "\n".join(L)
    print(table)
    (out / "table_r1_concordance.md").write_text(table + "\n")
    print(f"\nwrote {out / 'table_r1_concordance.md'}", file=sys.stderr)

    if check:
        got = {
            "matched bagpiper|bambu": mb[("bagpiper|bambu", "isoquant_annotation")][1],
            "matched bagpiper|isoquant_wa": mb[("bagpiper|isoquant_wa", "isoquant_annotation")][1],
            "matched isoquant_wa|bambu (reference)": mb[("isoquant_wa|bambu",
                                                        "isoquant_annotation")][1],
            "transcripts_detected bagpiper (k)": com[("bagpiper|common", "transcripts_detected")] / 1e3,
            "transcripts_detected bambu (k)": com[("bambu|common", "transcripts_detected")] / 1e3,
            "wall bagpiper (min)": wall["bagpiper"] / 60,
            "wall isoquant (h)": wall["isoquant"] / 3600,
            "wall bambu (h)": wall["bambu"] / 3600,
            "peak_rss bagpiper (GB)": mem["bagpiper"],
            "peak_rss isoquant (GB)": mem["isoquant"],
            "peak_rss isoquant_wa (GB)": mem["isoquant_wa"],
            "peak_rss bambu (GB)": mem["bambu"],
            "speedup vs isoquant (x)": wall["isoquant"] / wall["bagpiper"],
            "speedup vs bambu (x)": wall["bambu"] / wall["bagpiper"],
            "memory ratio vs isoquant (x)": mem["isoquant"] / mem["bagpiper"],
            "memory ratio vs bambu (x)": mem["bambu"] / mem["bagpiper"],
        }
        print("\n== check: every number the manuscript quotes ==", file=sys.stderr)
        bad = 0
        for label, want, tol in CLAIMS:
            have = got[label]
            ok = abs(have - want) <= tol
            bad += not ok
            print(f"{'PASS' if ok else 'FAIL'}  {label:<38} prose={want:<8} computed={have:.4f}",
                  file=sys.stderr)
        print(f"\n{len(CLAIMS) - bad}/{len(CLAIMS)} claims reproduce", file=sys.stderr)
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
