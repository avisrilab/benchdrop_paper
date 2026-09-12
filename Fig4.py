"""Fig. 4 and Supplementary Fig. 6, rendered from the staircase surfaces (Methods 2.5).

Stages: load the staircase surfaces, sensitivity, exemplar genes, render.
"""
from __future__ import annotations
import argparse
import json
import matplotlib
import numpy as np
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter


REV = Path(os.path.expandvars("$BENCHDROP_RESULTS"))
STA = REV / "staircase"
GEOM = REV / "pseudotime" / "geom"
FLOOR_DTU = REV / "floor_dtu" / "floor_summary.json"
SILO = REV / "fig4"

SWEEP = STA / "sweep_results.tsv"
TEMPLATE_FA = STA / "template.fa"
TEMPLATE_META = STA / "template_meta.tsv"
CELLS = STA / "cells.rung1.tsv"
EXEMPLARS = GEOM / "exemplar_window_check.tsv"
GEOMETRY = GEOM / "geometry.json"

DEPTHS = [1, 2, 4, 8, 16, 32]
PERCELLGENE = {d: STA / f"percellgene.depth{d}.tsv" for d in DEPTHS}
PERCELLTX = {d: STA / f"percelltx.depth{d}.tsv" for d in DEPTHS}

# The criterion, as in staircase/identifiability.py: K = 31; a transcript's unique
# fraction = the share of its 31-mers absent from every ANNOTATED sibling in the template; the
# GENE measure = the MIN across the gene's isoforms. Detection = estimated within-gene fraction
# >= 5% (identifiability.py's scale-invariant rule).
K = 31
FLOOR = 0.10
THRESHOLDS = [0.05, 0.10, 0.20]
DETECT_FRAC = 0.05
BINS = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, 0.25), (0.25, 0.50), (0.50, 1.01)]
WINDOW_BP = 382                      # measured PBMC terminal span (Fig6.py, geometry.json)
EXEMPLAR_GENES = ["CD8A", "IGHM"]    # a gene below the criterion and one above it

BARS = {
    # ---- input identity
    "sweep_rows": 168,               # 6 depths x 4 (level, metric) surfaces x 7 bins
    "template_records": 1037,        # template.fa records == template_meta rows
    "n_cells": 10,

    # ---- s05: the numbers the text quotes must rederive from the surfaces the figure is
    # drawn from; if any does not, figure and text disagree.
    #   (1) "dominant-isoform detection rising from 0.63 at one read per gene per cell to 0.98
    #        at 32"                       -> tx / detect_dom / ALL at depth1 and depth32
    #   (2) "at 32 reads per gene per cell, dominant-isoform choice is 0.98-1.00 above the 10%
    #        floor but 0.90-0.95 below it" -> gene / dom_choice / depth32 per bin
    #   (3) "2,250 of 3,000 multi-isoform cell-genes (75%) fall below the 10% floor"
    #                                        -> gene-level bin counts at depth32
    #   (4) "median absolute error in the isoform fraction at 32 reads was 0.019 above the
    #        threshold and 0.036 below" -> tx_frac_err above/below at the 10% threshold, depth32
    "quote_tol": 0.006,              # values are quoted to two decimals
    "quote_detect_dom_d1": 0.63,
    "quote_detect_dom_d32": 0.98,
    "quote_choice_d32_below": (0.90, 0.95),
    "quote_choice_d32_above": (0.98, 1.00),
    "quote_below_floor_n": 2250,
    "quote_below_floor_total": 3000,
    "quote_frac_err_d32_above": 0.019,
    "quote_frac_err_d32_below": 0.036,
}


# ---------- stage machinery ----------

@dataclass
class StageResult:
    metrics: dict = field(default_factory=dict)
    passed: bool | None = None      # None for a report-only stage
    note: str = ""


class Stage:
    name: str = ""
    bar: str = ""

    def run(self) -> StageResult:
        raise NotImplementedError


# ---------- render.py ----------
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402


SURFACE = "#ffffff"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ABOVE, BELOW = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]   # ordinal, steps 250-700
BIN_LABELS = ["<1", "1-5", "5-10", "10-25", "25-50", "50-100"]              # % unique sequence

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.fontsize": 6.8, "legend.frameon": False, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "pdf.fonttype": 42,
})


def _bin_key(lo, hi):
    return f"{lo:.2f}-{hi:.2f}"


def _style(ax, ygrid=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if ygrid:
        ax.grid(axis="y", linestyle="-")
        ax.set_axisbelow(True)


def _line(ax, x, y, color, label=None, lw=2.0, ms=4.2, z=3):
    ax.plot(x, y, color=color, lw=lw, solid_joinstyle="round", solid_capstyle="round",
            label=label, zorder=z)
    ax.plot(x, y, linestyle="none", marker="o", ms=ms + 1.6, color=SURFACE, zorder=z + 1)
    ax.plot(x, y, linestyle="none", marker="o", ms=ms, color=color, zorder=z + 2)


def _panel_letter(ax, letter):
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=INK, ha="left", va="bottom")


def render_main(surf, sens, ex, out_pdf, out_png):
    d = DEPTHS
    fig, axes = plt.subplots(2, 2, figsize=(6.85, 5.4))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.10, wspace=0.42, hspace=0.62)
    (a, b), (c, e) = axes
    s10 = {int(k): v["0.10"] for k, v in sens["by_depth"].items()}

    # (A) the depth axis: designed-dominant isoform detection vs reads/gene/cell
    _panel_letter(a, "A")
    _style(a)
    all_line = [surf["tx/detect_dom"][str(x)]["ALL"]["value"] for x in d]
    _line(a, d, [s10[x]["tx_detect_dom_above"] for x in d], ABOVE, "above 10% unique sequence")
    _line(a, d, [s10[x]["tx_detect_dom_below"] for x in d], BELOW, "below 10% unique sequence")
    a.plot(d, all_line, color=MUTED, lw=1.2, label="all transcripts", zorder=2)
    a.set_xscale("log", base=2)
    a.set_xticks(d)
    a.set_xticklabels([str(x) for x in d])
    a.set_xlim(0.85, 40)
    a.set_ylim(0.4, 1.02)
    a.set_yticks([0.4, 0.6, 0.8, 1.0])
    a.set_xlabel("reads per gene per cell")
    a.set_ylabel("dominant isoform detected")
    a.legend(loc="lower right", handlelength=1.6)

    # (B) the structural axis: gene-level dominant-isoform choice by unique-sequence bin
    _panel_letter(b, "B")
    _style(b)
    xs = list(range(len(BINS)))
    for i, x in enumerate(d):
        ys = [surf["gene/dom_choice"][str(x)][_bin_key(lo, hi)]["value"] for lo, hi in BINS]
        _line(b, xs, ys, RAMP[i], f"{x} reads", lw=1.6, ms=3.4)
    b.axvline(2.5, color=MUTED, lw=0.8, zorder=1)
    b.text(2.44, 0.562, "10% threshold", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    b.set_xticks(xs)
    b.set_xticklabels(BIN_LABELS)
    b.set_ylim(0.55, 1.02)
    b.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    b.set_xlabel("unique sequence per isoform (% of 31-mers)")
    b.set_ylabel("correct dominant isoform")
    b.legend(loc="lower right", ncol=2, handlelength=1.4, columnspacing=0.8,
             title="reads / gene / cell", title_fontsize=6.5)

    # (C) isoform-fraction error vs depth, above vs below the threshold
    _panel_letter(c, "C")
    _style(c)
    _line(c, d, [s10[x]["tx_frac_err_above"] for x in d], ABOVE, "above 10% unique sequence")
    _line(c, d, [s10[x]["tx_frac_err_below"] for x in d], BELOW, "below 10% unique sequence")
    c.set_xscale("log", base=2)
    c.set_xticks(d)
    c.set_xticklabels([str(x) for x in d])
    c.set_xlim(0.85, 40)
    c.set_ylim(0, 0.16)
    c.set_yticks([0, 0.05, 0.10, 0.15])
    c.set_xlabel("reads per gene per cell")
    c.set_ylabel("median |isoform-fraction error|")
    c.legend(loc="upper right", handlelength=1.6)

    # (D) the two exemplars inside the measured 382 bp terminal window
    _panel_letter(e, "D")
    _style(e, ygrid=False)
    e.grid(axis="x", linestyle="-")
    e.set_axisbelow(True)
    rows, ylabels, colors, groups = [], [], [], []
    y = 0
    for g in EXEMPLAR_GENES:
        iso = ex["genes"][g]["isoforms"]
        groups.append((g, y, y + len(iso) - 1))
        for r in iso:
            rows.append((y, r["window_unique"]))
            ylabels.append(r["tx"].replace("ENST", "ENST "))
            colors.append(ABOVE if r["window_unique"] >= FLOOR else BELOW)
            y += 1
        y += 0.8
    for (yy, v), col, lab in zip(rows, colors, ylabels):
        e.barh(yy, max(v, 0.004), height=0.62, color=col, linewidth=0, zorder=3)
        inside = v > 0.6
        e.text(0.985 if inside else max(v, 0.004) + 0.012, yy, lab,
               transform=e.get_yaxis_transform() if inside else e.transData, fontsize=5.8,
               color=SURFACE if inside else INK2, ha="right" if inside else "left",
               va="center", zorder=4, clip_on=False)
    e.axvline(FLOOR, color=MUTED, lw=0.8, zorder=2)
    gap_y = groups[0][2] + 0.9          # the empty row between the two gene groups
    e.text(FLOOR + 0.012, gap_y, "10% threshold", fontsize=6.5, color=MUTED, ha="left",
           va="center")
    for g, y0, y1 in groups:
        e.text(-0.02, (y0 + y1) / 2, g, transform=e.get_yaxis_transform(), fontsize=7.5,
               color=INK, ha="right", va="center", fontweight="bold")
    e.set_yticks([])
    e.set_ylim(y - 0.8 - 0.5, -0.9)
    e.set_xlim(0, 1.0)
    e.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    e.set_xticklabels(["0", "25", "50", "75", "100"])
    e.set_xlabel(f"unique sequence in the terminal {ex['window_bp']} bp window (%)")
    e.legend(handles=[Line2D([], [], color=ABOVE, lw=6, label="clears the threshold"),
                      Line2D([], [], color=BELOW, lw=6, label="below the threshold")],
             loc="upper right", handlelength=1.2)

    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def render_supp(sens, out_pdf, out_png):
    d = DEPTHS
    fig, axes = plt.subplots(2, 3, figsize=(6.85, 4.2), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.12,
                        wspace=0.35, hspace=0.45)
    for j, T in enumerate(THRESHOLDS):
        key = f"{T:.2f}"
        top, bot = axes[0][j], axes[1][j]
        rec = {x: sens["by_depth"][str(x)][key] for x in d}
        for ax in (top, bot):
            _style(ax)
            ax.set_xscale("log", base=2)
            ax.set_xticks(d)
            ax.set_xticklabels([str(x) for x in d])
            ax.set_xlim(0.85, 40)
        _line(top, d, [rec[x]["gene_choice_above"] for x in d], ABOVE, "above", lw=1.8, ms=3.6)
        _line(top, d, [rec[x]["gene_choice_below"] for x in d], BELOW, "below", lw=1.8, ms=3.6)
        top.set_ylim(0.5, 1.02)
        top.set_yticks([0.5, 0.75, 1.0])
        top.set_title(f"{int(T * 100)}% threshold", loc="left", color=INK, fontsize=6.8)
        _line(bot, d, [rec[x]["tx_frac_err_above"] for x in d], ABOVE, "above", lw=1.8, ms=3.6)
        _line(bot, d, [rec[x]["tx_frac_err_below"] for x in d], BELOW, "below", lw=1.8, ms=3.6)
        bot.set_ylim(0, 0.16)
        bot.set_yticks([0, 0.05, 0.10, 0.15])
        bot.set_xlabel("reads per gene per cell")
        if j == 0:
            top.set_ylabel("correct dominant isoform")
            bot.set_ylabel("median |fraction error|")
            top.legend(loc="lower right", handlelength=1.4, title="unique sequence",
                       title_fontsize=6.5)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


class S05Figure(Stage):
    name, bar = "s05_figure", ("Fig. 4 and Supplementary Fig. 6 written; the numbers the text "
                               "quotes rederive from the drawn surfaces within quote_tol; the "
                               "re-binned 10% gene share at depth 32 equals 1 - below share")
    asks = "Draw the figures; do the numbers the text quotes rederive from what is drawn?"

    def run(self) -> StageResult:
        surf = json.load(open(SILO / "surfaces.json"))
        sens = json.load(open(SILO / "sensitivity.json"))
        ex = json.load(open(SILO / "exemplars.json"))
        outs = {k: SILO / k for k in ("Fig4.pdf", "Fig4.png", "Supp.Fig6.pdf", "Supp.Fig6.png",
                                        "figure_numbers.json")}
        for p in outs.values():
            p.unlink(missing_ok=True)
        render_main(surf, sens, ex, outs["Fig4.pdf"], outs["Fig4.png"])
        render_supp(sens, outs["Supp.Fig6.pdf"], outs["Supp.Fig6.png"])

        b = BARS
        tol = b["quote_tol"]
        d1 = surf["tx/detect_dom"]["1"]["ALL"]["value"]
        d32 = surf["tx/detect_dom"]["32"]["ALL"]["value"]
        g32 = surf["gene/dom_choice"]["32"]
        below_vals = [g32[_bin_key(lo, hi)]["value"] for lo, hi in BINS if hi <= FLOOR]
        above_vals = [g32[_bin_key(lo, hi)]["value"] for lo, hi in BINS if lo >= FLOOR]
        below_n = sum(g32[_bin_key(lo, hi)]["n"] for lo, hi in BINS if hi <= FLOOR)
        total_n = g32["ALL"]["n"]
        share_rebinned = sens["by_depth"]["32"]["0.10"]["gene_share_clearing"]
        err_above = sens["by_depth"]["32"]["0.10"]["tx_frac_err_above"]
        err_below = sens["by_depth"]["32"]["0.10"]["tx_frac_err_below"]
        checks = {
            "q1_detect_dom_d1": abs(d1 - b["quote_detect_dom_d1"]) <= tol,
            "q1_detect_dom_d32": abs(d32 - b["quote_detect_dom_d32"]) <= tol,
            "q2_choice_d32_below_in_band": all(b["quote_choice_d32_below"][0] - tol <= v
                                               <= b["quote_choice_d32_below"][1] + tol
                                               for v in below_vals),
            "q2_choice_d32_above_in_band": all(b["quote_choice_d32_above"][0] - tol <= v
                                               <= b["quote_choice_d32_above"][1] + tol
                                               for v in above_vals),
            "q3_below_floor_count": (below_n == b["quote_below_floor_n"]
                                     and total_n == b["quote_below_floor_total"]),
            "q4_frac_err_d32_above": abs(err_above - b["quote_frac_err_d32_above"]) <= tol,
            "q4_frac_err_d32_below": abs(err_below - b["quote_frac_err_d32_below"]) <= tol,
            "rebinned_share_matches": abs(share_rebinned - (1 - below_n / total_n)) <= 0.001,
            "files_written": all(p.exists() and p.stat().st_size > 0 for p in outs.values()
                                 if p.suffix != ".json"),
        }
        numbers = {
            "quotes": {"detect_dom_all_d1": d1, "detect_dom_all_d32": d32,
                       "choice_d32_below_bins": below_vals, "choice_d32_above_bins": above_vals,
                       "below_floor_cellgenes_d32": [below_n, total_n],
                       "frac_err_d32": [err_above, err_below],
                       "rebinned_gene_share_clearing_10pct_d32": share_rebinned},
            "panelA_10pct": {str(x): {k: sens["by_depth"][str(x)]["0.10"][k]
                                      for k in ("tx_detect_dom_above", "tx_detect_dom_below")}
                             for x in DEPTHS},
            "panelC_10pct": {str(x): {k: sens["by_depth"][str(x)]["0.10"][k]
                                      for k in ("tx_frac_err_above", "tx_frac_err_below")}
                             for x in DEPTHS},
            "panelD": {g: ex["genes"][g] for g in EXEMPLAR_GENES},
            "annotation_only": sens["annotation_only"],
            "real_data_reach": ex["real_data_reach"],
            "checks": checks,
        }
        json.dump(numbers, open(outs["figure_numbers.json"], "w"), indent=1)
        ok = all(checks.values())
        failed = [k for k, v in checks.items() if not v]
        return StageResult(metrics={"checks": checks, "detect_dom_all": [d1, d32],
                                    "below_floor": [below_n, total_n]},
                           passed=ok, note=("all rederive" if ok else f"FAILED: {failed}"))


def _fresh(*paths: Path):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _made(*paths: Path):
    for p in paths:
        if not Path(p).exists() or Path(p).stat().st_size == 0:
            raise RuntimeError(f"stage exited but produced no {p}")


def _read_template():
    """template.fa + template_meta.tsv -> (tx->seq, gene->[tx in file order])."""
    seqs, cur, buf = {}, None, []
    with open(TEMPLATE_FA) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = ln[1:].strip(), []
            else:
                buf.append(ln.strip())
    if cur:
        seqs[cur] = "".join(buf)
    gene_tx = defaultdict(list)
    with open(TEMPLATE_META) as fh:
        assert fh.readline().startswith("transcript_id")
        for ln in fh:
            t, g, i, l = ln.rstrip("\n").split("\t")
            gene_tx[g].append(t)
    return seqs, gene_tx


def _tx_unique(seqs: dict, gene_tx: dict) -> dict:
    """Per-transcript unique 31-mer fraction against ALL annotated siblings (identifiability.py)."""
    km = {t: {s[i:i + K] for i in range(len(s) - K + 1)} for t, s in seqs.items()}
    out = {}
    for g, txs in gene_tx.items():
        for t in txs:
            if not km.get(t):
                out[t] = 0.0
                continue
            others = [x for x in txs if x != t and km.get(x)]
            if not others:
                out[t] = 1.0
                continue
            out[t] = len(km[t] - set().union(*[km[o] for o in others])) / len(km[t])
    return out


# ---------- s01: are the staircase artifacts all here, at the size they should be? ----------

class S01Inputs(Stage):
    name, bar = "s01_inputs", ("sweep rows == 168 over 6 depths; per-depth tables present for all "
                               "6 depths; template records == meta rows == 1037; 10 cells; "
                               "exemplar table carries CD8A + IGHM; floor_summary.json present")
    asks = "Are the staircase artifacts all present and the size they should be?"

    def run(self) -> StageResult:
        SILO.mkdir(parents=True, exist_ok=True)
        rows = [ln for ln in open(SWEEP) if ln.strip()]
        depths = sorted({int(r.split("\t")[0].replace("depth", "")) for r in rows})
        n_fa = sum(1 for ln in open(TEMPLATE_FA) if ln.startswith(">"))
        n_meta = sum(1 for ln in open(TEMPLATE_META)) - 1
        n_cells = sum(1 for ln in open(CELLS)) - 1
        missing = [str(p) for d in DEPTHS for p in (PERCELLGENE[d], PERCELLTX[d])
                   if not p.exists()]
        ex_genes = {ln.split("\t")[0] for ln in open(EXEMPLARS) if not ln.startswith("gene")}
        m = {"sweep_rows": len(rows), "depths": depths, "template_fa": n_fa,
             "template_meta": n_meta, "n_cells": n_cells, "missing_tables": missing,
             "exemplar_genes": sorted(ex_genes), "floor_summary": FLOOR_DTU.exists()}
        ok = (len(rows) == BARS["sweep_rows"] and depths == DEPTHS and not missing
              and n_fa == n_meta == BARS["template_records"]
              and n_cells == BARS["n_cells"]
              and set(EXEMPLAR_GENES) <= ex_genes and FLOOR_DTU.exists())
        return StageResult(metrics=m, passed=ok,
                           note=f"{len(rows)} sweep rows, depths {depths}, {n_fa} template tx")


# ---------- s02: what do the four depth x unique-fraction surfaces say? ----------

class S02Surfaces(Stage):
    name, bar = "s02_surfaces", "report-only: the four surfaces, long form + nested json"
    asks = "What do the four depth x unique-fraction surfaces say?"

    def run(self) -> StageResult:
        out_tsv, out_json = SILO / "surfaces.tsv", SILO / "surfaces.json"
        _fresh(out_tsv, out_json)
        surf = defaultdict(lambda: defaultdict(dict))   # "level/metric" -> depth -> bin -> (v, n)
        with open(out_tsv, "w") as fo:
            fo.write("level\tmetric\tdepth\tbin\tn\tvalue\n")
            for ln in open(SWEEP):
                lab, lev, met, b, n, v = ln.rstrip("\n").split("\t")
                d = int(lab.replace("depth", ""))
                surf[f"{lev}/{met}"][d][b] = (float(v), int(n))
                fo.write(f"{lev}\t{met}\t{d}\t{b}\t{n}\t{v}\n")
        json.dump({k: {str(d): {b: {"value": v, "n": n} for b, (v, n) in bins.items()}
                       for d, bins in ds.items()} for k, ds in surf.items()},
                  open(out_json, "w"), indent=1)
        _made(out_tsv, out_json)
        g32 = surf["gene/dom_choice"][32]
        below_n = sum(n for b, (v, n) in g32.items() if b != "ALL" and float(b.split("-")[0]) < FLOOR)
        total_n = g32["ALL"][1]
        m = {"surfaces": sorted(surf), "tx_detect_dom_ALL_d1": surf["tx/detect_dom"][1]["ALL"][0],
             "tx_detect_dom_ALL_d32": surf["tx/detect_dom"][32]["ALL"][0],
             "gene_dom_choice_d32": {b: v for b, (v, n) in g32.items()},
             "gene_below_floor_n_d32": below_n, "gene_total_n_d32": total_n,
             "gene_below_floor_share_d32": round(below_n / total_n, 4)}
        return StageResult(metrics=m, passed=None,
                           note=f"below-floor cell-genes at depth 32: {below_n}/{total_n}")


# ---------- s03: what changes if the floor is 5% or 20% instead of 10%? ----------

class S03Sensitivity(Stage):
    name, bar = "s03_sensitivity", ("report-only: re-bin the persisted per-depth tables at each "
                                    "threshold; denominators as in identifiability.py")
    asks = "What changes if the floor is 5% or 20% instead of 10%?"

    def run(self) -> StageResult:
        out_tsv, out_json = SILO / "sensitivity.tsv", SILO / "sensitivity.json"
        _fresh(out_tsv, out_json)
        seqs, gene_tx = _read_template()
        tx_u = _tx_unique(seqs, gene_tx)
        gene_u = {g: min(tx_u[t] for t in txs) for g, txs in gene_tx.items() if len(txs) >= 2}
        cells, state = [], {}
        with open(CELLS) as fh:
            fh.readline()
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                cells.append(f[0])
                state[f[0]] = int(f[5])
        genes = sorted(gene_u)

        # annotation-only shares (no data at all): what a reader screens from the reference
        ann = {f"{T:.2f}": {"tx_share_clearing": round(float(np.mean([u >= T for u in tx_u.values()])), 4),
                            "gene_share_clearing": round(float(np.mean([u >= T for u in gene_u.values()])), 4)}
               for T in THRESHOLDS}

        rows, nested = [], defaultdict(dict)
        for d in DEPTHS:
            dom = {}
            with open(PERCELLGENE[d]) as fh:
                fh.readline()
                for ln in fh:
                    cid, g, n_iso, reads, mae, ok = ln.rstrip("\n").split("\t")
                    dom[(cid, g)] = int(ok)
            # DESIGNED cell-gene denominator: every cell x every multi-isoform gene; a cell-gene
            # absent from the table drew no reads and is a call not made (scores 0).
            cg = [(gene_u[g], dom.get((c, g), 0)) for c in cells for g in genes]
            # transcript level, identifiability.py's rows_tx: (own u, |frac err|, detected at
            # est_frac >= 5%, is the designed dominant). Detection is scored over the DESIGNED
            # dominant transcripts; fraction error over EVERY transcript row, as the surface does.
            tx_rows = []
            with open(PERCELLTX[d]) as fh:
                fh.readline()
                for ln in fh:
                    cid, g, t, tc, ec, tf, ef = ln.rstrip("\n").split("\t")
                    txs = gene_tx.get(g, [])
                    s = state.get(cid, 0)
                    is_dom = int(bool(txs) and txs[s if s < len(txs) else 0] == t)
                    if t in tx_u:
                        tx_rows.append((tx_u[t], abs(float(tf) - float(ef)),
                                        int(float(ef) >= DETECT_FRAC), is_dom))
            for T in THRESHOLDS:
                above = [ok for u, ok in cg if u >= T]
                below = [ok for u, ok in cg if u < T]
                ta = [r for r in tx_rows if r[0] >= T]
                tb = [r for r in tx_rows if r[0] < T]
                da = [r[2] for r in ta if r[3]]
                db = [r[2] for r in tb if r[3]]
                rec = {
                    "gene_share_clearing": round(len(above) / len(cg), 4),
                    "gene_choice_above": round(float(np.mean(above)), 4) if above else None,
                    "gene_choice_below": round(float(np.mean(below)), 4) if below else None,
                    "tx_detect_dom_above": round(float(np.mean(da)), 4) if da else None,
                    "tx_detect_dom_below": round(float(np.mean(db)), 4) if db else None,
                    "tx_frac_err_above": round(float(np.median([r[1] for r in ta])), 4) if ta else None,
                    "tx_frac_err_below": round(float(np.median([r[1] for r in tb])), 4) if tb else None,
                    "n_cellgenes": len(cg), "n_tx_rows": len(tx_rows),
                    "n_dom_tx": int(sum(r[3] for r in tx_rows)),
                }
                nested[str(d)][f"{T:.2f}"] = rec
                rows.append((d, T, rec))
        with open(out_tsv, "w") as fo:
            keys = ["gene_share_clearing", "gene_choice_above", "gene_choice_below",
                    "tx_detect_dom_above", "tx_detect_dom_below", "tx_frac_err_above",
                    "tx_frac_err_below", "n_cellgenes", "n_dom_tx"]
            fo.write("depth\tthreshold\t" + "\t".join(keys) + "\n")
            for d, T, rec in rows:
                fo.write(f"{d}\t{T:.2f}\t" + "\t".join(str(rec[k]) for k in keys) + "\n")
        json.dump({"annotation_only": ann, "by_depth": nested,
                   "n_multi_isoform_genes": len(genes), "n_transcripts": len(tx_u)},
                  open(out_json, "w"), indent=1)
        _made(out_tsv, out_json)
        m = {"annotation_only": ann, "n_genes": len(genes),
             "d32": nested["32"]}
        return StageResult(metrics=m, passed=None,
                           note=("share of multi-isoform genes clearing 5/10/20%: "
                                 + " / ".join(str(ann[k]["gene_share_clearing"]) for k in ann)))


# ---------- s04: what do CD8A and IGHM look like inside the 382 bp window? ----------

class S04Exemplars(Stage):
    name, bar = "s04_exemplars", "report-only: per-isoform window-unique fraction for the two exemplars"
    asks = "What do CD8A and IGHM look like inside the measured terminal window?"

    def run(self) -> StageResult:
        out = SILO / "exemplars.json"
        _fresh(out)
        geom = json.load(open(GEOMETRY))
        if geom["window_bp"] != WINDOW_BP:
            return StageResult(metrics={"window_bp": geom["window_bp"]}, passed=False,
                               note="geometry.json window != config WINDOW_BP")
        ex = defaultdict(list)
        with open(EXEMPLARS) as fh:
            assert fh.readline().startswith("gene")
            for ln in fh:
                g, gid, t, L, u = ln.rstrip("\n").split("\t")
                ex[g].append({"tx": t, "len": int(L), "window_unique": float(u)})
        floor_dtu = json.load(open(FLOOR_DTU))
        rec = {"window_bp": geom["window_bp"], "span": geom["span"],
               "genes": {g: {"isoforms": sorted(ex[g], key=lambda r: -r["len"]),
                             "min_window_unique": min(r["window_unique"] for r in ex[g]),
                             "max_window_unique": max(r["window_unique"] for r in ex[g]),
                             "clears_floor": max(r["window_unique"] for r in ex[g]) >= FLOOR}
                         for g in EXEMPLAR_GENES},
               "real_data_reach": floor_dtu}
        json.dump(rec, open(out, "w"), indent=1)
        _made(out)
        m = {g: {"max_window_unique": rec["genes"][g]["max_window_unique"],
                 "clears_floor": rec["genes"][g]["clears_floor"]} for g in EXEMPLAR_GENES}
        return StageResult(metrics=m, passed=None,
                           note=", ".join(f"{g} max u@{WINDOW_BP}={m[g]['max_window_unique']}"
                                          for g in EXEMPLAR_GENES))


STAGES = [S01Inputs(), S02Surfaces(), S03Sensitivity(), S04Exemplars(), S05Figure()]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from-stage", metavar="NAME", help="resume at this stage")
    args = ap.parse_args()
    names = [s.name for s in STAGES]
    stages = STAGES[names.index(args.from_stage):] if args.from_stage else STAGES
    SILO.mkdir(parents=True, exist_ok=True)
    metrics_path = SILO / "stage_metrics.json"
    metrics = json.loads(metrics_path.read_text()) if metrics_path.exists() else {}
    for stage in stages:
        print(f"[{stage.name}] {stage.asks}\n[{stage.name}] bar: {stage.bar}")
        t0 = perf_counter()
        result = stage.run()
        status = "ok" if result.passed is None else ("pass" if result.passed else "FAIL")
        print(f"[{stage.name}] {status} ({perf_counter() - t0:.1f}s) {result.note}".rstrip())
        metrics[stage.name] = result.metrics
        metrics_path.write_text(json.dumps(metrics, indent=2, default=str))
        if result.passed is False:
            sys.exit(f"{stage.name} failed its check; stopping")


if __name__ == "__main__":
    main()
