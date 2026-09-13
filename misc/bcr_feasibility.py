#!/usr/bin/env python
"""The two B-cell checks on the PBMC long reads (Methods 2.7): light-chain exclusion and the IgM
membrane-versus-secreted form against cell state.

Step 1 writes bcr_feasibility_per_cell.tsv, the immunoglobulin heavy- and light-chain read counts
per barcode. Steps 2 to 4 join it to the Seurat anchor-transfer labels (transfer_labels.R) and
count the light-chain violations and the IgM form per B state.

Usage:  python bcr_feasibility.py
"""
import os
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

BASE = Path(os.path.expandvars("$BENCHDROP_FEED/pbmc"))
ANN = Path(os.path.expandvars("$BENCHDROP_FEED/annotation_maps/gencode_v32.tx.tsv"))
MATS = BASE / "count_lr"                              # the PBMC EM transcript matrix
WL = BASE / "count_sr/Gene/filtered/barcodes.tsv"     # the short-read called cells
OUT = Path(os.path.expandvars("$BENCHDROP_RESULTS"))
LABELS = OUT / "pseudotime" / "seurat_anchor_labels.tsv"   # written by transfer_labels.R
OUT.mkdir(parents=True, exist_ok=True)

MEMBRANE = "ENST00000637539"  # IGHM membrane form (BCR)
SECRETED = "ENST00000390559"  # IGHM secreted form (antibody)
LOCUS_BY_CHROM = {"chr14": "IGH", "chr2": "IGK", "chr22": "IGL"}

B_STATES = ["B naive", "B intermediate", "B memory", "Plasmablast"]
# Light-chain check: a cell enters with IGH_C >= 1 and light_C >= 1; a violation is the minor
# light-chain class at >= 2 reads and > 20% of light_C.
PURITY_MIN_READS = 1
VIOLATION_MINOR_READS = 2
VIOLATION_MINOR_FRAC = 0.20
# IgM form: a cell enters with mem + sec >= 3 (EM-allocated counts, not rounded); the form with
# more reads is called, and a tie is not called membrane.
IGM_MIN_SUM = 3


def load_ig_columns():
    """Map 0-based transcript-column index -> (group, transcript_id) for IG tx."""
    feats = [ln.split("\t")[0].strip() for ln in
             gzip.open(MATS / "features.tsv.gz", "rt")]
    col_of = {t: i for i, t in enumerate(feats)}
    ann = pd.read_csv(ANN, sep="\t").rename(
        columns={"tx_id": "transcript_id", "gene_biotype": "gene_type"})
    is_ig = ann["gene_type"].str.startswith("IG_", na=False) | \
        ann["gene_name"].str.match(r"^IG[HKL]", na=False)
    ig = ann[is_ig].copy()
    col2grp = {}
    for _, r in ig.iterrows():
        t = r["transcript_id"]
        if t not in col_of:
            continue
        locus = LOCUS_BY_CHROM.get(str(r["chrom"]))
        if locus is None:
            continue
        gt = str(r["gene_type"])
        seg = "C" if "_C_" in gt else "V" if "_V_" in gt else None
        if seg is None:  # fall back to name for the C genes (IGHM/IGKC/IGLC..)
            seg = "C" if r["gene_name"] in {
                "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",
                "IGHE", "IGHD", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC4",
                "IGLC5", "IGLC6", "IGLC7"} else None
        if seg is None:
            continue
        col2grp[col_of[t]] = (f"{locus}_{seg}", r["gene_name"])
    mem = col_of.get(MEMBRANE)
    sec = col_of.get(SECRETED)
    return col2grp, mem, sec, len(feats)


def per_cell_table():
    col2grp, mem_col, sec_col, n_feat = load_ig_columns()
    keep = set(col2grp) | {c for c in (mem_col, sec_col) if c is not None}
    ngrp = defaultdict(int)
    for _, (g, _) in col2grp.items():
        ngrp[g] += 1
    print(f"IG transcript columns indexed: {len(col2grp)} across "
          f"{n_feat:,} features", file=sys.stderr)
    for g in sorted(ngrp):
        print(f"  {g}: {ngrp[g]} transcripts", file=sys.stderr)
    print(f"  IGHM membrane col={mem_col}, secreted col={sec_col}",
          file=sys.stderr)

    # cell(row) -> {group: reads}; membrane/secreted counts
    cell = defaultdict(lambda: defaultdict(float))
    mtx = MATS / "matrix.mtx.gz"
    with gzip.open(mtx, "rt") as fh:
        for ln in fh:
            if ln.startswith("%"):
                continue
            break  # first non-comment line is the dims header; skip it
        n = 0
        for ln in fh:
            r, c, v = ln.split()
            c0 = int(c) - 1
            if c0 not in keep:
                continue
            row = int(r) - 1
            val = float(v)
            grp = col2grp.get(c0)
            if grp is not None:
                cell[row][grp[0]] += val
            if c0 == mem_col:
                cell[row]["_MEM"] += val
            elif c0 == sec_col:
                cell[row]["_SEC"] += val
            n += 1
    print(f"IG matrix entries streamed: {n:,}; barcodes with any IG read: "
          f"{len(cell):,}", file=sys.stderr)

    # barcode <-> row, and the short-read called-cell whitelist
    bcs = [ln.strip() for ln in gzip.open(MATS / "barcodes.tsv.gz", "rt")]
    wl = {ln.strip() for ln in open(WL)}
    # try direct match; if near-zero, strip a trailing -1
    ov = sum(1 for b in bcs if b in wl)
    if ov < 100:
        wl2 = {b.split("-")[0] for b in wl}
        ov = sum(1 for b in bcs if b.split("-")[0] in wl2)
        in_wl = lambda i: bcs[i].split("-")[0] in wl2
    else:
        in_wl = lambda i: bcs[i] in wl
    print(f"short-read whitelist cells: {len(wl):,}; overlap with LR barcodes: "
          f"{ov:,}", file=sys.stderr)

    rows = []
    for row, d in cell.items():
        hC = d.get("IGH_C", 0.0)
        kC, lC = d.get("IGK_C", 0.0), d.get("IGL_C", 0.0)
        hV = d.get("IGH_V", 0.0)
        kV, lV = d.get("IGK_V", 0.0), d.get("IGL_V", 0.0)
        rows.append(dict(
            barcode=bcs[row], in_whitelist=in_wl(row),
            IGH_C=hC, IGK_C=kC, IGL_C=lC, IGH_V=hV, IGK_V=kV, IGL_V=lV,
            light_C=kC + lC, light_V=kV + lV, ig_total=hC + kC + lC + hV + kV + lV,
            mem=d.get("_MEM", 0.0), sec=d.get("_SEC", 0.0)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "bcr_feasibility_per_cell.tsv", sep="\t", index=False)

    def report(sub, tag):
        print(f"\n===== {tag} (n barcodes with IG signal = {len(sub):,}) =====")
        for thr in (1, 2):
            heavy = sub.IGH_C >= thr
            light = sub.light_C >= thr
            hv = sub.IGH_V >= thr
            lv = sub.light_V >= thr
            isotyped = heavy & light                     # confident H + L (the bar)
            paired_vc = heavy & light & hv & lv          # full V+C pairing
            both_light = (sub.IGK_C >= thr) & (sub.IGL_C >= thr)
            print(f" thr>={thr} reads:")
            print(f"   has IGH_C:            {int(heavy.sum()):>6}")
            print(f"   isotyped (H_C+L_C):   {int(isotyped.sum()):>6}")
            print(f"   paired V+C (H+L):     {int(paired_vc.sum()):>6}")
            print(f"   kappa+lambda both:    {int(both_light.sum()):>6}   "
                  f"(allelic-exclusion-violation candidates)")
            if isotyped.sum():
                med = sub.loc[isotyped, "ig_total"].median()
                nk = int((sub.loc[isotyped, "IGK_C"] >
                          sub.loc[isotyped, "IGL_C"]).sum())
                nl = int(isotyped.sum()) - nk
                m = sub.loc[isotyped, "mem"].sum()
                s = sub.loc[isotyped, "sec"].sum()
                print(f"   median IG reads/cell (isotyped): {med:.0f}")
                print(f"   kappa:lambda among isotyped:     {nk}:{nl}")
                print(f"   IGHM membrane:secreted reads:    "
                      f"{int(m)}:{int(s)}")

    report(df, "ALL barcodes with IG signal")
    report(df[df.in_whitelist], "WHITELIST (short-read called) cells only")
    return df


def _strip_suffix(col):
    return col.astype(str).str.replace(r"-1$", "", regex=True)


def join_labels(df):
    """Step 2: join the per-cell table to the anchor-transfer labels on barcode."""
    lab = pd.read_csv(LABELS, sep="\t")
    df = df.copy()
    df["barcode_key"] = _strip_suffix(df["barcode"])
    lab["barcode_key"] = _strip_suffix(lab["barcode"])
    m = df.merge(lab, on="barcode_key", how="inner", suffixes=("", "_lab"))
    m = m.drop_duplicates(subset="barcode_key")
    m["is_b"] = m["seurat_l2"].isin(B_STATES)
    m.to_csv(OUT / "bcr_joined.tsv", sep="\t", index=False)
    print(f"\njoined {len(m):,} barcodes ({int(m['is_b'].sum()):,} labelled B cells)")
    return m


def _purity_block(sub):
    entered = sub[(sub["IGH_C"] >= PURITY_MIN_READS) & (sub["light_C"] >= PURITY_MIN_READS)]
    minor = entered[["IGK_C", "IGL_C"]].min(axis=1)
    minor_frac = minor / entered["light_C"].where(entered["light_C"] > 0, 1.0)
    violation = (minor >= VIOLATION_MINOR_READS) & (minor_frac > VIOLATION_MINOR_FRAC)
    n, v = int(len(entered)), int(violation.sum())
    return n, v, (100 * v / n if n else float("nan"))


def light_chain_exclusion(joined):
    """Step 3: both light-chain classes in one barcode, B cells and non-B cells."""
    print("\n===== light-chain exclusion =====")
    for tag, sub in (("B cells", joined[joined["is_b"]]), ("non-B cells", joined[~joined["is_b"]])):
        n, v, pct = _purity_block(sub)
        print(f"  {tag}: {n:,} entered, {v:,} carry both classes ({pct:.1f}%)")


def igm_form(joined):
    """Step 4: the IgM form with more reads, per B state."""
    print("\n===== IgM form (membrane vs secreted) per B state =====")
    b = joined[joined["is_b"]]
    for state in B_STATES:
        sub = b[b["seurat_l2"] == state]
        entered = (sub["mem"] + sub["sec"]) >= IGM_MIN_SUM
        mem, sec = sub.loc[entered, "mem"], sub.loc[entered, "sec"]
        n = int(entered.sum())
        n_mem = int((mem > sec).sum())
        n_sec = n - n_mem
        print(f"  {state}: {n:,} entered; membrane {n_mem:,} "
              f"({100 * n_mem / n if n else float('nan'):.1f}%), secreted {n_sec:,}")


def main():
    df = per_cell_table()
    joined = join_labels(df)
    light_chain_exclusion(joined)
    igm_form(joined)


if __name__ == "__main__":
    main()
