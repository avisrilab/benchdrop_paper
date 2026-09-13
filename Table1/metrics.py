#!/usr/bin/env python3
"""Benchmark metrics: Bagpiper, IsoQuant and bambu on identical inputs (Table 1, Methods 2.9).

Usage:  python3 metrics.py <OUT_dir>
"""
import sys, re, glob, gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats, io as sio


def _first(d, pats):
    for p in pats:
        hit = sorted(glob.glob(str(Path(d) / p)))
        if hit:
            return hit[0]
    raise FileNotFoundError(f"none of {pats} under {d}")


def _lines(f):
    op = gzip.open if str(f).endswith(".gz") else open
    return [ln.rstrip("\n") for ln in op(f, "rt")]


def _orient(m, cells, tx):
    m = sparse.csr_matrix(m)
    if m.shape == (len(tx), len(cells)) and m.shape[0] != m.shape[1]:
        m = m.T.tocsr()
    assert m.shape == (len(cells), len(tx)), f"shape {m.shape} != {(len(cells), len(tx))}"
    return dict(mat=m, cells=list(cells), tx=list(tx))


# ---------- LOADERS: return {mat: cells x tx CSR, cells: [...], tx: [...]} ----------
# Confirm the exact output filename per installed tool version; only these globs change.

def load_bagpiper(d):
    m = sio.mmread(_first(d, ["quant.mtx", "matrix.mtx", "matrix.mtx.gz"]))
    cells = _lines(_first(d, ["barcodes.tsv", "barcodes.tsv.gz"]))
    tx = [ln.split("\t")[0] for ln in
          _lines(_first(d, ["features.tsv.gz", "features.tsv", "transcripts.tsv"]))]
    return _orient(m, cells, tx)


def load_isoquant(d):
    # IsoQuant --read_group tag:CB writes an MTX triplet under <output>/OUT/:
    #   OUT.transcript_grouped_counts.{matrix.mtx, barcodes.tsv, features.tsv}. The '.linear.tsv' is
    #   long-format (feature,group,count), not a matrix, so read the MTX. The leading '.' before
    #   'transcript' excludes the 'discovered_' novel-model variant (we run annotation-guided).
    m = sio.mmread(_first(d, ["**/*.transcript_grouped_counts.matrix.mtx",
                              "*.transcript_grouped_counts.matrix.mtx"]))
    cells = _lines(_first(d, ["**/*.transcript_grouped_counts.barcodes.tsv",
                              "*.transcript_grouped_counts.barcodes.tsv"]))
    tx = [ln.split("\t")[0] for ln in
          _lines(_first(d, ["**/*.transcript_grouped_counts.features.tsv",
                            "*.transcript_grouped_counts.features.tsv"]))]
    return _orient(sparse.csr_matrix(m), cells, tx)


def load_bambu(d):
    # bambu writeBambuOutput: counts_transcript.txt -- TXNAME, GENEID, then one column per sample
    #   (= a cell; the per-cell BAM basename is the CB). Drop the two annotation columns.
    f = _first(d, ["counts_transcript.txt", "counts_transcript.tsv"])
    df = pd.read_csv(f, sep="\t")
    df = df.set_index(df.columns[0])                 # TXNAME
    df = df.drop(columns=[c for c in ("GENEID",) if c in df.columns]).select_dtypes("number")
    return _orient(sparse.csr_matrix(df.values), list(df.columns), list(df.index))


# FLAMES + wf-single-cell are CITED, not run (see Table1/public_benchmark.sh header), so no loader.
# `isoquant_wa` is the SAME tool and loader under ambiguous-aware quantification instead of its
# `unique_only` default (see Table1/public_benchmark.sh step 3). Keyed separately because the dir and
# the time.<tool> log are separate; main() skips it cleanly when that arm has not been run.
LOADERS = {"bagpiper": load_bagpiper, "isoquant": load_isoquant,
           "isoquant_wa": load_isoquant, "bambu": load_bambu}


# ---------- metrics ----------

def recovery(M):
    m = M["mat"]
    counts = np.asarray(m.sum(1)).ravel()
    det = np.asarray((m > 0).sum(1)).ravel()
    return dict(cells=len(M["cells"]),
                median_counts_per_cell=float(np.median(counts)) if counts.size else np.nan,
                median_tx_per_cell=float(np.median(det)) if det.size else np.nan,
                transcripts_detected=int((np.asarray((m > 0).sum(0)).ravel() > 0).sum()))


def _sub(M, cells, tx):
    ci = {c: i for i, c in enumerate(M["cells"])}
    ti = {t: i for i, t in enumerate(M["tx"])}
    return M["mat"][[ci[c] for c in cells]][:, [ti[t] for t in tx]].tocsr()


def concordance(A, B):
    bc = set(B["cells"]); bt = set(B["tx"])
    cells = [c for c in A["cells"] if c in bc]
    tx = [t for t in A["tx"] if t in bt]
    out = dict(shared_cells=len(cells), shared_tx=len(tx),
               pb_spearman=np.nan, cell_spearman=np.nan, jaccard=np.nan)
    if len(cells) < 20 or len(tx) < 50:
        return out
    a, b = _sub(A, cells, tx), _sub(B, cells, tx)
    pa, pb = np.asarray(a.sum(0)).ravel(), np.asarray(b.sum(0)).ravel()
    out["pb_spearman"] = float(stats.spearmanr(pa, pb).correlation)   # pseudobulk agreement
    idx = np.arange(len(cells))
    sel = idx if len(cells) <= 500 else np.random.RandomState(0).choice(idx, 500, replace=False)
    ss, jj = [], []
    for i in sel:
        av, bv = a[i].toarray().ravel(), b[i].toarray().ravel()
        if (av > 0).sum() < 3 or (bv > 0).sum() < 3:
            continue
        ss.append(stats.spearmanr(av, bv).correlation)
        sa, sb = set(np.nonzero(av)[0]), set(np.nonzero(bv)[0])
        jj.append(len(sa & sb) / len(sa | sb))
    out["cell_spearman"] = float(np.nanmean(ss)) if ss else np.nan
    out["jaccard"] = float(np.nanmean(jj)) if jj else np.nan
    return out


def recovery_on(M, cells):
    # recovery restricted to a fixed cell subset, so per-arm counts are comparable across tools
    # (the raw recovery() spans each tool's own cell set: raw barcodes vs bambu's called top-N).
    ci = {c: i for i, c in enumerate(M["cells"])}
    idx = [ci[c] for c in cells if c in ci]
    m = M["mat"][idx]
    counts = np.asarray(m.sum(1)).ravel()
    det = np.asarray((m > 0).sum(1)).ravel()
    return dict(cells=len(idx),
                median_counts_per_cell=float(np.median(counts)) if counts.size else np.nan,
                median_tx_per_cell=float(np.median(det)) if det.size else np.nan,
                transcripts_detected=int((np.asarray((m > 0).sum(0)).ravel() > 0).sum()))


def concordance_common(A, B, cells):
    # Pairwise concordance on a FIXED common cell set (passed in), and after dropping transcripts
    # that are zero in BOTH arms across those cells. That removes the two raw-table artifacts:
    # differing cell sets, and Spearman inflated by tens of thousands of tied zero-padded tx.
    bt = set(B["tx"])
    tx = [t for t in A["tx"] if t in bt]
    out = dict(shared_cells=len(cells), shared_tx=len(tx), expressed_tx=0, codetected_tx=0,
               cells_scored=0, pb_spearman=np.nan, pb_spearman_codet=np.nan,
               cell_spearman=np.nan, jaccard=np.nan)
    if len(cells) < 20 or len(tx) < 50:
        return out
    a, b = _sub(A, cells, tx), _sub(B, cells, tx)
    pa, pb = np.asarray(a.sum(0)).ravel(), np.asarray(b.sum(0)).ravel()
    keep = (pa > 0) | (pb > 0)                          # expressed by at least one arm here
    out["expressed_tx"] = int(keep.sum())
    if keep.sum() < 50:
        return out
    a, b, pa, pb = a[:, keep], b[:, keep], pa[keep], pb[keep]
    out["pb_spearman"] = float(stats.spearmanr(pa, pb).correlation)   # union: detection + quant
    both = (pa > 0) & (pb > 0)                          # co-detected: quantitative agreement alone
    out["codetected_tx"] = int(both.sum())
    if both.sum() >= 50:
        out["pb_spearman_codet"] = float(stats.spearmanr(pa[both], pb[both]).correlation)
    ss, jj = [], []
    for i in range(len(cells)):                         # <=~300 common cells: score them all
        av, bv = a[i].toarray().ravel(), b[i].toarray().ravel()
        if (av > 0).sum() < 3 or (bv > 0).sum() < 3:
            continue
        ss.append(stats.spearmanr(av, bv).correlation)
        sa, sb = set(np.nonzero(av)[0]), set(np.nonzero(bv)[0])
        jj.append(len(sa & sb) / len(sa | sb))
    out["cells_scored"] = len(ss)
    out["cell_spearman"] = float(np.nanmean(ss)) if ss else np.nan
    out["jaccard"] = float(np.nanmean(jj)) if jj else np.nan
    return out


def parse_time(f):
    # macOS `/usr/bin/time -l`: "<sec> real" and "<bytes> maximum resident set size"
    if not Path(f).exists():
        return dict(wall_s=np.nan, peak_rss_gb=np.nan)
    txt = Path(f).read_text()
    w = re.search(r"([\d.]+)\s+real", txt)
    r = re.search(r"(\d+)\s+maximum resident set size", txt)
    return dict(wall_s=float(w.group(1)) if w else np.nan,
                peak_rss_gb=int(r.group(1)) / 1e9 if r else np.nan)


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: metrics.py <OUT_dir>")
    out = Path(sys.argv[1])
    loaded = {}
    for tool, fn in LOADERS.items():
        try:
            loaded[tool] = fn(out / tool)
            print(f"[ok]   {tool}: {loaded[tool]['mat'].shape[0]} cells x "
                  f"{loaded[tool]['mat'].shape[1]} tx", file=sys.stderr)
        except Exception as e:
            print(f"[skip] {tool}: {e}", file=sys.stderr)
    rows = []
    for tool, M in loaded.items():
        for k, v in {**recovery(M), **parse_time(out / f"time.{tool}")}.items():
            rows.append((tool, k, v))
    tl = list(loaded)
    for i in range(len(tl)):
        for j in range(i + 1, len(tl)):
            for k, v in concordance(loaded[tl[i]], loaded[tl[j]]).items():
                rows.append((f"{tl[i]}|{tl[j]}", k, v))
    # shared alignment wall time (pipelines, wall only): per-tool total = align + its quantifier
    for al in ("txome_align", "genome_align"):
        tf = out / f"time.{al}"
        if tf.exists():
            rows.append((al, "wall_s", parse_time(tf)["wall_s"]))
    df = pd.DataFrame(rows, columns=["tool", "metric", "value"])
    df.to_csv(out / "metrics.tsv", sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {out / 'metrics.tsv'}")

    # ---- common-basis: every arm + pair scored on the SAME cells, expressed tx only ----
    # The 3-way cell intersection = bambu's called cells (the only cell-called arm), so this is the
    # apples-to-apples view. Raw metrics.tsv above keeps the full picture (bagpiper's permissiveness,
    # the noise-barcode tail); metrics_common.tsv is the comparable one for the writeup.
    if len(loaded) >= 2:
        inter = set.intersection(*[set(M["cells"]) for M in loaded.values()])
        common = [c for c in loaded[tl[0]]["cells"] if c in inter]
        crows = [("common", "cells", float(len(common)))]
        for tool, M in loaded.items():
            for k, v in recovery_on(M, common).items():
                crows.append((f"{tool}|common", k, v))
        for i in range(len(tl)):
            for j in range(i + 1, len(tl)):
                for k, v in concordance_common(loaded[tl[i]], loaded[tl[j]], common).items():
                    crows.append((f"{tl[i]}|{tl[j]}", k, v))
        cdf = pd.DataFrame(crows, columns=["tool", "metric", "value"])
        cdf.to_csv(out / "metrics_common.tsv", sep="\t", index=False)
        print("\n== common-basis (cells in ALL arms; expressed tx only) ==")
        print(cdf.to_string(index=False))
        print(f"wrote {out / 'metrics_common.tsv'}")


if __name__ == "__main__":
    main()
