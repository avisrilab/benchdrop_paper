#!/usr/bin/env python3
"""Rescore Bagpiper's concordance on the transcript set IsoQuant reports, so the comparison is on a
matched annotation basis (Supplementary Table 2).

Usage:  python3 matched_basis.py <OUT_dir>
"""
import gzip
import sys
from pathlib import Path


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def _first(d, pats):
    for p in pats:
        hit = sorted(Path(d).glob(p))
        if hit:
            return hit[0]
    raise FileNotFoundError(f"none of {pats} under {d}")


def _mtx(path, cell_axis, keep_rows):
    """Stream a MatrixMarket triplet file -> {cell_idx: {tx_idx: value}} for kept cells.

    cell_axis says which of the two coordinates is the cell: bagpiper writes cells x tx,
    IsoQuant writes tx x cells. Streaming keeps peak memory at the size of the KEPT
    submatrix, not the full 1.26M-cell table.
    """
    out = {}
    with _open(path) as fh:
        for ln in fh:
            if ln.startswith("%"):
                continue
            break
        for ln in fh:  # first non-comment line was the dims header, consumed above
            a, b, v = ln.split()
            ci, ti = (int(a) - 1, int(b) - 1) if cell_axis == 0 else (int(b) - 1, int(a) - 1)
            if ci in keep_rows:
                row = out.setdefault(ci, {})
                row[ti] = row.get(ti, 0.0) + float(v)
    return out


def load_bagpiper(d, want):
    cells = [ln.rstrip("\n") for ln in _open(_first(d, ["barcodes.tsv.gz", "barcodes.tsv"]))]
    tx = [ln.split("\t")[0].rstrip("\n")
          for ln in _open(_first(d, ["features.tsv.gz", "features.tsv"]))]
    idx = {c: i for i, c in enumerate(cells)}
    keep = {idx[c]: c for c in want if c in idx}
    per = _mtx(_first(d, ["matrix.mtx.gz", "matrix.mtx"]), 0, set(keep))
    return tx, {keep[i]: v for i, v in per.items()}


def load_isoquant(d, want):
    stem = "**/*.transcript_grouped_counts"
    cells = [ln.rstrip("\n") for ln in _open(_first(d, [stem + ".barcodes.tsv"]))]
    tx = [ln.split("\t")[0].rstrip("\n") for ln in _open(_first(d, [stem + ".features.tsv"]))]
    idx = {c: i for i, c in enumerate(cells)}
    keep = {idx[c]: c for c in want if c in idx}
    per = _mtx(_first(d, [stem + ".matrix.mtx"]), 1, set(keep))
    return tx, {keep[i]: v for i, v in per.items()}


def load_bambu(d, want):
    """counts_transcript.txt: TXNAME, GENEID, then one column per cell. Dense, 300 cells."""
    f = _first(d, ["counts_transcript.txt", "counts_transcript.tsv"])
    with _open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = [(j, c) for j, c in enumerate(hdr) if c in want]
        tx, per = [], {c: {} for _, c in cols}
        for ti, ln in enumerate(fh):
            p = ln.rstrip("\n").split("\t")
            tx.append(p[0])
            for j, c in cols:
                v = float(p[j])
                if v:
                    per[c][ti] = v
    return tx, per


LOADERS = {"bagpiper": load_bagpiper, "isoquant": load_isoquant, "bambu": load_bambu}

# Which IsoQuant arm supplies the basis. `isoquant` is its unique_only DEFAULT, whose 25,346-tx
# output annotation is itself a product of that mode (only transcripts that got a unique read got a
# row), so a matched basis built from it is contaminated. Once the ambiguous-aware arm exists, pass
# `isoquant_wa` and rescore against that instead. Set from argv[2]; see Table1/public_benchmark.sh step 3.
IQ_DIR = "isoquant"
DIR_OF = {"bagpiper": "bagpiper", "isoquant": IQ_DIR, "bambu": "bambu"}


def cells_of(d, tool):
    """Cell barcodes each tool emitted, read without loading any counts."""
    if tool == "bambu":
        with _open(_first(Path(d) / "bambu", ["counts_transcript.txt"])) as fh:
            return set(fh.readline().rstrip("\n").split("\t")[2:])
    if tool == "bagpiper":
        return {ln.rstrip("\n") for ln in
                _open(_first(Path(d) / "bagpiper", ["barcodes.tsv.gz", "barcodes.tsv"]))}
    return {ln.rstrip("\n") for ln in
            _open(_first(Path(d) / IQ_DIR, ["**/*.transcript_grouped_counts.barcodes.tsv"]))}


def ranks(xs):
    """Average ranks, so ties are handled the way scipy's spearmanr does."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va <= 0 or vb <= 0:
        return float("nan")
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    return cov / (va ** 0.5 * vb ** 0.5)


def spearman(a, b):
    return pearson(ranks(a), ranks(b))


def pseudobulk(tx, per, cells):
    """Per-transcript sum over a fixed cell set, as {tx_name: total}."""
    pb = {}
    for c in cells:
        for ti, v in per.get(c, {}).items():
            pb[tx[ti]] = pb.get(tx[ti], 0.0) + v
    return pb


def pair(pa_tx, pa, pb_tx, pb, basis=None):
    """Reproduce metrics.concordance_common's pseudobulk legs on an optional basis."""
    if basis is None:
        bset = set(pb_tx)                       # hoisted: rebuilding it per element is O(n*m)
        shared = [t for t in pa_tx if t in bset]
    else:
        shared = list(basis)
    va = [pa.get(t, 0.0) for t in shared]
    vb = [pb.get(t, 0.0) for t in shared]
    expressed = [i for i in range(len(shared)) if va[i] > 0 or vb[i] > 0]
    codet = [i for i in expressed if va[i] > 0 and vb[i] > 0]
    rho = spearman([va[i] for i in codet], [vb[i] for i in codet]) if len(codet) >= 50 else float("nan")
    return dict(shared_tx=len(shared), expressed_tx=len(expressed),
                codetected_tx=len(codet), pb_spearman_codet=rho)


def _key(a, b):
    """Pair key as metrics.py writes it: the internal 'isoquant' slot carries the arm's dir name,
    so the emitted TSV lines up with metrics_common.tsv and table.py."""
    return "|".join(IQ_DIR if t == "isoquant" else t for t in (a, b))


def main():
    global IQ_DIR
    d = Path(sys.argv[1])
    if len(sys.argv) > 2:
        IQ_DIR = sys.argv[2]
        DIR_OF["isoquant"] = IQ_DIR
    print(f"isoquant arm: {IQ_DIR}", file=sys.stderr)
    common = sorted(set.intersection(*[cells_of(d, t) for t in LOADERS]))
    print(f"common cells: {len(common)}", file=sys.stderr)

    tx, per, pb = {}, {}, {}
    for t in LOADERS:
        print(f"loading {t} ...", file=sys.stderr)
        tx[t], per[t] = LOADERS[t](d / DIR_OF[t], set(common))
        pb[t] = pseudobulk(tx[t], per[t], common)

    iq_annot = set(tx["isoquant"])
    iq_bambu = pair(tx["isoquant"], pb["isoquant"], tx["bambu"], pb["bambu"])
    iq_bambu_codet = {t for t in tx["isoquant"]
                      if pb["isoquant"].get(t, 0) > 0 and pb["bambu"].get(t, 0) > 0}

    print("\n== BASELINE: each pair on its own shared annotation ==")
    base = {}
    for a, b in [("bagpiper", "bambu"), ("bagpiper", "isoquant"), ("isoquant", "bambu")]:
        base[(a, b)] = pair(tx[a], pb[a], tx[b], pb[b])
        r = base[(a, b)]
        print(f"{a:>9}|{b:<9} shared={r['shared_tx']:>7,} expressed={r['expressed_tx']:>7,} "
              f"codetected={r['codetected_tx']:>7,} pb_spearman_codet={r['pb_spearman_codet']:.4f}")

    print("\n== MATCHED-A: all pairs on IsoQuant's annotation "
          f"({len(iq_annot):,} tx) ==")
    matched = {}
    for a, b in [("bagpiper", "bambu"), ("bagpiper", "isoquant"), ("isoquant", "bambu")]:
        bset = set(tx[b])
        basis = [t for t in tx[a] if t in iq_annot and t in bset]
        matched[(a, b)] = r = pair(tx[a], pb[a], tx[b], pb[b], basis=basis)
        print(f"{a:>9}|{b:<9} shared={r['shared_tx']:>7,} expressed={r['expressed_tx']:>7,} "
              f"codetected={r['codetected_tx']:>7,} pb_spearman_codet={r['pb_spearman_codet']:.4f}")

    print(f"\n== MATCHED-B: Bagpiper|bambu on the exact {len(iq_bambu_codet):,} tx "
          f"IsoQuant and bambu co-detect ==")
    r_b = pair(tx["bagpiper"], pb["bagpiper"], tx["bambu"], pb["bambu"],
               basis=sorted(iq_bambu_codet))
    print(f"{'bagpiper':>9}|{'bambu':<9} shared={r_b['shared_tx']:>7,} "
          f"expressed={r_b['expressed_tx']:>7,} codetected={r_b['codetected_tx']:>7,} "
          f"pb_spearman_codet={r_b['pb_spearman_codet']:.4f}")
    print(f"{'isoquant':>9}|{'bambu':<9} (reference, same basis) "
          f"pb_spearman_codet={iq_bambu['pb_spearman_codet']:.4f}")

    # Written for table.py, so the table cannot drift from what was scored here.
    with open(d / "matched_basis.tsv", "w") as fh:
        fh.write("pair\tbasis\tcodetected_tx\tpb_spearman_codet\n")
        for (a, b), r in matched.items():
            fh.write(f"{_key(a, b)}\tisoquant_annotation\t{r['codetected_tx']}\t"
                     f"{r['pb_spearman_codet']:.6f}\n")
        fh.write(f"bagpiper|bambu\tisoquant_bambu_codetected\t{r_b['codetected_tx']}\t"
                 f"{r_b['pb_spearman_codet']:.6f}\n")
    print(f"wrote {d / 'matched_basis.tsv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
