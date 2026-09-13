"""Export the Azimuth PBMC reference (h5seurat) into plain arrays that transfer_labels.R reads natively.
"""
from __future__ import annotations

import os

import h5py
import numpy as np

REF = os.path.expandvars("$BENCHDROP_FEED/azimuth_pbmc/pbmc_multimodal.h5seurat")
OUT = os.path.expandvars("$BENCHDROP_FEED/azimuth_pbmc/export_r")
CHUNK = 4000


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    f = h5py.File(REF, "r")
    sct = f["assays"]["SCT"]
    feats = [x.decode() for x in sct["features"][:]]
    spca_feats = [x.decode() for x in f["reductions"]["spca"]["features"][:]]
    keep = np.array([feats.index(g) for g in spca_feats])
    print(f"reference genes {len(feats)}, keeping {len(keep)} spca features")

    d = sct["data"]
    data, indices, indptr = d["data"], d["indices"], d["indptr"][:]
    ncell = len(indptr) - 1
    print("cells:", ncell, "| nnz total:", data.shape[0])

    sel = np.zeros(len(feats), dtype=bool)
    sel[keep] = True
    newrow = -np.ones(len(feats), dtype=np.int32)
    newrow[keep] = np.arange(len(keep))

    out_i, out_d = [], []
    out_p = np.zeros(ncell + 1, dtype=np.int64)
    for s in range(0, ncell, CHUNK):
        e = min(s + CHUNK, ncell)
        lo, hi = indptr[s], indptr[e]
        idx, val = indices[lo:hi], data[lo:hi]
        m = sel[idx]
        cuts = indptr[s:e + 1] - lo
        for j in range(e - s):
            out_p[s + j + 1] = out_p[s + j] + int(m[cuts[j]:cuts[j + 1]].sum())
        out_i.append(newrow[idx[m]])
        out_d.append(val[m].astype(np.float32))

    ind = np.concatenate(out_i)
    dat = np.concatenate(out_d)
    print(f"kept nnz: {len(dat)} ({100 * len(dat) / data.shape[0]:.1f}% of total)")

    # Row indices must be INCREASING WITHIN EACH COLUMN or R's dgCMatrix validity check rejects
    # the object ("'i' slot is not increasing within columns"). They are not, because the SPCA
    # feature order is not monotonic in the reference's original gene order, so remapping through
    # `newrow` permutes them. scipy sorts them in place; doing it here rather than in R keeps the
    # binary files a directly-loadable CSC.
    import scipy.sparse as sp
    m = sp.csc_matrix((dat, ind, out_p), shape=(len(keep), ncell))
    m.sort_indices()
    assert m.has_sorted_indices
    ind, dat, out_p = m.indices, m.data, m.indptr

    ind.astype(np.int32).tofile(f"{OUT}/indices.i32")
    dat.astype(np.float32).tofile(f"{OUT}/data.f32")
    out_p.astype(np.int64).tofile(f"{OUT}/indptr.i64")

    cells = [x.decode() for x in f["cell.names"][:]]
    open(f"{OUT}/genes.txt", "w").write("\n".join(spca_feats) + "\n")
    open(f"{OUT}/cells.txt", "w").write("\n".join(cells) + "\n")

    md = f["meta.data"]

    def col(name):
        g = md[name]
        if isinstance(g, h5py.Group):
            lv = [x.decode() for x in g["levels"][:]]
            return [lv[i - 1] for i in g["values"][:]]
        return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in g[:]]

    with open(f"{OUT}/meta.tsv", "w") as fh:
        fh.write("cell\tl1\tl2\n")
        for c, a, b in zip(cells, col("celltype.l1"), col("celltype.l2")):
            fh.write(f"{c}\t{a}\t{b}\n")

    # Transposed on write: h5seurat keeps R's column-major layout, so these read as
    # (50 PCs, cells) and (50 PCs, features) and R wants observations in rows.
    np.savetxt(f"{OUT}/spca_embeddings.tsv",
               f["reductions"]["spca"]["cell.embeddings"][:].T, delimiter="\t", fmt="%.6g")
    np.savetxt(f"{OUT}/spca_loadings.tsv",
               f["reductions"]["spca"]["feature.loadings"][:].T, delimiter="\t", fmt="%.6g")
    print(f"wrote export to {OUT}: {len(spca_feats)} genes x {ncell} cells")


if __name__ == "__main__":
    main()
