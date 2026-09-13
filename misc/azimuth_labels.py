"""Cell labels for the PBMC pseudotime vignette, transferred from the Azimuth PBMC reference
(Methods 2.8).
"""
from __future__ import annotations
import os

from pathlib import Path

import h5py
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

REF = Path(os.path.expandvars("$BENCHDROP_FEED/azimuth_pbmc/pbmc_multimodal.h5seurat"))  # 2026-08-14: references/ layer removed
SR = Path(os.path.expandvars("$BENCHDROP_FEED/pbmc/count_sr/Gene/filtered"))
OUT = Path(os.path.expandvars("$BENCHDROP_RESULTS/pseudotime"))

K_NEIGHBOURS = 50
BARS = {"min_l1_agreement": 0.85, "min_subtype_cells": 30, "min_median_vote": 0.50}


def _decode(arr):
    return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]


def read_reference(path: Path = REF) -> dict:
    """Pull labels and the SPCA loadings/embeddings from the h5seurat.

    ORIENTATION, verified against the file 2026-08-13 rather than assumed: h5seurat keeps R's
    column-major layout, so `cell.embeddings` reads as (50 PCs, 161,764 cells) and
    `feature.loadings` as (50 PCs, 5,000 features). Both are transposed here to the
    rows-are-observations convention the rest of this module uses.
    """
    with h5py.File(path, "r") as f:
        meta = f["meta.data"]

        def col(name):
            d = meta[name]
            if isinstance(d, h5py.Group):          # factor: values + levels
                levels = _decode(d["levels"][:])
                return np.array([levels[i - 1] for i in d["values"][:]], dtype=object)
            return np.array(_decode(d[:]), dtype=object)

        red = f["reductions"]["spca"]
        return {"l1": col("celltype.l1"), "l2": col("celltype.l2"),
                "ref_emb": red["cell.embeddings"][:].T,        # -> (cells, 50)
                "loadings": red["feature.loadings"][:].T,      # -> (5000, 50)
                "load_features": _decode(red["features"][:])}


def load_query():
    """Matched PBMC short-read gene matrix (same barcodes as the long-read cells)."""
    X = sio.mmread(SR / "matrix.mtx").tocsc().T.tocsr()      # genes x cells -> cells x genes
    genes = [ln.split("\t")[1] for ln in (SR / "features.tsv").read_text().splitlines()]
    bcs = (SR / "barcodes.tsv").read_text().split()
    return X, genes, bcs


def transfer(ref: dict, X, genes, k: int = K_NEIGHBOURS):
    """Project the query into the reference SPCA space and vote over kNN labels."""
    load_feats = ref["load_features"]
    gi = {g: j for j, g in enumerate(genes)}
    shared = [(i, gi[g]) for i, g in enumerate(load_feats) if g in gi]
    ri = np.array([i for i, _ in shared])
    qi = np.array([j for _, j in shared])

    Q = X[:, qi].astype(np.float64)
    lib = np.asarray(Q.sum(1)).ravel()
    lib[lib == 0] = 1
    Q = sp.diags(1e4 / lib) @ Q
    Q = Q.toarray()
    np.log1p(Q, out=Q)
    Q -= Q.mean(0)
    sd = Q.std(0)
    sd[sd == 0] = 1
    Q /= sd
    np.clip(Q, -10, 10, out=Q)

    emb_q = Q @ ref["loadings"][ri]                # (cells, 5000_shared) @ (5000_shared, 50)
    emb_r = ref["ref_emb"]                          # (161764, 50)

    # cosine kNN in the shared space (Azimuth scores neighbours by distance in spca space too)
    def norm(a):
        n = np.linalg.norm(a, axis=1, keepdims=True)
        n[n == 0] = 1
        return a / n

    qn, rn = norm(emb_q), norm(emb_r)
    labs_l1, labs_l2, votes = [], [], []
    step = 512
    for s in range(0, qn.shape[0], step):
        sim = qn[s:s + step] @ rn.T
        idx = np.argpartition(-sim, k, axis=1)[:, :k]
        for row, nb in enumerate(idx):
            l2 = ref["l2"][nb]
            vals, cnt = np.unique(l2, return_counts=True)
            win = int(np.argmax(cnt))
            labs_l2.append(vals[win])
            votes.append(cnt[win] / k)
            l1 = ref["l1"][nb]
            v1, c1 = np.unique(l1, return_counts=True)
            labs_l1.append(v1[int(np.argmax(c1))])
    return np.array(labs_l1, dtype=object), np.array(labs_l2, dtype=object), np.array(votes)


# Azimuth l1 -> the vignette's coarse vocabulary, for the A1 concordance bar.
L1_TO_COARSE = {"CD4 T": "t_cell", "CD8 T": "t_cell", "other T": "t_cell",
                "B": "b_cell", "NK": "nk", "Mono": "mono"}

# Azimuth l2 -> the five subtypes the canonical orderings need.
L2_TO_SUBTYPE = {
    "CD4 Naive": "t_naive", "CD8 Naive": "t_naive",
    "CD4 TCM": "t_memory", "CD4 TEM": "t_memory", "CD8 TCM": "t_memory", "CD8 TEM": "t_memory",
    "CD4 CTL": "t_memory", "Treg": "t_memory", "MAIT": "t_memory",
    "B naive": "b_naive", "B intermediate": "b_memory", "B memory": "b_memory",
    "Plasmablast": "plasmablast",
}
