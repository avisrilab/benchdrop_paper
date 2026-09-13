# Anchor-based cell-type label transfer from the Azimuth PBMC reference onto the PBMC short-read
# matrix, in R with Seurat (Methods 2.8). The reference is exported once by export_reference.py.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

# SeuratDisk's LoadH5Seurat fails on this reference under SeuratObject 5.x, so the reference is
# exported once from Python (export_reference.py) to binary CSC arrays plus TSVs, covering the
# 5,000 SPCA features, and the Seurat object is rebuilt here.
REFDIR <- file.path(Sys.getenv("BENCHDROP_FEED"), "azimuth_pbmc/export_r")
SR     <- file.path(Sys.getenv("BENCHDROP_FEED"), "pbmc/count_sr/Gene/filtered")
OUT    <- file.path(Sys.getenv("BENCHDROP_RESULTS"), "pseudotime")
NDIMS  <- 50

message("[1/5] rebuilding reference from the exported CSC arrays ...")
rgenes <- readLines(file.path(REFDIR, "genes.txt"))
rcells <- readLines(file.path(REFDIR, "cells.txt"))
ip <- readBin(file.path(REFDIR, "indptr.i64"), "integer", n = length(rcells) + 1L, size = 8)
nnz <- ip[length(ip)]
ii <- readBin(file.path(REFDIR, "indices.i32"), "integer", n = nnz, size = 4)
xx <- readBin(file.path(REFDIR, "data.f32"), "numeric", n = nnz, size = 4)
refm <- new("dgCMatrix", i = ii, p = as.integer(ip), x = as.numeric(xx),
            Dim = c(length(rgenes), length(rcells)),
            Dimnames = list(rgenes, rcells))
meta <- read.delim(file.path(REFDIR, "meta.tsv"), stringsAsFactors = FALSE, row.names = 1)
ref <- CreateSeuratObject(counts = refm, assay = "SCT", meta.data = meta)
ref <- SetAssayData(ref, layer = "data", new.data = refm, assay = "SCT")
VariableFeatures(ref) <- rgenes
emb  <- as.matrix(read.delim(file.path(REFDIR, "spca_embeddings.tsv"), header = FALSE))
load <- as.matrix(read.delim(file.path(REFDIR, "spca_loadings.tsv"),   header = FALSE))
rownames(emb) <- rcells; colnames(emb) <- paste0("SPCA_", seq_len(ncol(emb)))
rownames(load) <- rgenes; colnames(load) <- colnames(emb)
ref[["spca"]] <- CreateDimReducObject(embeddings = emb, loadings = load,
                                      key = "SPCA_", assay = "SCT")
message("      reference: ", ncol(ref), " cells, ", nrow(ref), " features, spca ", ncol(emb), " dims")

message("[2/5] loading query (matched PBMC short read) ...")
m <- readMM(file.path(SR, "matrix.mtx"))
genes <- read.delim(file.path(SR, "features.tsv"), header = FALSE)[[2]]
bcs   <- readLines(file.path(SR, "barcodes.tsv"))
rownames(m) <- make.unique(as.character(genes)); colnames(m) <- bcs
q <- CreateSeuratObject(counts = m, project = "PBMC")
# LogNormalize, not SCTransform: FindTransferAnchors(normalization.method = "SCT") demands the
# REFERENCE carry an SCTModel.list, which the rebuilt assay cannot (the model lives in the
# h5seurat's SCTModel.list and reconstructing an SCTAssay faithfully is a bigger lift than the
# transfer itself). The reference `data` exported here is already SCT-log-normalised, so the
# query is put on a comparable log scale and the anchor path runs in LogNormalize mode.
q <- NormalizeData(q, verbose = FALSE)
q <- FindVariableFeatures(q, verbose = FALSE)
message("      query: ", ncol(q), " cells")

message("[3/5] finding transfer anchors (spca, ", NDIMS, " dims) ...")
anchors <- FindTransferAnchors(
  reference = ref, query = q,
  normalization.method = "LogNormalize",
  reference.reduction = "spca",
  dims = 1:NDIMS, verbose = FALSE
)
message("      anchors: ", nrow(slot(anchors, "anchors")))

message("[4/5] transferring celltype.l1 / celltype.l2 ...")
q <- TransferData(anchorset = anchors, reference = ref, query = q,
                  # the rebuilt reference's metadata columns are l1/l2 (export_reference.py's
                  # meta.tsv header), not the h5seurat's celltype.l1/celltype.l2
                  refdata = list(azl1 = "l1", azl2 = "l2"),
                  dims = 1:NDIMS, verbose = FALSE)

res <- data.frame(
  barcode    = colnames(q),
  seurat_l1  = as.character(q$predicted.azl1),
  seurat_l1_score = as.numeric(q$predicted.azl1.score),
  seurat_l2  = as.character(q$predicted.azl2),
  seurat_l2_score = as.numeric(q$predicted.azl2.score),
  stringsAsFactors = FALSE
)
message("[5/5] writing ", file.path(OUT, "seurat_anchor_labels.tsv"))
write.table(res, file.path(OUT, "seurat_anchor_labels.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("TRANSFER-DONE\n")
