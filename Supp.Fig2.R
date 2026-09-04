# Supp. Figure 2: panels S2A, S2D, S2E. Inputs: read.lens.k562.txt, LR.K562.S3.rds, alevin/, kb/
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))

#### S2A: long-read length distribution ####
save.panel(read.length.density("read.lens.k562.txt", "SuppFig2A"), "SuppFig2A", width = 8, height = 5)

#### S2D-S2E: BenchDrop-seq vs short-read quantifiers, same cells ####
lr.k562 <- readRDS("LR.K562.S3.rds")
lr.pb <- pseudobulk(lr.k562)

read.quant <- function(mtx, cells, features) {
  m <- t(readMM(mtx))
  dimnames(m) <- list(readLines(features), readLines(cells))
  m[, intersect(colnames(m), colnames(lr.k562))]
}

alevin <- read.quant("alevin/quants_mat.mtx.gz", "alevin/quants_mat_rows.txt", "alevin/quants_mat_cols.txt")
save.panel(cor.hex(rowSums(alevin), lr.pb, "Salmon-Alevin", "BenchDrop-seq", "SuppFig2D"), "SuppFig2D")

kb <- read.quant("kb/cells_x_genes.mtx", "kb/cells_x_genes.barcodes.txt", "kb/cells_x_genes.genes.names.txt")
save.panel(cor.hex(rowSums(kb), lr.pb, "Kallisto-Bustools", "BenchDrop-seq", "SuppFig2E"), "SuppFig2E")
