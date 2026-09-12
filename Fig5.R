# Figure 4: panels 4A-4F. Inputs: SR/LR.PBMC.S3.rds, t2gnames.txt
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))
library(patchwork)
suppressPackageStartupMessages({
  library(Signac)
  library(EnsDb.Hsapiens.v86)
})

genes <- c(A = "IL7R", B = "CD8A", C = "NKG7", D = "SERPINA1", E = "CST3", F = "CD79A")
celltypes <- c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK")

sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")
t2g <- read.t2g("t2gnames.txt")

sr.pbmc <- subset(sr.pbmc, subset = predicted.celltype.l1 %in% celltypes)
lr.pbmc <- subset(lr.pbmc, cells = intersect(colnames(lr.pbmc), colnames(sr.pbmc)))
lr.pbmc$sr.celltype <- sr.pbmc$predicted.celltype.l1[match(colnames(lr.pbmc), colnames(sr.pbmc))]
sr.expr <- GetAssayData(sr.pbmc, assay = "RNA", layer = "data")
tx.expr <- GetAssayData(lr.pbmc, assay = "transcript", layer = "data")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
stub <- CreateChromatinAssay(
  counts = sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 2),
                        dimnames = list(c("chr1-1-100", "chr2-1-100"), c("c1", "c2"))),
  annotation = annotations, min.cells = 0, min.features = 0)
stub <- CreateSeuratObject(stub, assay = "peaks")

boxplot.by.celltype <- function(df, x, title = NULL) {
  ggplot(df, aes(x = .data[[x]], y = expression, fill = .data[[x]])) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
          strip.background = element_blank(), strip.text = element_text(face = "bold")) +
    labs(x = NULL, y = "Expression", title = title)
}

for (panel in names(genes)) {
  gname <- genes[[panel]]

  gene.plot <- boxplot.by.celltype(
    data.frame(expression = sr.expr[gname, ], celltype = sr.pbmc$predicted.celltype.l1),
    "celltype", title = gname)

  annotation.plot <- AnnotationPlot(stub, region = gname, mode = "transcript")

  tx.ids <- intersect(names(t2g)[t2g == gname], rownames(tx.expr))
  tx.ids <- tx.ids[Matrix::rowMeans(tx.expr[tx.ids, , drop = FALSE]) > 0.005]
  tx.df <- reshape2::melt(as.matrix(tx.expr[tx.ids, , drop = FALSE]),
                          varnames = c("transcript", "cell"), value.name = "expression")
  tx.df$celltype <- lr.pbmc$sr.celltype[match(tx.df$cell, colnames(lr.pbmc))]
  message(sprintf("Fig4%s %s: transcripts shown = %s", panel, gname, paste(tx.ids, collapse = ", ")))
  tx.plot <- boxplot.by.celltype(tx.df, "celltype") + facet_wrap(~ transcript, ncol = length(tx.ids))

  save.panel(gene.plot / annotation.plot / tx.plot, paste0("Fig4", panel), width = 10, height = 9)
}
