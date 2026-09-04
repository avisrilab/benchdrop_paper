# Supp. Figure 5: panels S5A-S5G. Inputs: SR/LR.PBMC.S3.rds
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))

sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

lr.pbmc$celltype.filtered <- group.celltypes(sr.pbmc$predicted.celltype.l2, fine = TRUE)[match(colnames(lr.pbmc), colnames(sr.pbmc))]
cells <- colnames(lr.pbmc)[!is.na(lr.pbmc$celltype.filtered)]
DefaultAssay(lr.pbmc) <- "transcript"
Idents(lr.pbmc) <- "celltype.filtered"
message(sprintf("SuppFig5A: %d of %d cells carry a plotted cell type", length(cells), ncol(lr.pbmc)))

#### S5A: transcript-level UMAP ####
p <- DimPlot(lr.pbmc, reduction = "transcript.umap", group.by = "celltype.filtered", cells = cells, pt.size = 1) +
  NoLegend() + xlab("UMAP_1") + ylab("UMAP_2") + theme(plot.title = element_blank())
save.panel(LabelClusters(p, id = "celltype.filtered", fontface = "bold", size = 5.5, repel = TRUE, color = "black"),
           "SuppFig5A", width = 9, height = 8)

#### S5B-S5G: transcript feature plots ####
features <- list(
  B = list("IL7R", "ENST00000303115", "CD4 T Cells"),
  C = list("CD8A", "ENST00000352580", "CD8 T Cells"),
  D = list("NKG7", "ENST00000221978", "NK"),
  E = list("SERPINA1", "ENST00000636712", c("CD14 Mono", "CD16 Mono")),
  F = list("CST3", "ENST00000376925", c("DC", "pDC")),
  G = list("CD79A", "ENST00000221972", "B Cells")
)
for (panel in names(features)) {
  f <- features[[panel]]
  p <- FeaturePlot(lr.pbmc, reduction = "transcript.umap", features = f[[2]], cells = cells, pt.size = 0.75) +
    xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(paste(f[[2]], "/", f[[1]])) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  p <- LabelClusters(p, id = "ident", clusters = f[[3]], size = 7, repel = TRUE, fontface = "bold", color = "black")
  save.panel(p, paste0("SuppFig5", panel), width = 8, height = 7)
}
