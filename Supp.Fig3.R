# Supp. Figure 3: panels S3A-S3E. Inputs: read.lens.pbmc.txt, SR/LR.PBMC.S3.rds, quant.sf, t2gnames.txt
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))
library(patchwork)

#### S3A: long-read length distribution ####
save.panel(read.length.density("read.lens.pbmc.txt", "SuppFig3A"), "SuppFig3A", width = 8, height = 5)

sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")
t2g <- read.t2g("t2gnames.txt")

#### S3B: UMIs, genes and transcripts per cell ####
violin <- function(y, ylab, fill) {
  ggplot(lr.pbmc@meta.data, aes(x = "BenchDrop-seq", y = .data[[y]])) +
    geom_violin(scale = "width", trim = TRUE, fill = fill) +
    theme_minimal(base_size = 25) +
    labs(x = NULL, y = ylab) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 20))
}
save.panel(violin("nCount_RNA", "UMIs per cell", "#F8766D") +
             violin("nFeature_RNA", "Genes per cell", "#00BFC4") +
             violin("nFeature_transcript", "Transcripts per cell", "gray"),
           "SuppFig3B", width = 14, height = 6)

#### S3C-S3E: pseudobulk correlations ####
bulk <- bulk.gene.tpm("quant.sf", t2g)
save.panel(cor.hex(bulk, pseudobulk(lr.pbmc, cpm = TRUE), "Bulk", "BenchDrop-seq", "SuppFig3C"), "SuppFig3C")
save.panel(cor.hex(bulk, pseudobulk(sr.pbmc, cpm = TRUE), "Bulk", "Short", "SuppFig3D"), "SuppFig3D")
save.panel(cor.hex(pseudobulk(sr.pbmc), pseudobulk(lr.pbmc), "Short", "BenchDrop-seq", "SuppFig3E", base.size = 18),
           "SuppFig3E")
