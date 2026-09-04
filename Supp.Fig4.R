# Supp. Figure 4: panel S4A. Inputs: LR.PBMC.S3.rds, quant.sf, masiso/iso/
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))
library(patchwork)

lr.pbmc <- readRDS("LR.PBMC.S3.rds")
lr.tpm <- pseudobulk(lr.pbmc, assay = "transcript", cpm = TRUE)

bulk <- read.table("quant.sf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bulk.tpm <- setNames(bulk$TPM, bulk$Name)

masiso <- NormalizeData(CreateSeuratObject(Read10X("masiso/iso", gene.column = 1)), verbose = FALSE)
masiso.tpm <- pseudobulk(masiso, cpm = TRUE)
names(masiso.tpm) <- sub("\\..*$", "", names(masiso.tpm))

txps <- Reduce(intersect, list(names(lr.tpm), names(bulk.tpm), names(masiso.tpm)))
df <- data.frame(masiso = masiso.tpm[txps], long = lr.tpm[txps], bulk = bulk.tpm[txps])
df$masiso <- df$masiso * 1e6 / sum(df$masiso)
print(cor(df, method = "spearman"))
df$mfc <- log10((df$masiso + 1) / (df$bulk + 1))
df$lfc <- log10((df$long + 1) / (df$bulk + 1))

radius <- 0.25
plotted <- df[sqrt(df$lfc^2 + df$mfc^2) > radius, ]
n.total <- nrow(plotted)
n.removed <- nrow(df) - nrow(plotted)
message(sprintf("SuppFig4A: N_total = %d transcripts, N_removed = %d transcripts", n.total, n.removed))

p <- ggplot(plotted, aes(x = mfc, y = lfc)) +
  geom_point() +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "count") +
  labs(x = "Log2FC_1", y = "Log2FC_2") +
  annotate("label", x = 0.5, y = -4, hjust = 0, size = 5, fontface = "bold",
           label = sprintf("N_total = %d transcripts\nN_removed = %d transcripts", n.total, n.removed)) +
  theme_minimal()
top <- ggplot(plotted, aes(mfc)) + geom_density(fill = "#D5722A", color = "black") + theme_void()
right <- ggplot(plotted, aes(lfc)) + geom_density(fill = "#5C1A8A", color = "black") + coord_flip() + theme_void()
save.panel(top + plot_spacer() + p + right + plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4)),
           "SuppFig4A", width = 8, height = 10)
