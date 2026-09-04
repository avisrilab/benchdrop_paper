# Figure 2: panels 2B, 2C, 2E-2H. Inputs: SR/LR.K562.S3.rds, lr.lens.txt, sr.lens.txt,
# k562.pe.bulk.quant.sf, t2gnames.txt, Downsampling/{10..90}/
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))
library(Matrix.utils)
library(patchwork)
library(viridis)

sr.k562 <- readRDS("SR.K562.S3.rds")
lr.k562 <- readRDS("LR.K562.S3.rds")
t2g <- read.t2g("t2gnames.txt")
message(sprintf("K562 cells: SR = %d, LR = %d", ncol(sr.k562), ncol(lr.k562)))

#### 2B: UMIs and genes per cell ####
meta <- bind_rows(
  data.frame(sr.k562@meta.data[, c("nCount_RNA", "nFeature_RNA")], sample = "Short"),
  data.frame(lr.k562@meta.data[, c("nCount_RNA", "nFeature_RNA")], sample = "BenchDrop-seq")
)
violin <- function(y, ylab) {
  ggplot(meta, aes(x = sample, y = .data[[y]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_y_log10() +
    theme_minimal(base_size = 25) +
    labs(x = "Sample", y = ylab) +
    theme(legend.position = "none", panel.grid = element_blank())
}
save.panel(violin("nCount_RNA", "UMIs per cell") + violin("nFeature_RNA", "Genes per cell"),
           "Fig2B", width = 10, height = 6)

#### 2C: equivalence-class lengths ####
eqc <- function(path, kind) {
  d <- as.data.frame(table(read.table(path)$V1))
  data.frame(length = d$Var1, percent = d$Freq * 100 / sum(d$Freq), kind = kind)
}
df <- rbind(eqc("lr.lens.txt", "BenchDrop-seq"), eqc("sr.lens.txt", "Short"))
p <- ggplot(df, aes(x = kind, y = length, color = kind, size = percent)) +
  geom_point() +
  theme_classic() +
  labs(title = "Equivalence Class Length", x = "Platform", y = "Length") +
  scale_size_continuous(range = c(2, 12)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))
save.panel(p, "Fig2C", width = 6, height = 8)

#### 2E-2G: pseudobulk correlations ####
bulk <- bulk.gene.tpm("k562.pe.bulk.quant.sf", t2g)
sr.pb <- pseudobulk(sr.k562, cpm = TRUE)
lr.pb <- pseudobulk(lr.k562, cpm = TRUE)
save.panel(cor.hex(bulk, lr.pb, "Bulk", "BenchDrop-seq", "Fig2E"), "Fig2E")
save.panel(cor.hex(bulk, sr.pb, "Bulk", "Short", "Fig2F"), "Fig2F")
save.panel(cor.hex(pseudobulk(sr.k562), pseudobulk(lr.k562), "Short", "BenchDrop-seq", "Fig2G"), "Fig2G")

#### 2H: per-cell correlation of full-depth vs subsampled long-read data ####
cells <- colnames(sr.k562)

gene.counts <- function(tx.mat) {
  tx.mat <- tx.mat[, cells]
  aggregate.Matrix(tx.mat, groupings = as.character(t2g[rownames(tx.mat)]), fun = "sum", MARGIN = 2)
}
read.mtx <- function(dir) {
  m <- t(readMM(file.path(dir, "matrix.mtx")))
  rownames(m) <- read.table(file.path(dir, "features.tsv"))$V1
  colnames(m) <- read.table(file.path(dir, "barcodes.tsv"))$V1
  m
}
levels <- c("1.0", paste0("0.", 9:1))
mats <- c(list(gene.counts(GetAssayData(lr.k562, assay = "transcript", layer = "counts"))),
          lapply(seq(90, 10, -10), function(f) gene.counts(read.mtx(file.path("Downsampling", f)))))
genes <- Reduce(union, lapply(mats, rownames))
mats <- lapply(mats, function(m) {
  full <- Matrix(0, length(genes), length(cells), sparse = TRUE, dimnames = list(genes, cells))
  full[rownames(m), ] <- m
  full
})
corrs <- t(vapply(cells, function(cn) {
  cor(do.call(cbind, lapply(mats, function(m) as.numeric(m[, cn]))), method = "pearson")[1, ]
}, numeric(length(mats))))
colnames(corrs) <- levels
df <- reshape2::melt(corrs, varnames = c("cell", "Subsampling"), value.name = "Correlation")
df$Subsampling <- factor(df$Subsampling, levels = levels)
print(tapply(df$Correlation, df$Subsampling, median))

p <- ggplot(df, aes(x = Subsampling, y = Correlation, fill = Subsampling)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal(base_size = 25) +
  labs(x = "Subsampling Level", y = "Correlation") +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.title = element_text(size = 20), axis.text = element_text(size = 20))
save.panel(p, "Fig2H", width = 14, height = 6)
