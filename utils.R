# Helpers shared by the figure scripts. Run scripts from the input directory; output goes to figures/.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

out.dir <- "figures"
dir.create(out.dir, showWarnings = FALSE)

save.panel <- function(p, name, width = 7, height = 6) {
  ggsave(file.path(out.dir, paste0(name, ".pdf")), p, width = width, height = height)
  invisible(p)
}

read.t2g <- function(path = "t2gnames.txt") {
  t2g <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  setNames(t2g$V2, t2g$V1)
}

to.gene <- function(x, t2g) {
  keep <- names(x)[names(x) %in% names(t2g)]
  tapply(x[keep], t2g[keep], sum)
}

bulk.gene.tpm <- function(path, t2g) {
  q <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  to.gene(setNames(q$TPM, q$Name), t2g)
}

pseudobulk <- function(obj, assay = "RNA", cpm = FALSE) {
  v <- rowSums(GetAssayData(obj, assay = assay, layer = "data"))
  if (cpm) v * 1e6 / sum(v) else v
}

align.union <- function(a, b) {
  g <- union(names(a), names(b))
  A <- setNames(numeric(length(g)), g)
  B <- A
  A[names(a)] <- a
  B[names(b)] <- b
  list(A, B)
}

cor.hex <- function(x, y, xlab, ylab, tag, base.size = 25) {
  v <- align.union(x, y)
  rho <- cor(v[[1]], v[[2]], method = "spearman")
  message(sprintf("%s: n = %d, Spearman = %.4f", tag, length(v[[1]]), rho))
  df <- data.frame(x = log1p(v[[1]]), y = log1p(v[[2]]))
  ggplot(df, aes(x, y)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(option = "plasma", trans = "log10") +
    labs(x = xlab, y = ylab, fill = "Density") +
    theme_minimal(base_size = base.size) +
    theme(axis.line = element_line(color = "black", linewidth = 0.8),
          panel.grid = element_blank()) +
    annotate("text", x = min(df$x), y = max(df$y), hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(rho, 3)))
}

read.length.density <- function(path, tag, short.read.length = 68) {
  df <- read.table(path, col.names = c("frequency", "length"))
  lens <- rep(df$length, df$frequency)
  wt <- as.numeric(df$length) * df$frequency
  o <- order(df$length, decreasing = TRUE)
  n50 <- df$length[o][which(cumsum(wt[o]) > sum(wt) / 2)[1]]
  message(sprintf("%s: reads = %d, mean = %.1f, median = %.0f, N50 = %d",
                  tag, length(lens), mean(lens), median(lens), n50))
  ggplot(data.frame(length = lens), aes(x = log2(length))) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 100) +
    geom_density(alpha = 0.2, fill = "#FF6666") +
    geom_vline(xintercept = log2(short.read.length), color = "blue", linetype = "dashed", linewidth = 1) +
    theme_classic()
}

group.celltypes <- function(l2, fine = FALSE) {
  l2 <- unname(l2)  # Seurat `$<-` matches named vectors by name
  out <- dplyr::case_when(
    l2 %in% c("cDC1", "cDC2", if (fine) "ASDC")                  ~ "DC",
    l2 %in% c("NK", "NK Proliferating", "NK_CD56bright")         ~ "NK",
    l2 %in% c("B memory", "B naive", "B intermediate")           ~ "B Cells",
    l2 %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL")      ~ "CD4 T Cells",
    l2 %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM")                 ~ "CD8 T Cells",
    l2 == "Eryth"                                                ~ "Erythroid",
    l2 %in% c("dnT", "ILC")                                      ~ NA_character_,
    TRUE                                                         ~ l2
  )
  if (fine) {
    keep <- c("CD14 Mono", "CD16 Mono", "DC", "Plasmablast", "B Cells", "HSPC", "Erythroid",
              "Platelet", "pDC", "CD4 T Cells", "Treg", "CD8 T Cells", "MAIT",
              "CD4 Proliferating", "NK")
    out[!out %in% keep] <- NA_character_
  }
  out
}
