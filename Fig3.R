# Figure 3: panels 3A-3C, 3E, 3F. Inputs: SR/LR.PBMC.S3.rds, quant.sf, t2gnames.txt, genes.gtf
script.dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(if (length(script.dir)) script.dir else ".", "utils.R"))
library(patchwork)
library(ggrepel)
suppressPackageStartupMessages(library(rtracklayer))

sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")
t2g <- read.t2g("t2gnames.txt")
message(sprintf("PBMC cells: SR = %d, LR = %d", ncol(sr.pbmc), ncol(lr.pbmc)))

#### 3A-3C: fold change vs bulk, short-read vs BenchDrop-seq ####
sr.cpm <- pseudobulk(sr.pbmc, cpm = TRUE)
lr.tpm <- pseudobulk(lr.pbmc, assay = "transcript", cpm = TRUE)
lr.cpm <- to.gene(lr.tpm, t2g)
bulk.cpm <- bulk.gene.tpm("quant.sf", t2g)

genes <- Reduce(intersect, list(names(lr.cpm), names(sr.cpm), names(bulk.cpm)))
df <- data.frame(gname = genes, long = lr.cpm[genes], short = sr.cpm[genes], bulk = bulk.cpm[genes])
print(cor(df[, c("long", "short", "bulk")], method = "spearman"))
df$sfc <- log10((df$short + 1) / (df$bulk + 1))
df$lfc <- log10((df$long + 1) / (df$bulk + 1))

gtf <- import("genes.gtf")
gene.len <- data.frame(gname = gtf$gene_name[gtf$type == "gene"],
                       gene_length = width(gtf[gtf$type == "gene"])) %>%
  distinct(gname, .keep_all = TRUE)
tx <- gtf[gtf$type == "transcript"]
exon.sum <- data.frame(transcript_id = gtf$transcript_id[gtf$type == "exon"],
                       exon_length = width(gtf[gtf$type == "exon"])) %>%
  group_by(transcript_id) %>% summarise(total_exon_length = sum(exon_length))
intron.len <- data.frame(transcript_id = tx$transcript_id, gname = tx$gene_name,
                         transcript_length = width(tx)) %>%
  left_join(exon.sum, by = "transcript_id") %>%
  group_by(gname) %>%
  summarise(mean_intron_length = mean(transcript_length - total_exon_length))
df <- df %>% left_join(gene.len, by = "gname") %>% left_join(intron.len, by = "gname")

# N_total counts genes outside the radius in all three fold-change spaces (as in the paper inset).
radius <- 0.25
df$fc <- log((df$short + 1) / (df$long + 1))
removed <- with(df, sqrt(lfc^2 + sfc^2) <= radius & sqrt(lfc^2 + fc^2) <= radius & sqrt(sfc^2 + fc^2) <= radius)
n.total <- sum(!removed)
n.removed <- sum(removed)
plotted <- df[sqrt(df$lfc^2 + df$sfc^2) > radius, ]
message(sprintf("Fig3A: N_total = %d genes, N_removed = %d genes (%d genes drawn)", n.total, n.removed, nrow(plotted)))

labels <- c("HIST1H4F", "HIST2H2AB", "RPS4Y1", "DDX3Y", "PTPRCAP", "PTPRD", "RBFOX1", "CORO1B")
fc.scatter <- function(fill.layer, label.genes) {
  ggplot(plotted, aes(x = sfc, y = lfc)) +
    geom_point(alpha = 0.5) +
    fill.layer +
    geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.6, linetype = "dashed") +
    geom_label_repel(aes(label = ifelse(gname %in% label.genes, gname, "")),
                     size = 7, fontface = "bold", box.padding = 0.35, point.padding = 0.5,
                     segment.color = "grey50", max.overlaps = Inf) +
    labs(x = "Log2FC_1", y = "Log2FC_2") +
    theme_minimal()
}
with.marginals <- function(p) {
  top <- ggplot(plotted, aes(sfc)) + geom_density(fill = "#D5722A", color = "black") + theme_void()
  right <- ggplot(plotted, aes(lfc)) + geom_density(fill = "#5C1A8A", color = "black") +
    coord_flip() + theme_void()
  top + plot_spacer() + p + right + plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
}

p3a <- fc.scatter(geom_hex(bins = 100), labels) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Count") +
  annotate("label", x = 1.5, y = -3, hjust = 0, size = 5, fontface = "bold",
           label = sprintf("N_total = %d genes\nN_removed = %d genes", n.total, n.removed))
save.panel(with.marginals(p3a), "Fig3A", width = 10, height = 10)

hex.by <- function(z) stat_summary_hex(aes(z = log10(.data[[z]])), fun = mean, bins = 100)
save.panel(fc.scatter(hex.by("mean_intron_length"), c("PTPRD", "RBFOX1")) +
             scale_fill_viridis_c(option = "plasma", name = "Mean Intron Length (Log)"), "Fig3B")
save.panel(fc.scatter(hex.by("gene_length"), c("PTPRD", "RBFOX1")) +
             scale_fill_viridis_c(option = "plasma", name = "Gene Length (Log)"), "Fig3C")

#### 3E: gene-level UMAP with transferred cell types ####
lr.pbmc$celltype.grouped <- group.celltypes(sr.pbmc$predicted.celltype.l2)[match(colnames(lr.pbmc), colnames(sr.pbmc))]
p <- DimPlot(lr.pbmc, reduction = "ref.umap", group.by = "celltype.grouped", pt.size = 1.5) +
  NoLegend() + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(NULL)
save.panel(LabelClusters(p, id = "celltype.grouped", fontface = "bold", size = 5.5, repel = TRUE, color = "black"),
           "Fig3E", width = 9, height = 8)

#### 3F: marker-gene dot plot ####
lr.pbmc$celltype.grouped <- group.celltypes(lr.pbmc$predicted.celltype.l2, fine = TRUE)
lr.pbmc <- subset(lr.pbmc, !is.na(celltype.grouped))
Idents(lr.pbmc) <- "celltype.grouped"
markers <- FindAllMarkers(lr.pbmc, assay = "RNA", only.pos = TRUE) %>%
  mutate(spec.score = avg_log2FC * (pct.1 - pct.2)) %>%
  filter(avg_log2FC > 0.4, pct.1 > 0.25, pct.2 < 0.15) %>%
  arrange(cluster, desc(spec.score))
top3 <- c()
for (cl in unique(markers$cluster)) {
  top3 <- c(top3, head(setdiff(markers$gene[markers$cluster == cl], top3), 3))
}
message("Fig3F markers: ", paste(top3, collapse = ", "))
lr.pbmc <- ScaleData(lr.pbmc, features = top3, assay = "RNA")
p <- DotPlot(lr.pbmc, features = top3, assay = "RNA", dot.scale = 10) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        plot.margin = unit(c(1, 1, 1, 2), "cm"))
save.panel(p, "Fig3F", width = 16, height = 7)
