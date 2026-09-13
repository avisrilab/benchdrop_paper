#!/usr/bin/env Rscript
# bambu single-cell quantification for the benchmark: annotation-guided against GENCODE v32,
# one sample per cell. Usage: Rscript bambu_sc.R <cells_bam_dir> <genome.fa> <genes.gtf> <outdir>
suppressMessages(library(bambu))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) stop("usage: bambu_sc.R <cells_bam_dir> <genome.fa> <genes.gtf> <outdir>")
cells_dir <- args[[1]]; genome <- args[[2]]; gtf <- args[[3]]; outdir <- args[[4]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bams <- sort(list.files(cells_dir, pattern = "\\.bam$", full.names = TRUE))
if (!length(bams)) stop("no per-cell BAMs in ", cells_dir)
message(sprintf("bambu: %d per-cell BAMs (samples)", length(bams)))

annotations <- prepareAnnotations(gtf)
# ncore = 1 (SerialParam), NOT parallel: bambu's MulticoreParam fork workers spuriously hit
# "reached CPU time limit" on macOS (ulimit -t is unlimited, so it is BiocParallel forking, not an
# OS cap). Serial is robust and bounded here, we cap at TOPN cells.
se <- bambu(reads = bams, annotations = annotations, genome = genome,
            discovery = FALSE, quant = TRUE, ncore = 1L, lowMemory = TRUE, verbose = TRUE)

# writeBambuOutput emits counts.txt (transcript x sample) + CPM + the annotation; sample columns are
# the per-cell BAM basenames (= the CB), so this is transcript x cell.
writeBambuOutput(se, path = outdir)
saveRDS(se, file = file.path(outdir, "bambu_se.rds"))
message(sprintf("bambu: wrote transcript x cell counts -> %s", outdir))
