# BenchDrop-seq analysis

Figure code for *BenchDrop-seq: a microfluidics-free platform for benchtop single-cell long-read RNA
sequencing* (bioRxiv, doi:10.64898/2026.03.12.706999).

Raw and processed BenchDrop-seq short-read and long-read data: GEO accession GSE318099.
Bulk K562 and PBMC RNA-seq used for benchmarking: ENCODE (SRX2370564).
Long reads were processed with Bagpiper v0.1.0 (https://github.com/avisrilab/bagpiper).

## Scripts

| Script | Panels |
|---|---|
| `Fig2.R` | 2B, 2C, 2E, 2F, 2G, 2H |
| `Fig3.R` | 3A, 3B, 3C, 3E, 3F |
| `Fig4.R` | 4A-4F |
| `Supp.Fig2.R` | S2A, S2D, S2E |
| `Supp.Fig3.R` | S3A-S3E |
| `Supp.Fig4.R` | S4A |
| `Supp.Fig5.R` | S5A-S5G |
| `coverage_tracks.py` | 2D, 3D, S2B, S2C (genome-browser tracks) |
| `utils.R` | helpers shared by the R scripts |

Figures 1, 2A, S1A and S1B are illustrations. S1C (pipeline runtime) is a bar chart of the wall-clock
times reported in Methods 2.4.

## Inputs

Place the files below in `input/` (git-ignored) and run every script from there. The Seurat objects are
those built for the paper (Seurat v5, Azimuth PBMC reference 2.10); the remaining files are Bagpiper,
PIPseeker, salmon or kallisto-bustools outputs and the GENCODE v32 annotation.

| File | Content |
|---|---|
| `SR.K562.S3.rds`, `LR.K562.S3.rds` | K562 short-read and BenchDrop-seq Seurat objects (2,468 matched cells) |
| `SR.PBMC.S3.rds`, `LR.PBMC.S3.rds` | PBMC short-read and BenchDrop-seq Seurat objects, Azimuth labels, `transcript` assay |
| `t2gnames.txt` | transcript ID to gene name, two tab-separated columns |
| `genes.gtf` | GENCODE v32 annotation |
| `k562.pe.bulk.quant.sf`, `quant.sf` | salmon quantification of bulk K562 and bulk PBMC RNA-seq |
| `lr.lens.txt`, `sr.lens.txt` | equivalence-class size per long read and per short read |
| `read.lens.k562.txt`, `read.lens.pbmc.txt` | `uniq -c` tables of long-read lengths |
| `Downsampling/{10..90}/` | Bagpiper transcript matrices from 10-90 % subsampled K562 long reads |
| `alevin/`, `kb/` | gene-level Salmon-Alevin and kallisto-bustools quantification of the K562 short reads |
| `masiso/iso/` | MAS-ISO-seq PBMC transcript counts (10x matrix layout) |
| `*.bam`, `genome.fa`, `genes.bed.gz` | genome alignments and reference for `coverage_tracks.py` |

## Running

```
cd input
Rscript ../Fig2.R                 # writes ../figures/Fig2B.pdf, Fig2C.pdf, ...
python ../coverage_tracks.py
```

Each script prints the numbers quoted in the paper (cell counts, Spearman correlations, N_total and
N_removed, read-length statistics) to the console.

## Environment

`renv.lock` records R 4.6.1, Bioconductor 3.23 and the exact version of every R package the scripts load
(Seurat 5.5.1, Signac 1.17.1, ggplot2 4.0.3, Matrix 1.7-5, ...) together with their dependencies.
Recreate the library with:

```
Rscript -e 'install.packages("renv"); renv::restore(lockfile = "renv.lock")'
```

`Matrix.utils` is no longer on CRAN; renv installs it from the CRAN archive.
Python: integrative_transcriptomics_viewer (for `coverage_tracks.py` only).
