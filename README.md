# BenchDrop-seq analysis

Code behind the figures and tables of *BenchDrop-seq: a microfluidics-free workflow for benchtop
single-cell long-read RNA sequencing* (bioRxiv, doi:10.64898/2026.03.12.706999). Figure and table
numbers follow the revised manuscript.

Raw and processed BenchDrop-seq short-read and long-read data: GEO accession GSE318099.
Bulk PBMC RNA-seq used for benchmarking: ENCODE SRX2370564. Bulk K562 RNA-seq: ENCODE (accession
given in the manuscript's Availability section). ENCODE sorted immune-population bulk RNA-seq used
for the isoform-usage silver standard (B cells, ENCFF795CTQ/ENCFF959BDM; T cells, ENCFF970FKA;
NK cells, ENCFF738YZJ; monocytes, ENCFF774DLJ/ENCFF831SVM).
Long reads were processed with Bagpiper v0.1.0 (https://github.com/avisrilab/bagpiper).

## Scripts

| Script | Panels |
|---|---|
| `Fig2.R` | 2B, 2C, 2E, 2F, 2G, 2H |
| `Fig3.R` | 3A, 3B, 3C, 3E, 3F |
| `Fig4.py`, `staircase/` | 4A-4D, S6 (simulation and the coverage criterion, Methods 2.5) |
| `Fig5.R` | 5A-5F |
| `Fig6.py` | 6A-6D, Supplementary Table 3 (Methods 2.8) |
| `Table1/` | Table 1, Supplementary Table 2 (Methods 2.9) |
| `Table2.py` | Table 2 (Methods 1.5) |
| `Supp.Fig2.R` | S2A, S2D, S2E |
| `Supp.Fig3.R` | S3A-S3E |
| `Supp.Fig4.R` | S4A |
| `Supp.Fig5.R` | S5A-S5G |
| `Supp.Fig7.py` | S7 (cell hashing, Methods 1.6) |
| `Supp.Fig8.py` | S8 (Scrublet doublets, Methods 2.1) |
| `misc/coverage_tracks.py` | the genome-browser track renderer behind 2D, 3D, S2B and S2C |
| `utils.R` | helpers shared by the R scripts |

Figures 1, 2A, S1A and S1B are illustrations. S1C (pipeline runtime) is a bar chart of the wall-clock
times reported in Methods 2.4. 2D, 3D, S1C, S2B, S2C and S2E ship as the published renders. The
Fig. 3F marker genes are listed in Supplementary Table 1.

Other analyses reported in the text live in `misc/`: `error_control.py` (same-length error control,
Methods 2.5), `floor_dtu.py` (which PBMC genes clear the coverage criterion), `auprc_bulk.py` and
`auprc_sweep.py` (matched-bulk isoform-switch detection and its threshold sweep, Methods 2.6),
`bcr_feasibility.py` (per-cell immunoglobulin chain calls, Methods 2.7), `barcode_collision.py`
(barcode collision bound, Methods 2.1), `phasing_coverage_profile.py` (read geometry, the 382 bp
window, Methods 2.5), `export_reference.py` and `transfer_labels.R` (Azimuth PBMC label transfer
onto the short-read matrix, Methods 2.8), and `azimuth_labels.py`, which `Fig6.py` calls.

## Inputs

R scripts: place the files below in `input/` (git-ignored) and run every script from there. The
Seurat objects are those built for the paper (Seurat v5, Azimuth PBMC reference 2.10); the remaining
files are Bagpiper, PIPseeker, salmon or kallisto-bustools outputs and the GENCODE v32 annotation.

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
| `*.bam`, `genome.fa`, `genes.bed.gz` | genome alignments and reference for `misc/coverage_tracks.py` |

Python analyses: paths come from environment variables. `BENCHDROP_FEED` is the data root (the
PBMC and K562 long-read matrices and reads from GSE318099, the GENCODE v32 genome and annotation,
the hashing libraries, the isonome reads SRR37807425, the ENCODE quantifications, the Azimuth PBMC
reference, and `pbmc/pbmc_celltypes.tsv`, the barcode-to-lineage table the matched-bulk scripts
read); `BENCHDROP_RESULTS` and `BENCHDROP_SILO` are where outputs go; `BENCHDROP_REPO` is this
directory; `BAGPIPER_BIN` and `BAGPIPER_DIR` point at Bagpiper; `BENCHDROP_READS_BIN` is a bin
directory with minimap2, samtools and a Python with pysam; `BENCHDROP_PYTHON` is a Python with
the scanpy stack; `NANOSIM_DIR` is a NanoSim checkout. Each Python analysis is one file: parameters
and paths at the top, the steps in order, a runner at the bottom (`--from-stage NAME` resumes).

## Running

```
cd input
Rscript ../Fig2.R                 # writes ../figures/Fig2B.pdf, Fig2C.pdf, ...
python ../misc/coverage_tracks.py
$BENCHDROP_PYTHON ../Fig6.py
```

Each R script prints the numbers quoted in the paper (cell counts, Spearman correlations, N_total and
N_removed, read-length statistics) to the console.

## Environment

`renv.lock` records R 4.6.1, Bioconductor 3.23 and the exact version of every R package the scripts load
(Seurat 5.5.1, Signac 1.17.1, ggplot2 4.0.3, Matrix 1.7-5, ...) together with their dependencies.
Recreate the library with:

```
Rscript -e 'install.packages("renv"); renv::restore(lockfile = "renv.lock")'
```

`Matrix.utils` is no longer on CRAN; renv installs it from the CRAN archive.
Python: integrative_transcriptomics_viewer (for `misc/coverage_tracks.py` only). The Python analyses ran
under Python 3.11 with scanpy 1.10.4, anndata 0.10.9, pandas 2.2.3, numpy 1.26.4, matplotlib 3.11.0,
pysam 0.22.1, minimap2 2.31, samtools 1.21 and NanoSim 3.2.3; the benchmark used IsoQuant 3.10.0 and
bambu 3.14.0 (its renv lock is `Table1/bambu_renv/renv.lock`).
