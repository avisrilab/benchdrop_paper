#!/usr/bin/env bash
# ============================================================================
# Public-data benchmark: Bagpiper, bambu and IsoQuant on the isonome PIP-seq+ONT PBMC reads
# (SRR37807425), each quantifier given the same per-cell barcode assignment (Table 1,
# Supplementary Table 2, Methods 2.9). FLAMES and wf-single-cell are cited, not run: neither
# processes the PIP-seq barcode. Inputs under BENCHDROP_FEED; outputs under BENCHDROP_RESULTS/table1.
# Timing uses macOS /usr/bin/time -l (wall clock and peak RSS), so the efficiency rows reproduce
# on macOS as run.
# ============================================================================
set -euo pipefail

# ---- inputs: BENCHDROP_FEED = raw inputs; GENCODE v32 ----------------------
FEED=${BENCHDROP_FEED}
GENOME=$FEED/genome/genome.fa                 # GRCh38, GENCODE v32 (genome.md)
GTF=$FEED/genome/genes.gtf                     # GENCODE v32
TXOME=$FEED/genome/transcriptome.fa            # GENCODE v32 transcripts (bagpiper aligns here)
WL=$FEED/barcodes/pip_v4_hybrid_barcodes.csv  # PIP-seq v4 barcode whitelist
READS="${1:-$FEED/public/isonome/SRR37807425.fastq.gz}"   # isonome subset (prefetch + fastq-dump -X)

# ---- tools ----------------------------------------------------------------
# Tools come from the environment (see README.md): BENCHDROP_READS_BIN holds minimap2,
# samtools and a python with pysam; BAGPIPER_BIN is the bagpiper 0.1.0 binary.
RBIN="${BENCHDROP_READS_BIN:?set BENCHDROP_READS_BIN}"
BAG="${BAGPIPER_BIN:?set BAGPIPER_BIN to the bagpiper 0.1.0 binary}"
MM2=$RBIN/minimap2
ST=$RBIN/samtools
PY=$RBIN/python
HERE=$(cd "$(dirname "$0")" && pwd)
# bambu is Bioconductor R; restore its library once with `Rscript -e 'renv::restore()'` inside
# bambu_renv/ (README). Rscript activates that library when run from there.
BAMBU_RENV="${BAMBU_RENV:-$HERE/bambu_renv}"

THREADS="${THREADS:-8}"
TOPN="${TOPN:-300}"                            # cells for the bambu per-sample arm (bounded compute)

# tunable alignment recipes -------------------------------------------------
TXOME_MM2=(-ax map-ont --for-only -N 200 -p 0.9)     # bagpiper's pinned transcriptome recipe
GENOME_MM2=(-ax splice -uf -k14)                     # ONT cDNA spliced genome align (IsoQuant/bambu)

OUT=${BENCHDROP_RESULTS}/table1
mkdir -p "$OUT"; cd "$OUT"
ulimit -n 4096                                 # the bambu split opens one writer per kept cell

# ---- preflight: fail in one second, not mid-run ---------------------------
for f in "$GENOME" "$GTF" "$TXOME" "$WL" "$READS" "$BAG"; do
  [ -s "$f" ] || { echo "FATAL: missing input $f" >&2; exit 1; }
done
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH (bambu arm needs system R + renv)" >&2; exit 1; }
TIME=/usr/bin/time     # macOS BSD time; `-l` prints '<sec> real' + '<bytes> maximum resident set size' (metrics.py parses both)
[ -x "$TIME" ] || { echo "FATAL: $TIME not found (timing wrappers need macOS /usr/bin/time -l)" >&2; exit 1; }

echo "== 0. shared barcode assignment (bagpiper barcode; CB+UMI -> read name) =="
if [ -s bc/passed.bcd.nanopore.fa.gz ]; then
  echo "   (reusing existing bc/ barcode assignment; delete bc/ to force a re-run)"
else
  "$BAG" barcode --r1 "$READS" --whitelist "$WL" --nanopore -o bc
fi
# Select the PASSED barcoded fasta ONLY. Trap: 'passed.bcd.nanopore.fa.gz' also matches the
# '*.nanopore.fa.gz' fallback pattern, and 'failed' sorts BEFORE 'passed', so a combined
# 'ls bc/passed* bc/*.nanopore* | head -1' silently returns FAILED (the barcodeless reads) and voids
# every arm. Prefer passed*; only then fall back to *.nanopore.fa.gz with failed excluded.
BCFA=$(ls bc/passed*.fa.gz 2>/dev/null | head -1)                       # >origid_CB_UMI + cDNA
[ -s "$BCFA" ] || BCFA=$(ls bc/*.nanopore.fa.gz 2>/dev/null | grep -v '/failed' | head -1)
[ -s "$BCFA" ] || { echo "FATAL: no barcoded fasta under bc/" >&2; exit 1; }

echo "== 1. bagpiper arm (native: transcriptome align + count EM) =="
t0=$(date +%s)
"$MM2" "${TXOME_MM2[@]}" -t "$THREADS" "$TXOME" "$BCFA" | "$ST" view -b -o aln.txome.bam -
printf '%s real\n' "$(( $(date +%s) - t0 ))" > time.txome_align    # bagpiper's own align (wall only; pipeline)
"$TIME" -l "$BAG" count --b1 aln.txome.bam -o bagpiper 2> >(tee time.bagpiper >&2)   # quantifier: wall + peak RSS

echo "== 2. shared genome alignment + CB tag (feeds IsoQuant + bambu) =="
t0=$(date +%s)
"$MM2" "${GENOME_MM2[@]}" -t "$THREADS" "$GENOME" "$BCFA" | "$ST" sort -@ "$THREADS" -o aln.genome.bam -
"$ST" index aln.genome.bam
"$PY" "$HERE/inject_cb.py" aln.genome.bam aln.genome.cb.unsorted.bam
"$ST" sort -@ "$THREADS" -o aln.genome.cb.bam aln.genome.cb.unsorted.bam
"$ST" index aln.genome.cb.bam
rm -f aln.genome.cb.unsorted.bam
printf '%s real\n' "$(( $(date +%s) - t0 ))" > time.genome_align   # shared align (IsoQuant + bambu; wall only)

echo "== 3. IsoQuant arms (native: genome, per-CB group; ephemeral pixi env) =="
pixi exec --spec isoquant -c bioconda -c conda-forge isoquant.py --version >> versions.txt
# /usr/bin/time wraps the whole pixi exec and rolls up the isoquant.py descendant's peak RSS.
#
# Two IsoQuant arms. IsoQuant's default --transcript_quantification is `unique_only`: a read is
# counted only where it is unambiguously assignable to one isoform, which on these reads discards
# 43.6% of assigned reads (ambiguous plus inconsistent_ambiguous). Bagpiper and bambu both allocate
# those reads by EM, so `with_ambiguous` (ambiguous reads on, inconsistent alignments still off;
# IsoQuant v3.10.0 src/long_read_counter.py) is the like-for-like mode. Both arms are reported:
# the default is what a reader gets out of the box, the ambiguous-aware arm is the comparison.
# --transcript_quantification is hidden from --help; see --full_help.
ISOQUANT_MODE="${ISOQUANT_MODE:-with_ambiguous}"
for arm in "unique_only:isoquant" "$ISOQUANT_MODE:isoquant_wa"; do
  mode=${arm%%:*}; dir=${arm##*:}
  if [ -s "$dir/OUT/OUT.transcript_grouped_counts.matrix.mtx" ]; then
    echo "   (reusing existing $dir/; delete it to force a re-run)"
    continue
  fi
  "$TIME" -l pixi exec --spec isoquant -c bioconda -c conda-forge isoquant.py \
      --reference "$GENOME" --genedb "$GTF" --bam aln.genome.cb.bam \
      --data_type nanopore --read_group tag:CB \
      --transcript_quantification "$mode" --output "$dir" 2> >(tee "time.$dir" >&2)
done

echo "== 4. bambu arm (the isonome paper's tool; per-cell, discovery=FALSE) =="
"$PY" "$HERE/split_cells.py" aln.genome.cb.bam cells "$TOPN"
( cd "$BAMBU_RENV" && "$TIME" -l Rscript "$HERE/bambu_sc.R" "$OUT/cells" "$GENOME" "$GTF" "$OUT/bambu" ) 2> >(tee "$OUT/time.bambu" >&2)

echo "== 5. metrics on identical inputs =="
"$PY" "$HERE/metrics.py" "$OUT"
echo "done -> $OUT/metrics.tsv"
