#!/usr/bin/env bash
# Rung 3: sweep per-cell depth at the empirical ONT error rate -> the identifiability floor.
# Usage: run_depth_sweep.sh 1 2 4 8 16 32      (ERR=1.0 by default)
CODE=${BENCHDROP_REPO}/staircase
. "$CODE/_env.sh"                  # set -euo pipefail, PATH, $PY, require/fresh/made

STA=${BENCHDROP_RESULTS}/staircase
B=${BAGPIPER_BIN}
WL=${BAGPIPER_DIR}/test/whitelist/pip.v5.csv
require minimap2 samtools "$PY" "$B"
cd "$STA"
ERR=${ERR:-1.0}
rm -f sweep_results.tsv

for D in "$@"; do
  echo "########## DEPTH=${D} reads/gene/cell (error=${ERR}x empirical) ##########"
  # The barcode stage's output is in this list because its `|| true` (grep exits 1 when it
  # matches nothing) would otherwise swallow a bagpiper failure and leave the PREVIOUS run's
  # barcode output in place for the aligner to consume. fresh + made close that.
  fresh sim.rung1.fastq.gz out/barcode/passed.bcd.nanopore.fa.gz aligned.bam eqclass.tsv
  $PY "$CODE/simulate_rung1.py" "$ERR" "$D" | tail -1
  made sim.rung1.fastq.gz

  $B barcode --r1 sim.rung1.fastq.gz --whitelist "$WL" --nanopore -o out/barcode -t 8 2>&1 \
    | grep -E 'Matched:' || true
  made out/barcode/passed.bcd.nanopore.fa.gz

  minimap2 -ax map-ont --for-only -N 200 -p 0.9 -t 8 txome.mmi \
    out/barcode/passed.bcd.nanopore.fa.gz | samtools sort -n -@4 -o aligned.bam
  made aligned.bam

  $B eqclass --b1 aligned.bam -o eqclass.tsv >/dev/null
  made eqclass.tsv
  $B count --eq eqclass.tsv -o counts -t 8 >/dev/null

  $PY "$CODE/evaluate_rung1.py" | tail -2
  # Keep the per-depth tables so the metrics can be re-derived without re-running the pipeline.
  cp -f percelltx.rung1.tsv "percelltx.depth${D}.tsv"
  cp -f percellgene.rung1.tsv "percellgene.depth${D}.tsv"
  $PY "$CODE/identifiability.py" "depth${D}" >/dev/null
  echo
done

echo "=================================================================="
$PY "$CODE/floor_table.py"
