#!/usr/bin/env bash
# Staircase rungs 1-2: run the full chain at one or more ONT error scales.
#   0.0 = rung 1 (clean ceiling), 1.0 = rung 2 (empirical K562 profile), >1 = stress.
# Usage: run_error_sweep.sh 0.0 1.0 2.0
CODE=${BENCHDROP_REPO}/staircase
. "$CODE/_env.sh"                  # set -euo pipefail, PATH, $PY, require/fresh/made

STA=${BENCHDROP_RESULTS}/staircase
B=${BAGPIPER_BIN}
WL=${BAGPIPER_DIR}/test/whitelist/pip.v5.csv
require minimap2 samtools "$PY" "$B"
cd "$STA"

for S in "$@"; do
  echo "########## ERROR_SCALE=${S} ##########"
  # The barcode stage's output is in this list because its `|| true` (grep exits 1 when it
  # matches nothing) would otherwise swallow a bagpiper failure and leave the PREVIOUS run's
  # barcode output in place for the aligner to consume. fresh + made close that.
  fresh sim.rung1.fastq.gz out/barcode/passed.bcd.nanopore.fa.gz aligned.bam eqclass.tsv
  $PY "$CODE/simulate_rung1.py" "$S" | tail -1
  made sim.rung1.fastq.gz

  $B barcode --r1 sim.rung1.fastq.gz --whitelist "$WL" --nanopore -o out/barcode -t 8 2>&1 \
    | grep -E 'Total:|Matched:|Mismatch:' || true
  made out/barcode/passed.bcd.nanopore.fa.gz

  minimap2 -ax map-ont --for-only -N 200 -p 0.9 -t 8 txome.mmi \
    out/barcode/passed.bcd.nanopore.fa.gz | samtools sort -n -@4 -o aligned.bam
  made aligned.bam
  echo "mapped: $(samtools view -c -F 4 aligned.bam)"

  $B eqclass --b1 aligned.bam -o eqclass.tsv >/dev/null
  made eqclass.tsv
  $B count --eq eqclass.tsv -o counts -t 8 >/dev/null

  $PY "$CODE/evaluate_rung1.py" | tail -4
  echo "--- identifiability ---"
  $PY "$CODE/identifiability.py" | tail -8
  echo
done
