# Sourced by the staircase runners: environment (BENCHDROP_READS_BIN, BAGPIPER_BIN, NanoSim) and
# the helper functions require/fresh/made.

set -euo pipefail

# The reads toolchain (minimap2, samtools, python with pysam, NanoSim deps) is supplied by the
# environment: BENCHDROP_READS_BIN is a bin directory holding them (see README.md).
READS_BIN="${BENCHDROP_READS_BIN:?set BENCHDROP_READS_BIN to the bin dir with minimap2, samtools and python}"
READS="$(dirname "$READS_BIN")"
export PATH="$READS_BIN:$PATH"
PY="$READS_BIN/python"

require() {   # preflight the whole toolchain: fail in one second, not in three hours
  local t
  for t in "$@"; do
    command -v "$t" >/dev/null || { echo "FATAL: '$t' is not on PATH ($READS/bin)" >&2; exit 1; }
  done
}

fresh() { rm -f "$@"; }   # call BEFORE the stage that writes these

made() {      # call AFTER: a stage that exits 0 with no output is still a failed stage
  local f
  for f in "$@"; do
    [ -s "$f" ] || { echo "FATAL: stage exited 0 but produced no '$f'" >&2; exit 1; }
  done
}
