#!/usr/bin/env bash
## check_results.sh — completeness report for reproduce/fits/ and reproduce/results_*.csv (or, under
## submission/simulation and reproduce/results/*.csv). Run standalone any time, or automatically at the end
## of run_fits_parallel.sh. Lists every file the canonical generators and summary builders are
## supposed to produce and marks it present/MISSING, so a server run's completeness can
## be checked at a glance instead of re-deriving "what should exist" from each script by hand.
##
##   cd reproduce && ./check_results.sh
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"
RESDIR="reproduce/results"; [ -d "$RESDIR" ] || RESDIR="reproduce"   # main project has no results/ subfolder

ok=0; missing=0
check() {
  if [ -f "$1" ]; then ok=$((ok+1)); else echo "  MISSING: $1"; missing=$((missing+1)); fi
}

NSEEDS="${SIM_NSEEDS:-100}"
echo "== fits/ (expected seeds 1..$NSEEDS) =="
check "reproduce/fits/sim_truth.rds"
check "reproduce/fits/sim_config.rds"
for s in $(seq 1 "$NSEEDS"); do for tag in H2 H3 R2 R3 R4 R5; do check "reproduce/fits/sim_s${s}_${tag}.rds"; done; done
echo "== paper-facing result tables =="
check "$RESDIR/results_negative_control.csv"
check "$RESDIR/results_robustness.csv"
check "$RESDIR/results_n_sweep.csv"
check "$RESDIR/results_distortion_raw.csv"
check "$RESDIR/results_distortion.csv"
check "../table/table_distortion.tex"
check "../table/table_robustness.tex"

echo
echo "$ok present, $missing missing (out of $((ok+missing)))"
if [ "$missing" -eq 0 ]; then
  SIM_NSEEDS="$NSEEDS" Rscript reproduce/validate_results.R
  echo "All expected outputs present and internally consistent."
  exit 0
fi
echo "Re-run the missing generator(s) individually, e.g.: Rscript reproduce/<script>.R" >&2
exit 1
