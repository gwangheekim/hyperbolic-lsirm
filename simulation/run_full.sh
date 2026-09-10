#!/usr/bin/env bash
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
export LSIRM_PRES="$HERE"
export LSIRM_WORKSPACE="$(cd "$HERE/.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export SIM_NSEEDS="${SIM_NSEEDS:-100}"
export SIM_SEEDS="$SIM_NSEEDS"
export LSIRM_MC_CORES="${LSIRM_MC_CORES:-4}"
cd "$HERE"
mkdir -p reproduce/fits reproduce/results logs "$LSIRM_WORKSPACE/figure" "$LSIRM_WORKSPACE/table" "$LSIRM_WORKSPACE/result/simulation/regenerated"

Rscript reproduce/preflight.R
Rscript reproduce/fit_sim.R
Rscript reproduce/table_distortion.R
ROB_NSEED="$SIM_NSEEDS" Rscript reproduce/sim_robustness.R
NS_SEEDS="$(seq -s, 1 "$SIM_NSEEDS")" Rscript reproduce/sim_n_sweep.R
NEG_SEEDS="$SIM_NSEEDS" Rscript reproduce/sim_negative_control.R
Rscript reproduce/table_robustness.R
Rscript reproduce/validate_results.R
cp reproduce/results/results_*.csv "$LSIRM_WORKSPACE/result/simulation/regenerated/"
echo "simulation reproduction complete"
