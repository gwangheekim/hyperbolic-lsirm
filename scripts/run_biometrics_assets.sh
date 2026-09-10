#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
Rscript "$ROOT/simulation/render_reported.R"
Rscript "$ROOT/simulation/render_biometrics_truths.R"
Rscript "$ROOT/R/figures/table_dementia_case.R"
Rscript "$ROOT/R/figures/figure_dementia_radial_tree.R"
Rscript "$ROOT/scripts/validate_biometrics_submission.R"
