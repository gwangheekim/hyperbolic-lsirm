#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/.."
# Rebuild the Appendix C artefacts from the stored subnetwork chains (no refit).
SUPP_SKIP_FIT=1 Rscript R/supplementary.R
