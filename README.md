# Hyperbolic LSIRM: Biometrics submission workspace

This package reproduces the analyses, tables, and figures in `lsirm_hyperbolic_biometrics.tex`.

```text
hyperboliclsirm/
├── R/                          main.R (Section 6), supplementary.R (Appendix C), util.R
│   └── figures/                dementia figure and case-table renderers
├── cpp/                        MCMC and alignment C++ engines
├── data/
│   ├── raw/                    supplied HPO and Mondo releases
│   └── analysis_data/          dementia and subnetwork response matrices
├── result/
│   ├── model_fits/             H3 chains, not tracked, rebuilt by R/main.R and R/supplementary.R
│   ├── dementia/reported/      Section 6 artefacts
│   ├── simulation/reported/    100-seed simulation summaries and Figure 1
│   └── subnetworks/            Appendix C artefacts (ranking_summaries/ not tracked)
├── simulation/                 simulation generator, fitter, and renderers
├── figure/                     the two manuscript figures
├── table/                      the three manuscript table fragments
└── scripts/                    run and validation commands
```

`result/model_fits/` (about 11 GB of MCMC chains) and `result/subnetworks/ranking_summaries/` are not under version control. `Rscript R/main.R` and `Rscript R/supplementary.R` recreate them from `data/`. The reported summaries under `result/dementia/reported/`, `result/subnetworks/reported/`, and `result/simulation/reported/` are tracked so that validation runs without a refit, and re-running the scripts overwrites them with identical content.

## Dementia application (Section 6)

```bash
Rscript R/main.R
```

`R/main.R` builds the 170 by 313 dementia response matrix from the HPO and Mondo releases (set `SKIP_MATRIX_BUILD=1` to reuse the stored one), compiles the C++ engine, fits three H3 chains (100,000 iterations, 10,000 burn-in, thinning by 10, fixed chain seeds, polar-split radial-ridge sampler), writes them to `result/model_fits/dementia/`, and then runs the Section 6 analysis. `S6_SKIP=1` fits only. Both the fit and the analysis are deterministic.

The analysis writes, under `result/dementia/reported/`:

- `dementia_f3_geometry.rds` — pooled 27,000-draw posterior-mean disease geometry (distances, radii, Gromov products).
- `rank_stability.rds`, `rank_stability_summary.csv` — cross-chain full-vocabulary top-ten stability, over all retained draws.
- `annotation_internal_disease_summary.csv`, `annotation_internal_consensus.csv`, `full_vocab_by_chain.csv` — the distance-versus-Gromov ranking comparison over the 134 eligible records.
- `application_summary.csv` — the scalar figures quoted in Section 6.
- `selected_case_table.csv` — the choreoacanthocytosis rows of Table 3.
- `s6_intermediate.rds` — the exact-profile bookkeeping shared by the stages.
- `dementia_groups.csv` — the eight-group record-level membership for Web Appendix D.

It then writes `table/table_dementia_case.tex` (Table 3) and `figure/fig_dementia_radial_tree.pdf` (Figure 2). Check the numbers against the manuscript with `Rscript scripts/check_section6_reproduction.R result/dementia/reported`.

## Supplementary Appendix C subnetworks

```bash
Rscript R/supplementary.R
```

Fits three H3 chains for each subnetwork in `scope_tags` (`R/util.R`) with the same sampler as `R/main.R` (~15 h; deterministic), then rebuilds, under `result/subnetworks/` and over all retained draws: per-tag scalar diagnostics, cross-chain rank stability (`ranking_summaries/`, `<tag>_rank_stability.csv`), the four annotation-internal case tables (`cases/<tag>.csv`), and the summary table (`appendix_c_summary.csv`). The dementia row is computed from the `R/main.R` fit. `SUPP_SKIP_FIT=1` reuses the stored chains; `APPC_RESUME=1` skips subnetworks whose outputs are already current.

## Simulation (Section 5, Tables 1 and 2, Figure 1)

```bash
bash simulation/run_full.sh
```

Generates the tree families, fits the hyperbolic and Euclidean embeddings over 100 seeds, and writes the summaries to `result/simulation/reported/`. `simulation/render_reported.R` and `simulation/render_biometrics_truths.R` produce `table/table_distortion.tex`, `table/table_robustness.tex`, and `figure/fig_sim_negcontrol_truths.pdf` from those summaries.

## Manuscript assets and validation

```bash
bash scripts/run_biometrics_assets.sh
```

Renders the simulation tables and Figure 1, the dementia Table 3 and Figure 2, and runs `scripts/validate_biometrics_submission.R`, which checks every displayed numerical result against the reported files and finishes with `BIOMETRICS SUBMISSION VALIDATION PASS`.
