# Hyperbolic LSIRM: Biometrics Submission Workspace

This package reproduces the analyses, tables, and figures reported in "Mapping Disease–Phenotype Associations with a Hyperbolic Latent Space Item Response Model."

```text
hyperboliclsirm/
├── R/                          main.R (Section 6), supplementary.R (Appendix C), util.R
│   └── figures/                renderers for the dementia figure and case table
├── cpp/                        C++ engines for MCMC sampling and alignment
├── data/
│   ├── raw/                    HPO and Mondo releases
│   └── analysis_data/          dementia and subnetwork response matrices
├── result/
│   ├── model_fits/             H3 chains (untracked; rebuilt by R/main.R and R/supplementary.R)
│   ├── dementia/reported/      Section 6 outputs
│   ├── simulation/reported/    100-seed simulation summaries and Figure 1
│   └── subnetworks/            Appendix C outputs (ranking_summaries/ is untracked)
├── simulation/                 simulation generators, fitting scripts, and renderers
├── figure/                     the two manuscript figures
├── table/                      the three manuscript table fragments
└── scripts/                    execution and validation scripts
```

The directories `result/model_fits/` (approximately 2.2 GB of MCMC chains) and `result/subnetworks/ranking_summaries/` are excluded from version control. Running `Rscript R/main.R` and `Rscript R/supplementary.R` recreates them from the contents of `data/`. The reported summaries in `result/dementia/reported/`, `result/subnetworks/reported/`, and `result/simulation/reported/` are tracked so that validation can be performed without refitting the models. Re-running `R/main.R` and `R/supplementary.R` overwrites the dementia and subnetwork summaries with identical content, whereas the simulation pipeline writes its regenerated summaries to a separate directory.

## Dementia Application (Section 6)

```bash
Rscript R/main.R
```

`R/main.R` first builds the 170 × 313 dementia response matrix from the HPO and Mondo releases; setting `SKIP_MATRIX_BUILD=1` reuses the stored matrix instead. The script then compiles the C++ engine and fits three H3 chains with 100,000 iterations, 10,000 burn-in iterations, thinning by 10, fixed chain seeds, and the polar-split sampler described in Section 4.2. The chains are written to `result/model_fits/dementia/`, after which the Section 6 analysis is run. Setting `S6_SKIP=1` performs the model fit without the analysis. Both the fit and the analysis are deterministic.

The analysis writes the following files to `result/dementia/reported/`:

- `dementia_f3_geometry.rds`: the posterior-mean disease geometry (distances, radii, and Gromov products) pooled over 27,000 draws.
- `rank_stability.rds` and `rank_stability_summary.csv`: the cross-chain stability of the full-vocabulary top-ten rankings, computed over all retained draws.
- `annotation_internal_disease_summary.csv`, `annotation_internal_consensus.csv`, and `full_vocab_by_chain.csv`: the comparison of distance and Gromov-product rankings over the 134 eligible disease records.
- `application_summary.csv`: the scalar values quoted in Section 6.
- `selected_case_table.csv`: the choreoacanthocytosis rows of Table 3.
- `s6_intermediate.rds`: the exact-profile bookkeeping shared across the analysis stages.
- `dementia_groups.csv`: the record-level membership of the eight disease groups reported in Web Appendix D.

The script then writes `table/table_dementia_case.tex` (Table 3) and `figure/fig_dementia_radial_tree.pdf` (Figure 2). The reproduced values can be checked against the manuscript with `Rscript scripts/check_section6_reproduction.R result/dementia/reported`.

## Supplementary Subnetworks (Web Appendix C)

```bash
Rscript R/supplementary.R
```

`R/supplementary.R` fits three H3 chains for each subnetwork listed in `scope_tags` (defined in `R/util.R`) with the same sampler as `R/main.R`. The fits take approximately 15 hours and are deterministic. Using all retained draws, the script then rebuilds the following outputs in `result/subnetworks/`: the scalar diagnostics for each subnetwork, the cross-chain rank stability (`ranking_summaries/` and `<tag>_rank_stability.csv`), the four annotation-internal case tables (`cases/<tag>.csv`), and the summary table (`appendix_c_summary.csv`). The dementia row is computed from the `R/main.R` fit. Setting `SUPP_SKIP_FIT=1` reuses the stored chains, and setting `APPC_RESUME=1` skips subnetworks whose outputs are already up to date.

## Simulation (Section 5, Tables 1 and 2, Figure 1)

```bash
bash simulation/run_full.sh
```

`simulation/run_full.sh` generates the tree families, fits the hyperbolic and Euclidean embeddings over 100 seeds, and validates the resulting summaries. The summaries are written to `simulation/reproduce/results/` and copied to `result/simulation/regenerated/`, so the reported files in `result/simulation/reported/` are not overwritten. `simulation/render_reported.R` and `simulation/render_biometrics_truths.R` produce `table/table_distortion.tex`, `table/table_robustness.tex`, and `figure/fig_sim_negcontrol_truths.pdf` from the reported summaries. Further details are given in `simulation/README.md`.
