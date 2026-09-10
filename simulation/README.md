# Simulation pipeline

The paper-facing simulation profile contains four components:

1. H2/H3 and R2--R5 fits on the baseline depth-4 binary tree.
2. Nine branching-tree settings, including two deterministic unbalanced trees.
3. Five respondent sample sizes.
4. Two nonbranching controls.

Run `bash run_full.sh` to regenerate all components. Reported 100-seed summaries
are stored in `../result/simulation/reported`; newly generated summaries are
stored in `reproduce/results` and copied to `../result/simulation/regenerated`
after validation. Figures are written to `../figure` and LaTex table fragments
to `../table`.

`cpp/lsirmhyper.cpp` and `cpp/hyperboloid_utils.h` are the Section 5 engine and
header. `preflight.R` recompiles and checks them before any fit runs.

The production wrapper defaults to the stable intrinsic-coordinate joint
tangent proposal and radius cap 13. `STABLE_INTRINSIC_PROPOSAL=false` selects
the legacy numerical realization of the same statistical proposal, and
`MAX_HYPERBOLIC_RADIUS=15` (or another positive value) requests an explicit
cap-sensitivity run. Inspect cap-hit rates and posterior radial tails before
treating a larger cap as a production setting.

`RADIAL_ASIS=true` enables the joint radius/main-effect update and is the
production default; `RADIAL_ASIS_JUMP` controls its initial log-radius step
(default `.15`). The Section 5 fits use this augmentation together with the
joint tangent proposal. They do not use the polar-split angular replacement.
