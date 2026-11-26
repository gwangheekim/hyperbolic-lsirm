# Simulation Pipeline

The simulation reported in the manuscript consists of four components:

1. H2, H3, and R2–R5 fits to the baseline binary tree of depth 4.
2. Nine branching-tree settings, including two deterministic unbalanced trees.
3. Five respondent sample sizes.
4. Two nonbranching controls.

Running `bash run_full.sh` regenerates all four components. The reported 100-seed summaries are stored in `../result/simulation/reported`. Newly generated summaries are written to `reproduce/results` and, after validation, copied to `../result/simulation/regenerated`. Figures are written to `../figure`, and LaTeX table fragments are written to `../table`.

`cpp/lsirmhyper.cpp` and `cpp/hyperboloid_utils.h` contain the engine and header used for Section 5. `preflight.R` recompiles and checks both files before any model is fitted.

The production wrapper uses the numerically stable intrinsic-coordinate proposal and caps the hyperbolic radius at 13. Setting `STABLE_INTRINSIC_PROPOSAL=false` selects the legacy numerical implementation of the same statistical proposal. Setting `MAX_HYPERBOLIC_RADIUS=15`, or any other positive value, requests an explicit sensitivity analysis of the radius cap. Cap-hit rates and the posterior radial tails should be inspected before a larger cap is adopted as a production setting.

`RADIAL_ASIS=true`, the production default, enables the coupled update of the log radius and the main effect. `RADIAL_ASIS_JUMP` sets the initial step size on the log-radius scale (default `.15`). `POLAR_SPLIT=true`, which is also the default, replaces the joint tangent proposal with a fixed-radius angular proposal, so that each hyperbolic position update consists of an angular move followed by the coupled radial move. This polar-split schedule is the sampler described in Section 4.2 and was used for the reported Section 5 fits. It requires `RADIAL_ASIS=true`. Setting `POLAR_SPLIT=false` runs the joint tangent schedule as a sensitivity check.
