## Fast server preflight: dependencies, parsing, engine compile, DGP invariants,
## and a tiny H^2/E^2 end-to-end fit under the production model settings.

PRES <- Sys.getenv("LSIRM_PRES")
if (!nzchar(PRES)) stop("LSIRM_PRES is not set")
required <- c("Rcpp", "RcppArmadillo", "igraph")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("missing R packages: ", paste(missing, collapse = ", "))

## The submission package contains only the paper-facing Section 5 profile;
## optional exploratory and sensitivity generators remain in the development
## workspace and must not be required by this self-contained preflight.
scheduled <- c(
  "fit_sim.R", "sim_n_sweep.R", "sim_negative_control.R",
  "sim_robustness.R", "table_distortion.R", "table_robustness.R",
  "validate_results.R"
)
for (f in scheduled) parse(file.path(PRES, "reproduce", f))

source(file.path(PRES, "reproduce", "lib_sim.R"))
## Guard the R-side geometry against regression to the cancellation-prone
## acosh(-<x,y>_L) formula at the sampler's peripheral radius.
r_test <- 13
angle_test <- 1e-8
P_test <- rbind(
  c(cosh(r_test), sinh(r_test), 0),
  c(
    cosh(r_test), sinh(r_test) * cos(angle_test),
    sinh(r_test) * sin(angle_test)
  )
)
d_expected <- 2 * asinh(sinh(r_test) * sin(angle_test / 2))
stopifnot(
  abs(hyperbolic_dist(P_test)[1, 2] - d_expected) < 1e-12,
  max(abs(radial_depth(P_test) - r_test)) < 1e-12
)

t <- build_tree()
stopifnot(
  length(t$items) == 30L, max(t$depi) == 4L,
  all(diag(t$G) == 0), all(t$G == t(t$G))
)
shallow <- build_shallow_unbalanced()
shallow_leaf_counts <- tabulate(shallow$branch[shallow$dep == 2L], nbins = 2L)
stopifnot(
  length(shallow$items) == 30L, max(shallow$depi) == 2L,
  identical(shallow_leaf_counts, c(5L, 23L)), all(shallow$G == t(shallow$G))
)
unequal_depth <- build_depth_unbalanced()
observed_branch_depths <- vapply(
  seq_along(unequal_depth$branch_depths),
  function(b) max(unequal_depth$dep[unequal_depth$branch == b]), integer(1)
)
stopifnot(
  length(unequal_depth$items) == 30L,
  identical(observed_branch_depths, c(2L, 3L, 4L, 6L)),
  all(unequal_depth$G == t(unequal_depth$G))
)
Y <- gen_binary(t, N = 24L, seed = 919L)
Y_shallow <- gen_binary(shallow, N = 24L, seed = 920L)
Y_depth <- gen_binary(unequal_depth, N = 24L, seed = 921L)
stopifnot(identical(dim(Y), c(24L, 30L)), all(Y %in% c(0, 1)))
stopifnot(
  identical(dim(Y_shallow), c(24L, 30L)), all(Y_shallow %in% c(0, 1)),
  identical(dim(Y_depth), c(24L, 30L)), all(Y_depth %in% c(0, 1))
)

H <- fit_embed(Y, TRUE, 2L, niter = 60L, param = "fixgamma")
E <- fit_embed(Y, FALSE, 2L, niter = 60L, param = "fixgamma")
B <- fit_embed(Y, TRUE, 2L, niter = 60L, param = "fixgamma", estimate_beta_sd = TRUE)
smoke_burn <- min(CFG$NBURN, 60L %/% 2L)
smoke_first_saved <- ceiling(smoke_burn / CFG$NTHIN) * CFG$NTHIN
smoke_saved_iterations <- seq.int(smoke_first_saved, 59L, by = CFG$NTHIN)
stopifnot(
  isTRUE(H$hyperbolic), identical(E$hyperbolic, FALSE),
  identical(CFG$NBURN, 10000L), identical(CFG$NTHIN, 10L),
  H$latent_dim == 2L, E$latent_dim == 2L,
  dim(H$w)[1L] == length(smoke_saved_iterations),
  identical(as.integer(H$retained_iteration), as.integer(smoke_saved_iterations)),
  identical(as.integer(E$retained_iteration), as.integer(smoke_saved_iterations)),
  identical(as.integer(B$retained_iteration), as.integer(smoke_saved_iterations)),
  dim(H$z)[3] == 3L, dim(E$z)[3] == 2L,
  H$fix_gamma, E$fix_gamma, H$sample_zw_sd == 2L,
  H$fix_beta_sd, E$fix_beta_sd, !H$fix_theta_sd, !E$fix_theta_sd,
  !B$fix_beta_sd, isTRUE(B$estimate_beta_sd), all(B$sigma_beta > 0),
  !H$use_gromov, !E$use_gromov,
  E$sample_zw_sd == 2L, all(H$gamma == 1), all(E$gamma == 1),
  identical(H$stable_intrinsic_proposal, CFG$STABLE_INTRINSIC_PROPOSAL),
  identical(H$max_hyperbolic_radius, CFG$MAX_HYPERBOLIC_RADIUS),
  identical(H$use_radial_asis, CFG$RADIAL_ASIS),
  !E$use_radial_asis,
  identical(H$polar_split, CFG$POLAR_SPLIT),
  !E$polar_split,
  all(H$postburn_radial_acceptance_rate_z >= 0 &
    H$postburn_radial_acceptance_rate_z <= 1),
  all(H$postburn_radial_acceptance_rate_w >= 0 &
    H$postburn_radial_acceptance_rate_w <= 1),
  all(H$postburn_radial_boundary_rate_z >= 0 &
    H$postburn_radial_boundary_rate_z <= 1),
  all(H$postburn_radial_boundary_rate_w >= 0 &
    H$postburn_radial_boundary_rate_w <= 1),
  all(H$sigma_z > 0), all(H$sigma_w > 0), all(E$sigma_z > 0), all(E$sigma_w > 0),
  is.finite(H$accept_sigma_z), is.finite(H$accept_sigma_w),
  is.finite(E$accept_sigma_z), is.finite(E$accept_sigma_w),
  all(is.finite(H$log_posterior)), all(is.finite(E$log_posterior))
)

D <- t$G[t$items, t$items]
stopifnot(abs(distortion(2 * D, D)) < 1e-12)
DH <- per_draw_dist(H$w, TRUE, 6L)
DE <- per_draw_dist(E$w, FALSE, 6L)
stopifnot(
  identical(dim(DH), c(30L, 30L)), identical(dim(DE), c(30L, 30L)),
  max(abs(DH - t(DH))) < 1e-10, max(abs(DE - t(DE))) < 1e-10,
  all(DH >= 0), all(DE >= 0), is.finite(distortion(DH, D)), is.finite(distortion(DE, D))
)
rec <- eval_recovery(H$w, list(items = t$items, depi = t$depi, tp = t$tp), "hyperbolic")$overall
stopifnot(is.finite(rec), rec >= 0, rec <= 1)
cat(sprintf(
  paste0(
    "PREFLIGHT PASS | MCMC default %d/%d/%d | ",
    "balanced/shallow/depth-unbalanced trees valid | proposal %s | ",
    "ASIS %s | cap %g | H^2/E^2 smoke test finite\n"
  ),
  CFG$NITER, CFG$NBURN, CFG$NTHIN,
  if (CFG$STABLE_INTRINSIC_PROPOSAL) "stable" else "legacy",
  if (CFG$RADIAL_ASIS) "on" else "off", CFG$MAX_HYPERBOLIC_RADIUS
))
