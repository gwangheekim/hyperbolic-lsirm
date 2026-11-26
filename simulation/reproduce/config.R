## config.R — SINGLE SOURCE OF TRUTH for the tree-simulation pipeline.
## Every canonical script sources this; change a setting here and re-run, never edit settings elsewhere.
## (PRES must be defined by the caller before sourcing.)

.env_positive_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum || as.character(value) != raw) {
    stop(sprintf("%s must be an integer >= %d (received '%s')", name, minimum, raw))
  }
  value
}
.env_flag <- function(name, default) {
  raw <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (raw %in% c("1", "true", "yes")) {
    return(TRUE)
  }
  if (raw %in% c("0", "false", "no")) {
    return(FALSE)
  }
  stop(sprintf("%s must be true/false, yes/no, or 1/0 (received '%s')", name, raw))
}

CFG <- list(
  ## ---- DGP (tree-distance generator; see gen_binary in lib_sim.R) ----
  DEPTH = 4L, # tree depth
  BF = 2L, # branching factor  -> depth4 BF2 = 30 items
  N = 1500L, # respondents (persons), placed uniformly over depth
  GAMMA = 1.5, # DGP distance weight in logit P = theta - gamma*d(p,i) + mu
  SPREAD = 2, # SD of person ability theta_k ~ N(0, SPREAD^2); beta omitted by design

  ## ---- fitting (our MCMC engine, R/setup.R::lsirm1pl_cpp) ----
  ## At least 60 iterations are required because short pilots use nburn=niter/2
  ## and the production adaptation interval is 25 (< nburn).
  NITER = .env_positive_integer("SIM_NITER", 50000L, 60L),
  NBURN = 10000L, # burn-in (capped at niter/2 for short investigational runs)
  NTHIN = 10L, # 40,000 post-burn iterations / 10 = 4,000 retained production draws
  PR_SD_BETA = 1, # item-easiness prior SD (fixed; estimate_beta_sd=FALSE)
  SUB = 120L, # per-draw distance subsample for distortion/decoder

  ## ---- THE MODEL (gamma fixed to one) ----
  ## logit P(Y_ij) = theta_i + beta_j - d_H(z_i, w_j).  gamma is EXCLUDED (fixed = 1): the latent scale,
  ## not an interaction weight, carries the magnitude. sigma_z, sigma_w are ESTIMATED with a
  ## HALF-NORMAL(0, PR_SD_Z) prior (Metropolis). theta-sd estimated (inverse-gamma);
  ## beta-sd fixed = PR_SD_BETA. Alternative gamma-estimating modes are not part
  ## of this submission workflow.
  PARAM = "fixgamma",
  PR_SD_Z = 2.5, # half-Normal prior SCALE on sigma_z, sigma_w (also the sigma init); weakly informative
  ## The intrinsic-coordinate proposal is algebraically equivalent near the
  ## origin and avoids Lorentz cancellation for peripheral positions.
  STABLE_INTRINSIC_PROPOSAL = .env_flag("STABLE_INTRINSIC_PROPOSAL", TRUE),
  MAX_HYPERBOLIC_RADIUS = as.numeric(Sys.getenv("MAX_HYPERBOLIC_RADIUS", "13")),
  ## Coupled log-radius/main-effect moves improve traversal of the radial ridge
  ## while leaving the likelihood and wrapped-normal target unchanged.
  RADIAL_ASIS = .env_flag("RADIAL_ASIS", TRUE),
  RADIAL_ASIS_JUMP = as.numeric(Sys.getenv("RADIAL_ASIS_JUMP", ".15")),
  ## Polar split replaces the joint tangent move with a fixed-radius angular
  ## move, followed by the coupled radial move (Section 4.2). The reported
  ## Section 5 results were generated with TRUE; POLAR_SPLIT=false selects the
  ## joint-tangent sensitivity schedule.
  POLAR_SPLIT = .env_flag("POLAR_SPLIT", TRUE),

  ## ---- embeddings for the distortion table: H^2 and R^2..R^5 ----
  ## tag -> (hyp, mdim = manifold dimension); engine ndim = mdim+1 (hyp) or mdim (euclid)
  EMBED = list(
    H2 = list(hyp = TRUE, mdim = 2L),
    H3 = list(hyp = TRUE, mdim = 3L), # dimension-scaling check: does H, unlike R, need no extra room?
    R2 = list(hyp = FALSE, mdim = 2L),
    R3 = list(hyp = FALSE, mdim = 3L),
    R4 = list(hyp = FALSE, mdim = 4L),
    R5 = list(hyp = FALSE, mdim = 5L)
  ),

  ## SEEDS: env-overridable so the whole pipeline can be re-run at replicate strength without editing
  ## this file (SIM_SEEDS=100 -> 1:100). Everything that iterates CFG$SEEDS -- fit_sim.R,
  ## fit_sim_betasd05/2.R and the table builders that read fits/ -- follows automatically.
  SEEDS = seq_len(.env_positive_integer("SIM_SEEDS", 100L))
)
rm(.env_positive_integer, .env_flag)
if (!is.finite(CFG$MAX_HYPERBOLIC_RADIUS) || CFG$MAX_HYPERBOLIC_RADIUS <= 0) {
  stop("MAX_HYPERBOLIC_RADIUS must be finite and positive")
}
if (!is.finite(CFG$RADIAL_ASIS_JUMP) || CFG$RADIAL_ASIS_JUMP <= 0) {
  stop("RADIAL_ASIS_JUMP must be finite and positive")
}
if (length(CFG$NTHIN) != 1L || !is.finite(CFG$NTHIN) || CFG$NTHIN < 1L) {
  stop("NTHIN must be a positive integer")
}
CFG$FITS <- file.path(PRES, "reproduce", "fits") # where W-draws are saved
WORKSPACE <- Sys.getenv("LSIRM_WORKSPACE", unset = normalizePath(file.path(PRES, ".."), mustWork = TRUE))
CFG$FIG <- file.path(WORKSPACE, "figure")
CFG$TABLE <- file.path(WORKSPACE, "table")
dir.create(CFG$FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$TABLE, recursive = TRUE, showWarnings = FALSE)
