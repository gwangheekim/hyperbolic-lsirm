## Validate the paper-facing simulation profile: settings, shapes, and
## arithmetic consistency, never whether a result has the desired direction.
PRES <- Sys.getenv("LSIRM_PRES")
if (!nzchar(PRES)) PRES <- normalizePath(".")
nseed <- as.integer(Sys.getenv("SIM_NSEEDS", "100"))
expected_niter <- max(60L, as.integer(Sys.getenv("SIM_NITER", "50000")))
stopifnot(length(nseed) == 1L, is.finite(nseed), nseed > 1L)
fits <- file.path(PRES, "reproduce", "fits")
results <- file.path(PRES, "reproduce", "results")

cfg <- readRDS(file.path(fits, "sim_config.rds"))
env_flag <- function(name, default) {
  raw <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (raw %in% c("1", "true", "yes")) {
    return(TRUE)
  }
  if (raw %in% c("0", "false", "no")) {
    return(FALSE)
  }
  stop(sprintf("invalid %s value: %s", name, raw))
}
stopifnot(
  identical(as.integer(cfg$SEEDS), seq_len(nseed)),
  identical(cfg$PARAM, "fixgamma"), identical(cfg$NITER, expected_niter),
  identical(cfg$NBURN, 10000L), identical(cfg$NTHIN, 10L),
  identical(cfg$PR_SD_BETA, 1),
  identical(cfg$PR_SD_Z, 2.5),
  identical(
    cfg$STABLE_INTRINSIC_PROPOSAL,
    env_flag("STABLE_INTRINSIC_PROPOSAL", TRUE)
  ),
  identical(
    cfg$MAX_HYPERBOLIC_RADIUS,
    as.numeric(Sys.getenv("MAX_HYPERBOLIC_RADIUS", "13"))
  ),
  identical(
    cfg$RADIAL_ASIS,
    env_flag("RADIAL_ASIS", TRUE)
  ),
  identical(
    cfg$RADIAL_ASIS_JUMP,
    as.numeric(Sys.getenv("RADIAL_ASIS_JUMP", ".15"))
  ),
  identical(
    cfg$POLAR_SPLIT,
    env_flag("POLAR_SPLIT", TRUE)
  )
)
rm(env_flag)

rd <- function(name) read.csv(file.path(results, name), stringsAsFactors = FALSE)
raw <- rd("results_distortion_raw.csv")
dist <- rd("results_distortion.csv")
rob <- rd("results_robustness.csv")
ns <- rd("results_n_sweep.csv")
neg <- rd("results_negative_control.csv")

emb <- c("H2", "H3", "R2", "R3", "R4", "R5")
first_saved <- ceiling(cfg$NBURN / cfg$NTHIN) * cfg$NTHIN
expected_ndraw <- length(seq.int(first_saved, cfg$NITER - 1L, by = cfg$NTHIN))
for (tag in emb) {
  first_fit <- readRDS(file.path(fits, sprintf("sim_s1_%s.rds", tag)))
  W <- if (is.list(first_fit) && !is.null(first_fit$w)) first_fit$w else first_fit
  stopifnot(
    length(dim(W)) == 3L, dim(W)[1L] == expected_ndraw,
    dim(W)[2L] == 30L,
    dim(W)[3L] == cfg$EMBED[[tag]]$mdim + as.integer(cfg$EMBED[[tag]]$hyp)
  )
}
stopifnot(
  nrow(raw) == length(emb) * nseed,
  nrow(unique(raw[c("seed", "embed")])) == length(emb) * nseed,
  identical(sort(unique(raw$embed)), sort(emb)),
  identical(sort(unique(as.integer(raw$seed))), seq_len(nseed)),
  nrow(dist) == length(emb), identical(sort(dist$embed), sort(emb)),
  nrow(rob) == 9L, nrow(ns) == 5L, nrow(neg) == 2L,
  identical(sort(as.integer(ns$N)), c(500L, 1000L, 1500L, 2000L, 4000L)),
  identical(sort(neg$truth), sort(c("flat 2D plane", "flat 1D line")))
)

for (x in list(raw, dist, rob, ns, neg)) {
  stopifnot(all(x$niter == expected_niter))
  z <- unlist(x[vapply(x, is.numeric, logical(1))], use.names = FALSE)
  if (any(!is.finite(z[!is.na(z)]))) stop("non-finite numeric result")
}
tol <- 1e-12
stopifnot(
  max(abs(rob$gap - (rob$R2 - rob$H2))) < tol,
  max(abs(ns$gap - (ns$dist_r - ns$dist_h))) < tol,
  max(abs(neg$gap - (neg$R2 - neg$H2))) < tol,
  max(abs(neg$gap_add - (neg$R2_add - neg$H2_add))) < tol
)

## The rows described as the common baseline must be numerically identical to
## the canonical per-seed fit, not an independent refit with another RNG stream.
h <- raw[raw$embed == "H2", c("seed", "distortion")]
r <- raw[raw$embed == "R2", c("seed", "distortion")]
b <- merge(h, r, by = "seed", suffixes = c("_H", "_R"))
bh <- mean(b$distortion_H)
br <- mean(b$distortion_R)
bg <- mean(b$distortion_R - b$distortion_H)
rb <- rob[rob$family == "balanced BF2 depth4 (baseline)", ]
nb <- ns[ns$N == 1500L, ]
stopifnot(
  nrow(rb) == 1L, nrow(nb) == 1L,
  max(abs(c(
    rb$H2 - bh, rb$R2 - br, rb$gap - bg,
    nb$dist_h - bh, nb$dist_r - br, nb$gap - bg
  ))) < tol
)

cat(sprintf(
  "RESULT VALIDATION PASS | %d seeds | MCMC %d/%d/%d | %d retained\n",
  nseed, expected_niter, cfg$NBURN, cfg$NTHIN, expected_ndraw
))
