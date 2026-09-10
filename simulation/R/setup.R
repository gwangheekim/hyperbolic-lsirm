## setup.R — portable setup for the simulation workflow.
##
## The original presentation/R/00_setup.R hardcodes an absolute local path (PKG <- "/Volumes/...")
## at the top, which OVERRIDES any LSIRM_PRES the caller passed in, so sourcing it on a server (or any
## machine without that exact path) fails. This file provides only what the simulation pipeline needs
## (the compiled MCMC engine + per_draw_dist), and resolves its own location from LSIRM_PRES exclusively
## (no hardcoded fallback path), so the whole sim_pipeline_server/ folder is relocatable: copy it
## anywhere, set LSIRM_PRES to that location (the run_*.sh wrappers do this for you), and it works.
##
## Deliberately excludes the CDI/HPO real-data loaders in the original 00_setup.R: this package is
## simulation-only and does not ship or need that data.

PRES <- Sys.getenv("LSIRM_PRES")
if (!nzchar(PRES)) stop("LSIRM_PRES is not set. Run via run_fits.sh / run_fits_parallel.sh, or `export LSIRM_PRES=/path/to/sim_pipeline_server` first.")
SRC <- file.path(PRES, "cpp")

suppressMessages(library(Rcpp))
if (!isTRUE(getOption("ENGINE_READY"))) {
  engine_text <- readLines(file.path(SRC, "lsirmhyper.cpp"), warn = FALSE)
  header_text <- readLines(file.path(SRC, "hyperboloid_utils.h"), warn = FALSE)
  polar_helpers <- c(
    "propose_angular_fixed_radius_into",
    "geodesic_distance_cached_polar"
  )
  missing_helpers <- polar_helpers[
    vapply(polar_helpers, function(name) {
      any(grepl(name, engine_text, fixed = TRUE)) &&
        !any(grepl(name, header_text, fixed = TRUE))
    }, logical(1))
  ]
  if (length(missing_helpers)) {
    stop(
      "lsirmhyper.cpp and hyperboloid_utils.h are from incompatible versions; ",
      "copy the current top-level source and header together. Missing header helpers: ",
      paste(missing_helpers, collapse = ", ")
    )
  }
  sourceCpp(file.path(SRC, "lsirmhyper.cpp")) # provides lsirm1pl_cpp()
  options(ENGINE_READY = TRUE)
}

## ---- how many workers an mclapply here should use -------------------------
## Every script used to carry its own hardcoded cap: min(5L, ...), min(4L, ...), a bare 2L.  Those
## numbers had no stated basis, they did not agree with each other, and they are why 100 seeds went
## through 20 sequential rounds instead of using the cores that were there.  One place decides it now,
## and it decides by asking how many cores this job was actually given.
##
## In order: LSIRM_MC_CORES (what run_fits_parallel.sh hands each job), then SLURM_CPUS_PER_TASK (what
## the scheduler allocated), then the machine.  A job never asks for more workers than it has tasks.
sim_cores <- function(n_tasks, quiet = FALSE) {
  n_tasks <- max(1L, as.integer(n_tasks))
  env_int <- function(k) {
    v <- suppressWarnings(as.integer(Sys.getenv(k, "")))
    if (is.na(v) || v < 1L) NA_integer_ else v
  }
  budget <- env_int("LSIRM_MC_CORES")
  src <- "LSIRM_MC_CORES"
  if (is.na(budget)) {
    budget <- env_int("SLURM_CPUS_PER_TASK")
    src <- "SLURM_CPUS_PER_TASK"
  }
  if (is.na(budget)) {
    budget <- max(1L, parallel::detectCores())
    src <- "detectCores"
  }
  n <- as.integer(max(1L, min(n_tasks, budget)))
  if (!quiet) {
    message(sprintf("[sim_cores] %d workers for %d tasks (budget %d from %s)", n, n_tasks, budget, src))
  }
  n
}

## ---- per-draw, gauge-invariant pairwise distance (items) -------------------
## per-draw distances are on-manifold and invariant to rotation/translation/boost;
## averaging them is the correct summary (NOT the distance of mean coordinates).
per_draw_dist <- function(W, hyperbolic, ndraws = 200) {
  nd <- dim(W)[1]
  nw <- dim(W)[2]
  ds <- round(seq(1, nd, length.out = min(ndraws, nd)))
  D <- matrix(0, nw, nw)
  for (d in ds) {
    P <- W[d, , ]
    if (hyperbolic) {
      ## tree_utils.R supplies the same non-cancelling polar identity used by
      ## src/hyperboloid_utils.h; lib_sim.R loads it before this is evaluated.
      D <- D + hyperbolic_dist(P)
    } else {
      D <- D + as.matrix(dist(P))
    }
  }
  D / length(ds)
}
