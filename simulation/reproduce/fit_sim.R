## fit_sim.R — fit the tree DGP in H^2, H^3 and R^2..R^5 for every seed.
## By default only item-position draws are retained; they are all the summary builders use.
## ALL settings come from config.R. Readers use getW() so both full-fit and bare-W files work.
## Outputs: fits/sim_truth.rds, fits/sim_config.rds, fits/sim_s{seed}_{tag}.rds  (tag in H2,R2,R3,R4,R5).
PRES <- Sys.getenv("LSIRM_PRES") # override: export LSIRM_PRES=/path/to/presentation
source(file.path(PRES, "reproduce/lib_sim.R"))
suppressMessages(library(parallel))
dir.create(CFG$FITS, showWarnings = FALSE, recursive = TRUE)

t <- build_tree() # CFG$DEPTH, CFG$BF
saveRDS(list(par = t$par, dep = t$dep, items = t$items, depi = t$depi, tp = t$tp), file.path(CFG$FITS, "sim_truth.rds"))
saveRDS(CFG, file.path(CFG$FITS, "sim_config.rds")) # provenance: what produced these fits
cat(sprintf(
  paste0(
    "tree depth%d BF%d: %d items | N=%d | MCMC=%d/%d/%d (iter/burn/thin) | ",
    "PARAM=%s (pr_sd_z=%g) | proposal=%s, radial-ASIS=%s | embeddings: %s\n"
  ),
  CFG$DEPTH, CFG$BF, length(t$items), CFG$N, CFG$NITER, CFG$NBURN, CFG$NTHIN, CFG$PARAM, CFG$PR_SD_Z,
  if (CFG$STABLE_INTRINSIC_PROPOSAL) "stable" else "legacy",
  if (CFG$RADIAL_ASIS) "on" else "off",
  paste(names(CFG$EMBED), collapse = ",")
))

jobs <- mclapply(CFG$SEEDS, function(s) {
  Y <- gen_binary(t, seed = s)
  for (tag in names(CFG$EMBED)) {
    e <- CFG$EMBED[[tag]]
    .f <- fit_embed(Y, e$hyp, e$mdim, CFG$NITER, CFG$PARAM)
    ## SAVE_W_ONLY=1 (default): store only the column draws.  Every reader goes through getW(), which
    ## accepts a bare-W array, and the row draws z are ~98% of the file (4000 x N x d vs 4000 x M x d).
    ## At 100 seeds the full objects would be about 160 GB.
    if (Sys.getenv("SAVE_W_ONLY", "1") == "1") .f <- .f$w
    saveRDS(
      .f,
      file.path(CFG$FITS, sprintf("sim_s%d_%s.rds", s, tag))
    )
  }
  s
}, mc.cores = sim_cores(length(CFG$SEEDS)), mc.preschedule = FALSE)
bad <- vapply(jobs, inherits, logical(1), "try-error")
if (any(bad)) {
  stop(sprintf(
    "fit_sim failed for seeds: %s",
    paste(CFG$SEEDS[bad], collapse = ",")
  ))
}
cat(sprintf(
  "ALL FITS DONE  (%d seeds x %d embeddings, MCMC=%d/%d/%d, PARAM=%s)\n",
  length(CFG$SEEDS), length(CFG$EMBED), CFG$NITER, CFG$NBURN, CFG$NTHIN, CFG$PARAM
))
