## sim_n_sweep.R
## ---------------------------------------------------------------------------
## Reviewer-defense proposal #4: does the H-vs-E distortion advantage depend on the respondent
## sample size N, or is it a finite-sample artifact of the paper's default N=1500 (CFG$N)? Holds the
## SAME main-text DGP fixed (tree depth 4, BF 2, M=30 items, gamma=1.5, spread=2 -- everything in
## CFG except N) and refits H^2/R^2 at N = 500, 1000, 1500 (baseline), 2000, 4000, over 100 seeds.
## Steps are kept close together (roughly 2x end-to-end, not the earlier 300->6000 20x range) so the
## curve is read off adjacent, comparable points rather than two widely separated regimes.
## A reviewer could otherwise argue the gap is a small-N quirk that closes as more respondents let
## flat space "fill in" the tree's exponential branching -- or, in the other direction, an artifact
## that ONLY appears once N is large enough to pin down the geometry precisely. This sweep checks both.
##
## The N=1500 row reuses the canonical fit_sim.R baseline; the other N values
## are independently fitted. Writes reproduce/results/results_n_sweep.csv.
##   cd reproduce && LSIRM_PRES=.. Rscript sim_n_sweep.R
##   env overrides: NS_N (comma-separated, default 500,1000,1500,2000,4000), NS_SEEDS (default 1:5),
##                  NS_NITER (default CFG$NITER = 50000)
## ---------------------------------------------------------------------------
PRES <- Sys.getenv("LSIRM_PRES")
source(file.path(PRES, "reproduce/lib_sim.R"))
suppressMessages(library(parallel))

N_GRID <- {
  e <- Sys.getenv("NS_N", "")
  if (nzchar(e)) as.integer(strsplit(e, ",")[[1]]) else c(500L, 1000L, 1500L, 2000L, 4000L)
}
SEEDS <- {
  e <- Sys.getenv("NS_SEEDS", "")
  if (nzchar(e)) as.integer(strsplit(e, ",")[[1]]) else 1:100
}
NITER <- as.integer(Sys.getenv("NS_NITER", as.character(CFG$NITER)))

t <- build_tree() # depth 4, BF 2, M = 30 -- identical DGP tree to Section 5
tg <- t$G[t$items, t$items]

cat(sprintf(
  "== N-scaling sweep (tree depth%d BF%d M=%d, NITER=%d, seeds %s, N in {%s}) ==\n",
  CFG$DEPTH, CFG$BF, length(t$items), NITER, paste(SEEDS, collapse = ","), paste(N_GRID, collapse = ",")
))

run_n <- function(N) {
  m <- mc_cols(SEEDS, function(s) {
    Y <- gen_binary(t, N = N, seed = s)
    fH <- fit_embed(Y, TRUE, 2L, NITER, "fixgamma")
    fR <- fit_embed(Y, FALSE, 2L, NITER, "fixgamma")
    c(
      dist_h = distortion(per_draw_dist(fH$w, TRUE, CFG$SUB), tg),
      dist_r = distortion(per_draw_dist(fR$w, FALSE, CFG$SUB), tg)
    )
  }, cores = sim_cores(length(SEEDS)), lab = paste0("N=", N))
  m <- t(m)
  se <- function(x) sd(x) / sqrt(length(x))
  row <- data.frame(
    N = N, dist_h = mean(m[, "dist_h"]), dist_h_se = se(m[, "dist_h"]),
    dist_r = mean(m[, "dist_r"]), dist_r_se = se(m[, "dist_r"]),
    gap = mean(m[, "dist_r"] - m[, "dist_h"]), gap_se = se(m[, "dist_r"] - m[, "dist_h"]),
    truth_delta = 0, nseed = nrow(m), niter = NITER,
, stringsAsFactors = FALSE
  )
  cat(sprintf(
    "  N=%-5d  H2=%.3f(%.3f)  R2=%.3f(%.3f)  gap=%+.3f(%.3f)\n",
    N, row$dist_h, row$dist_h_se, row$dist_r, row$dist_r_se, row$gap, row$gap_se
  ))
  row
}

base <- read_baseline_pair()
se <- function(x) sd(x) / sqrt(length(x))
base_row <- data.frame(
  N = CFG$N, dist_h = mean(base$H2), dist_h_se = se(base$H2),
  dist_r = mean(base$R2), dist_r_se = se(base$R2), gap = mean(base$gap), gap_se = se(base$gap),
  truth_delta = 0, nseed = nrow(base), niter = CFG$NITER,
  stringsAsFactors = FALSE
)
other <- setdiff(N_GRID, CFG$N)
res <- rbind(if (length(other)) do.call(rbind, lapply(other, run_n)) else NULL, base_row)
res <- res[order(res$N), ]
write.csv(res, file.path(PRES, "reproduce/results", "results_n_sweep.csv"), row.names = FALSE)
cat("\nsaved reproduce/results/results_n_sweep.csv\nDONE sim_n_sweep\n")
