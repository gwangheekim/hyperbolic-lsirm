## sim_robustness.R
## ---------------------------------------------------------------------------
## Robustness of the H-vs-E distortion advantage across tree families (Section 5,
## "Robustness across tree families"). Refits H^2 and R^2 on nine generators:
## balanced branching/noise variants, two deterministic 30-node imbalances
## (shallow branch size and mixed branch depth), and two random-growth trees.
## It reports multiplicative distortion over the
## configured response seeds. The baseline row is reused from fit_sim.R.
##
## Reads only lib_sim/config; writes ONLY reproduce/results/results_robustness.csv.
## ---------------------------------------------------------------------------
PRES <- Sys.getenv("LSIRM_PRES") # override: export LSIRM_PRES=/path/to/presentation
source(file.path(PRES, "reproduce/lib_sim.R"))
suppressMessages(library(parallel))
## canonical strength, matches fit_sim.R (was NSEED=3L, NITER=12000L -- weaker than the headline
## table this gets compared against; also a 3-seed SE estimate is fragile, e.g. one row previously
## rounded to R2_se=0.000). mclapply over seeds keeps wall time down at the higher NITER.
NSEED <- as.integer(Sys.getenv("ROB_NSEED", "100"))
NITER <- as.integer(Sys.getenv("ROB_NITER", as.character(CFG$NITER)))

## graph-distance matrix for an arbitrary parent/depth vector (as in build_tree)
gdmat <- function(par, dep) {
  NN <- length(par)
  gd <- function(u, v) {
    du <- dep[u]
    dv <- dep[v]
    a <- u
    c <- v
    while (du > dv) {
      a <- par[a]
      du <- du - 1
    }
    while (dv > du) {
      c <- par[c]
      dv <- dv - 1
    }
    while (a != c) {
      a <- par[a]
      c <- par[c]
      du <- du - 1
    }
    (dep[u] - du) + (dep[v] - du)
  }
  outer(1:NN, 1:NN, Vectorize(gd))
}
## unbalanced random-growth tree: each new node attaches to a random node of depth < maxdepth
build_unbal <- function(nitems = 30, maxdepth = 5, seed = 1) {
  set.seed(seed)
  par <- c(0)
  dep <- c(0)
  while (length(par) - 1 < nitems) {
    cand <- which(dep < maxdepth)
    p <- cand[sample(length(cand), 1)]
    par <- c(par, p)
    dep <- c(dep, dep[p] + 1)
  }
  items <- which(dep > 0)
  list(
    par = par, dep = dep, items = items, depi = dep[items], G = gdmat(par, dep),
    tp = sapply(items, function(v) {
      p <- par[v]
      if (p == 1) 0L else match(p, items)
    })
  )
}
## item-item distances via per_draw_dist() (R/setup.R): distance is formed PER DRAW and
## averaged, the manifold-aware estimator used by the main table (table_distortion.R). We do NOT
## average the ambient Lorentz coordinates first: on the hyperboloid the arithmetic mean lands OFF
## the manifold and shrinks distances asymmetrically (Euclidean, being flat, is unaffected), which
## would inflate the H-vs-E gap in hyperbolic's favour. Using the same estimator as the headline
## table keeps the robustness numbers directly comparable to Table~\ref{tab:distortion}.

run <- function(t, spread, lab) {
  items <- t$items
  Dt <- t$G[items, items]
  m <- mc_cols(1:NSEED, function(s) {
    Y <- gen_binary(t, N = CFG$N, seed = s, spread = spread)
    fH <- fit_embed(Y, TRUE, 2L, NITER, "fixgamma")
    fE <- fit_embed(Y, FALSE, 2L, NITER, "fixgamma")
    c(
      dH = distortion(per_draw_dist(fH$w, TRUE, CFG$SUB), Dt),
      dE = distortion(per_draw_dist(fE$w, FALSE, CFG$SUB), Dt)
    )
  }, cores = sim_cores(NSEED), lab = lab)
  dH <- m["dH", ]
  dE <- m["dE", ]
  se <- function(x) sd(x) / sqrt(length(x))
  cat(sprintf(
    "  %-30s M=%-3d  H2=%.3f(%.3f)  R2=%.3f(%.3f)  gap(R-H)=%+.3f\n",
    lab, length(items), mean(dH), se(dH), mean(dE), se(dE), mean(dE) - mean(dH)
  ))
  g <- dE - dH
  data.frame(
    family = lab, M = length(items), truth_delta = 0,
    H2 = mean(dH), H2_se = se(dH), R2 = mean(dE), R2_se = se(dE),
    gap = mean(g), gap_se = se(g), nseed = length(dH), niter = NITER,
, stringsAsFactors = FALSE
  )
}

cat(sprintf("== H^2 vs R^2 distortion across tree families (mean(SE) over %d seeds, NITER=%d) ==\n", NSEED, NITER))
base <- read_baseline_pair()
se <- function(x) sd(x) / sqrt(length(x))
baseline <- data.frame(
  family = "balanced BF2 depth4 (baseline)", M = 30L, truth_delta = 0,
  H2 = mean(base$H2), H2_se = se(base$H2), R2 = mean(base$R2), R2_se = se(base$R2),
  gap = mean(base$gap), gap_se = se(base$gap), nseed = nrow(base), niter = CFG$NITER,
, stringsAsFactors = FALSE
)
res <- rbind(
  baseline,
  run(build_tree(DEPTH = 3, BF = 3), 2, "balanced BF3 depth3"),
  run(build_tree(DEPTH = 2, BF = 5), 2, "balanced BF5 depth2 (bushy)"),
  run(build_shallow_unbalanced(30, 5), 2, "unbalanced shallow 5/23"),
  run(build_depth_unbalanced(30, c(2, 3, 4, 6)), 2, "unbalanced mixed depth 2/3/4/6"),
  run(build_tree(DEPTH = 4, BF = 2), 1, "balanced BF2 row SD 1"),
  run(build_tree(DEPTH = 4, BF = 2), 4, "balanced BF2 row SD 4"),
  run(build_unbal(30, 5, 7), 2, "unbalanced random"),
  run(build_unbal(30, 6, 11), 2, "unbalanced random deeper")
)
write.csv(res, file.path(PRES, "reproduce/results/results_robustness.csv"), row.names = FALSE)
cat("\nsaved reproduce/results/results_robustness.csv\nDONE sim_robustness\n")
