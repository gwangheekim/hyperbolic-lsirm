## sim_negative_control.R — NEGATIVE CONTROL for the H-vs-E distortion advantage.
##
## The tree DGP shows H beats E when the truth is a hierarchy. It cannot show the
## comparison is DISCRIMINATING rather than hyperbolic-biased. This script generates
## from NON-tree truths (2D Euclidean plane, 1D line) with everything else held fixed,
## and checks that the hyperbolic fit does NOT beat Euclidean there.
##
## Design: identical to the tree DGP (same N, gamma, spread, mu-centering, H^2/R^2 fits,
## per-draw distortion + optimal c*) — the change is the ground-truth metric.
## The flat truths are rescaled so their pairwise-distance SD matches the tree's, so a
## distortion difference reflects the TRUTH's geometry, not signal strength.
## Output: console table + reproduce/results/results_negative_control.csv (writes nothing to figure/).
##
##   cd reproduce && LSIRM_PRES=.. Rscript sim_negative_control.R   (or via env NEG_SEEDS / NEG_NITER)

PRES <- Sys.getenv("LSIRM_PRES")
source(file.path(PRES, "reproduce/lib_sim.R")) # fit_embed, per_draw_dist, distortion, getW, build_tree, gen_binary, CFG
suppressMessages(library(parallel))

SEEDS <- as.integer(Sys.getenv("NEG_SEEDS", "100"))
NITER <- as.integer(Sys.getenv("NEG_NITER", as.character(CFG$NITER)))
SUB <- CFG$SUB

## signal-strength anchor: tree Dtrue off-diagonal SD (all flat truths matched to this)
tt <- build_tree()
tg <- tt$G[tt$items, tt$items]
M <- length(tt$items)
sd_tree <- sd(tg[upper.tri(tg)])

## flat Euclidean truth (dim=2 plane, dim=1 line), signal-matched to the tree
build_flat <- function(M, dim, seed) {
  set.seed(1000L + seed)
  pos <- matrix(runif(M * dim), M, dim)
  D <- as.matrix(dist(pos))
  D <- D * (sd_tree / sd(D[upper.tri(D)])) # match tree pairwise-distance SD
  list(items = 1:M, G = D)
}
## responses from a flat truth — identical form to gen_binary (persons seated at item positions)
gen_flat <- function(t, N = CFG$N, seed, gamma = CFG$GAMMA, spread = CFG$SPREAD,
                     seat_prob = NULL) {
  M <- length(t$items)
  set.seed(seed)
  pn <- sample(t$items, N, replace = TRUE, prob = seat_prob)
  Dpi <- t$G[pn, t$items]
  th <- rnorm(N, 0, spread)
  mu <- mean(gamma * Dpi)
  Y <- matrix(rbinom(N * M, 1, 1 / (1 + exp(-(th - gamma * Dpi + mu)))), N, M)
  storage.mode(Y) <- "double"
  Y
}

## Match the deep tree's equal-per-depth row placement in the nonbranching
## controls. Item labels carry the same seat probabilities; only their metric
## changes. Flat layouts are regenerated across seeds, so no single assignment
## of these weights to flat coordinates is privileged.
deep_prob <- ave(rep(1, M), tt$depi, FUN = function(x) 1 / length(x))
deep_prob <- deep_prob / sum(deep_prob)

## one nonbranching truth condition: H^2 and R^2 across seeds, distortion vs its OWN truth
additive_error <- function(Dhat, Dtrue) {
  ut <- upper.tri(Dtrue)
  d <- Dhat[ut]
  g <- Dtrue[ut]
  o <- g > 0 & is.finite(g) & is.finite(d) & d >= 0
  den <- sum(d[o]^2)
  if (!any(o) || den <= 0) stop("additive error is undefined")
  cc <- sum(d[o] * g[o]) / den
  mean(abs(cc * d[o] - g[o])) / mean(g[o])
}
run_cond <- function(lab, truth_fn, gen_fn) {
  m <- mc_cols(1:SEEDS, function(s) {
    t <- truth_fn(s)
    Dt <- t$G[t$items, t$items]
    Y <- gen_fn(t, seed = s)
    DH <- per_draw_dist(getW(fit_embed(Y, TRUE, 2L, NITER, "fixgamma")), TRUE, SUB)
    DE <- per_draw_dist(getW(fit_embed(Y, FALSE, 2L, NITER, "fixgamma")), FALSE, SUB)
    v <- c(
      H2 = distortion(DH, Dt), R2 = distortion(DE, Dt),
      H2_add = additive_error(DH, Dt), R2_add = additive_error(DE, Dt),
      truth_delta = delta_over_diam_exact(Dt)
    )
    cat(sprintf(
      "  [%s seed %d] H2=%.3f R2=%.3f gap(R-H)=%+.3f\n",
      lab, s, v["H2"], v["R2"], v["R2"] - v["H2"]
    ))
    v
  }, cores = sim_cores(SEEDS), lab = lab)
  n_ok <- ncol(m)
  g <- m["R2", ] - m["H2", ]
  ga <- m["R2_add", ] - m["H2_add", ]
  data.frame(
    truth = lab, H2 = mean(m["H2", ]), R2 = mean(m["R2", ]),
    gap = mean(g), gap_se = sd(g) / sqrt(n_ok),
    H2_add = mean(m["H2_add", ]), R2_add = mean(m["R2_add", ]),
    gap_add = mean(ga), gap_add_se = sd(ga) / sqrt(n_ok),
    truth_delta = mean(m["truth_delta", ]), truth_delta_se = sd(m["truth_delta", ]) / sqrt(n_ok),
    nseed = n_ok, niter = NITER,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("== NEGATIVE CONTROL (mean over %d seeds, NITER=%d, N=%d, M=%d) ==\n", SEEDS, NITER, CFG$N, M))
cat("gap = distortion(R2) - distortion(H2); POSITIVE => hyperbolic more faithful.\n")
cat("Expectation: negative gaps for nonbranching truths.\n\n")
res <- rbind(
  run_cond(
    "flat 2D plane", function(s) build_flat(M, 2, s),
    function(t, seed) gen_flat(t, seed = seed, seat_prob = deep_prob)
  ),
  run_cond(
    "flat 1D line", function(s) build_flat(M, 1, s),
    function(t, seed) gen_flat(t, seed = seed, seat_prob = deep_prob)
  )
)
cat("\n")
print(format(res, digits = 3), row.names = FALSE)
write.csv(res, file.path(PRES, "reproduce/results/results_negative_control.csv"), row.names = FALSE)
cat("\nsaved reproduce/results/results_negative_control.csv\n")
