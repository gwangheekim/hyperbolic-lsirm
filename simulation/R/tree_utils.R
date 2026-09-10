## ============================================================================
## tree_utils.R — hyperbolic tree-reconstruction primitives and lexicographic
##  two-Gromov decoder (decode_lex), ported verbatim from the rollcall reproduce
##  pipeline (simulation_recovery/decode_tree.R).
##
##  Per draw s and node v, among shallower candidates u (depth gate r_u < r_v):
##    primary : maximise  (u·v)_o = ½(r_u + r_v − D_uv)        [shared-lineage depth, Gromov]
##    tie band: candidates within 1 sigma_v of the top G, where sigma_v = across-draw
##              SD of node v's best (max) Gromov product  [v-centric; more stable than
##              the SD of one fixed (u,v) pair]
##    tiebreak: minimise  e_uv = D_uv − (r_v − r_u) ≥ 0        [off-line gap: ancestor vs sibling]
##  Vote per-draw edges -> PIP_uv (edge posterior inclusion); consensus = each
##  node's max-PIP shallower parent (depth gate => acyclic; polytomies survive).
## ============================================================================

lorentz_gram <- function(W) {
  J <- diag(c(-1, rep(1, ncol(W) - 1L)))
  W %*% J %*% t(W)
}

## Stable polar form, matching src/hyperboloid_utils.h. Directly evaluating
## acosh(-<x,y>_L) subtracts ambient terms of order exp(2r), which loses
## precision for nearby peripheral points.
radial_depth <- function(W) asinh(sqrt(rowSums(W[, -1, drop = FALSE]^2)))
hyperbolic_dist <- function(W) {
  spatial <- W[, -1, drop = FALSE]
  s <- sqrt(rowSums(spatial^2))
  r <- asinh(s)
  u <- spatial / pmax(s, 1e-300)
  chord2 <- matrix(0, nrow(u), nrow(u))
  for (k in seq_len(ncol(u))) chord2 <- chord2 + outer(u[, k], u[, k], "-")^2
  half2 <- pmin(chord2 / 4, 1)
  radial <- sinh(outer(r, r, "-") / 2)
  D <- 2 * asinh(sqrt(radial^2 + outer(s, s, "*") * half2))
  diag(D) <- 0
  D
}
euclidean_dist <- function(W) as.matrix(stats::dist(W, method = "euclidean"))
euclidean_radial_depth <- function(W) sqrt(rowSums(W^2))

## add a virtual root (index 1) at the origin; return distance + radial-depth matrices.
.augment_with_root <- function(W, metric = "hyperbolic") {
  metric <- match.arg(metric, c("hyperbolic", "euclidean"))
  Dp <- ncol(W)
  if (metric == "hyperbolic") {
    Wa <- rbind(c(1, rep(0, Dp - 1)), W)
    list(W_aug = Wa, D = hyperbolic_dist(Wa), r = radial_depth(Wa))
  } else {
    Wa <- rbind(rep(0, Dp), W)
    list(W_aug = Wa, D = euclidean_dist(Wa), r = euclidean_radial_depth(Wa))
  }
}

## decode_core: the parent-rule + consensus-vote engine, factored out of decode_lex so ANY per-draw
## source of (radius, pairwise distance) -- not just the fitted hyperbolic/Euclidean geometry -- can be
## decoded with the identical decision rule. Used by decode_lex itself and by lib_baseline_trees.R's
## NJ/MST/hclust decoders, so a method comparison varies ONLY the upstream distance/tree reconstruction,
## never the parent-assignment logic.
## Dl, rl: length-n lists of per-draw (D, r), each already root-augmented (index 1 = root, r[1] = 0).
decode_core <- function(Dl, rl, k = 1, gap_floor = 0, root = c("virtual", "observed")) {
  root <- match.arg(root)
  n <- length(Dl)
  N <- length(rl[[1]])
  rs <- numeric(N)
  Gstar <- matrix(NA_real_, n, N) # Gstar[s, v] = per-draw max Gromov over shallower
  for (i in seq_len(n)) { # pass 1: per-draw best-Gromov for v
    r <- rl[[i]]
    rs <- rs + r
    for (v in 2:N) {
      cand <- which(r < r[v] - 1e-9)
      cand <- cand[cand != v]
      if (length(cand)) Gstar[i, v] <- max(0.5 * (r[cand] + r[v] - Dl[[i]][cand, v]))
    }
  }
  rmean <- rs / n
  gfloor <- gap_floor * stats::median(rmean[rmean > 1e-6]) # scale-adaptive band floor (Gromov units)
  ## v-centric tie tolerance: across-draw SD of each node's best (max) Gromov product
  ## -- more stable than the SD of one fixed (u,v) pair (matches the G*_v / sigma_v def).
  sigv <- apply(Gstar, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1) sd(x) else 0
  })
  Pm <- matrix(0, N, N)
  for (i in seq_len(n)) { # pass 2: per-draw lex decode, vote
    D <- Dl[[i]]
    r <- rl[[i]]
    for (v in 2:N) {
      cand <- which(r < r[v] - 1e-9)
      cand <- cand[cand != v]
      if (!length(cand)) {
        Pm[1L, v] <- Pm[1L, v] + 1
        next
      }
      G <- 0.5 * (r[cand] + r[v] - D[cand, v])
      e <- D[cand, v] - (r[v] - r[cand])
      top <- which.max(G)
      band <- G[top] - (k * sigv[v] + gfloor)
      tied <- which(G >= band)
      Pm[cand[tied[which.min(e[tied])]], v] <- Pm[cand[tied[which.min(e[tied])]], v] + 1
    }
  }
  Pm <- Pm / n # PIP (edge posterior inclusion)
  par <- rep(NA_integer_, N) # consensus: max-PIP shallower parent
  for (v in 2:N) {
    adm <- which(rmean < rmean[v] - 1e-9)
    if (!length(adm) || all(Pm[adm, v] == 0)) par[v] <- 1L else par[v] <- adm[which.max(Pm[adm, v])]
  }
  if (root == "observed") {
    rk <- which(par == 1L)
    if (length(rk) > 1L) {
      rr <- rk[which.min(rmean[rk])]
      par[rk[rk != rr]] <- rr
    }
  }
  attr(par, "pip") <- Pm
  attr(par, "rmean") <- rmean
  par
}

## gap_floor: widen the Gromov tie-band by a FIXED floor (fraction of median radius) so the off-line-gap
## tiebreak fires even when across-draw variance sigma_v ~ 0. gap_floor=0 -> original sigma_v-gated behaviour.
decode_lex <- function(W, metric = "hyperbolic", sub = 180, k = 1,
                       root = c("virtual", "observed"), gap_floor = 0) {
  root <- match.arg(root)
  use <- round(seq(1, dim(W)[1], length.out = min(sub, dim(W)[1])))
  n <- length(use)
  Dl <- vector("list", n)
  rl <- vector("list", n)
  for (i in seq_len(n)) {
    aug <- .augment_with_root(W[use[i], , ], metric = metric)
    Dl[[i]] <- aug$D
    rl[[i]] <- aug$r
  }
  decode_core(Dl, rl, k = k, gap_floor = gap_floor, root = root)
}
