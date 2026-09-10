## lib_sim.R — tree DGP + fitting + per-level recovery, self-contained server package version.
## Settings live in reproduce/config.R (CFG); this file only holds the functions.
## Reuses the engine setup and tree utilities.
PRES <- Sys.getenv("LSIRM_PRES")
if (!nzchar(PRES)) stop("LSIRM_PRES is not set. Run via run_fits.sh / run_fits_parallel.sh, or `export LSIRM_PRES=/path/to/sim_pipeline_server` first.")
suppressMessages({
  source(file.path(PRES, "R/setup.R"))
  source(file.path(PRES, "R/tree_utils.R"))
})
source(file.path(PRES, "reproduce/config.R"))

## Convert a parent vector (node 1 is the root) into the common tree object used
## by every DGP and decoder. Parents must precede children, which also makes the
## depth calculation deterministic and cheap.
.tree_from_parent <- function(par, branch = NULL, shape = "custom") {
  par <- as.integer(par)
  if (!length(par) || par[1] != 0L ||
    (length(par) > 1L && any(par[-1L] < 1L | par[-1L] >= seq.int(2L, length(par))))) {
    stop("invalid parent vector: root must be first and every parent must precede its child")
  }
  dep <- integer(length(par))
  if (length(par) > 1L) for (v in 2:length(par)) dep[v] <- dep[par[v]] + 1L
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
  items <- which(dep > 0)
  out <- list(
    par = par, dep = dep, items = items, depi = dep[items],
    G = outer(seq_len(NN), seq_len(NN), Vectorize(gd)),
    tp = vapply(items, function(v) {
      p <- par[v]
      if (p == 1L) 0L else match(p, items)
    }, integer(1)),
    shape = shape
  )
  if (!is.null(branch)) {
    if (length(branch) != NN) stop("branch labels must have one entry per node")
    out$branch <- as.integer(branch)
  }
  out
}

## Balanced branching tree; items are all non-root nodes.
build_tree <- function(DEPTH = CFG$DEPTH, BF = CFG$BF) {
  par <- 0L
  bld <- function(node, lv) {
    if (lv >= DEPTH) {
      return()
    }
    for (k in seq_len(BF)) {
      par <<- c(par, node)
      bld(length(par), lv + 1L)
    }
  }
  bld(1L, 0L)
  .tree_from_parent(par, shape = sprintf("balanced-depth%d-BF%d", DEPTH, BF))
}

## Shallow, deliberately unbalanced tree with the SAME 30 non-root nodes as the
## baseline. Counting the root as L1, two L2 branches carry 5 versus 23 L3
## leaves. The large side is 23 (rather than the illustrative 20) because
## 2 + 5 + 23 = 30 exactly.
build_shallow_unbalanced <- function(nitems = 30L, small_leaves = 5L) {
  nitems <- as.integer(nitems)
  small_leaves <- as.integer(small_leaves)
  large_leaves <- nitems - 2L - small_leaves
  if (nitems < 4L || small_leaves < 1L || large_leaves < 1L) {
    stop("shallow unbalanced tree needs two branches and at least one leaf per branch")
  }
  par <- c(
    0L, 1L, 1L,
    rep(2L, small_leaves), rep(3L, large_leaves)
  )
  branch <- c(
    0L, 1L, 2L,
    rep(1L, small_leaves), rep(2L, large_leaves)
  )
  out <- .tree_from_parent(
    par, branch,
    sprintf("shallow-unbalanced-%d-%d", small_leaves, large_leaves)
  )
  out$leaf_counts <- c(small = small_leaves, large = large_leaves)
  out
}

## Unequal-depth tree with a fixed node budget. Four root branches have maximum
## depths 2, 3, 4, and 6 by default. A spine first guarantees each requested
## depth; remaining nodes are added breadth-first inside the branches without
## exceeding their caps, so local branching also varies while M stays fixed.
build_depth_unbalanced <- function(nitems = 30L, branch_depths = c(2L, 3L, 4L, 6L)) {
  nitems <- as.integer(nitems)
  branch_depths <- as.integer(branch_depths)
  if (!length(branch_depths) || any(branch_depths < 1L)) {
    stop("branch depths must be positive integers")
  }
  minimum_items <- sum(branch_depths)
  if (nitems < minimum_items) {
    stop(sprintf("nitems must be at least %d to realize all branch depths", minimum_items))
  }

  par <- 0L
  dep <- 0L
  branch <- 0L
  roots <- integer(length(branch_depths))
  branch_nodes <- vector("list", length(branch_depths))
  for (b in seq_along(branch_depths)) {
    par <- c(par, 1L)
    dep <- c(dep, 1L)
    branch <- c(branch, b)
    roots[b] <- length(par)
    branch_nodes[[b]] <- roots[b]
  }
  for (b in seq_along(branch_depths)) {
    parent <- roots[b]
    if (branch_depths[b] >= 2L) {
      for (d in 2:branch_depths[b]) {
        par <- c(par, parent)
        dep <- c(dep, d)
        branch <- c(branch, b)
        parent <- length(par)
        branch_nodes[[b]] <- c(branch_nodes[[b]], parent)
      }
    }
  }
  child_count <- integer(length(par))
  if (length(par) > 1L) for (v in 2:length(par)) child_count[par[v]] <- child_count[par[v]] + 1L
  next_branch <- 1L
  while (length(par) - 1L < nitems) {
    b <- next_branch
    cand <- branch_nodes[[b]][dep[branch_nodes[[b]]] < branch_depths[b]]
    if (!length(cand)) stop("internal error: no eligible parent below branch cap")
    ord <- order(child_count[cand], dep[cand], cand)
    parent <- cand[ord[1L]]
    par <- c(par, parent)
    dep <- c(dep, dep[parent] + 1L)
    branch <- c(branch, b)
    child_count <- c(child_count, 0L)
    child_count[parent] <- child_count[parent] + 1L
    branch_nodes[[b]] <- c(branch_nodes[[b]], length(par))
    next_branch <- if (b == length(branch_depths)) 1L else b + 1L
  }
  out <- .tree_from_parent(
    par, branch,
    paste0("depth-unbalanced-", paste(branch_depths, collapse = "-"))
  )
  out$branch_depths <- branch_depths
  out
}

## binary responses from person-item graph distance; theta-tiling (spread) + gamma.
## placement="perLevel" (default, unchanged paper behaviour): persons placed UNIFORMLY OVER DEPTH (equal
##   count per level) so shallow internal nodes are well identified -- but since node count grows with
##   depth (2,4,8,16 for depth4/BF2), this gives deep leaves far FEWER local anchors per node (~23/node at
##   depth4 vs ~187/node at depth1 for N=1500).
## placement="perItem": persons placed UNIFORMLY OVER ITEMS (equal count per node), so deep leaves get
##   proportionally more anchors, matching their greater fan-out. See
##   presentation/reproduce/RADIAL_COMPRESSION_DIAG.md: this was tested as a possible fix for the fitted
##   radial-compression bias and found NOT to move it (ratio-by-depth barely changes) -- kept here as a
##   config option for the separate perItem-placement production comparison
##   (presentation/reproduce/fit_sim_perItem.R), not because it fixes that specific bias.
gen_binary <- function(t, N = CFG$N, seed, gamma = CFG$GAMMA, spread = CFG$SPREAD,
                       placement = c("perLevel", "perItem")) {
  placement <- match.arg(placement)
  items <- t$items
  M <- length(items)
  di <- t$depi
  dl <- sort(unique(di))
  set.seed(seed)
  pn <- if (placement == "perLevel") {
    unlist(lapply(dl, function(d) {
      it <- items[di == d]
      sample(it, ceiling(N / length(dl)), replace = TRUE)
    }))[1:N]
  } else {
    sample(items, N, replace = TRUE)
  }
  Dpi <- t$G[pn, items]
  th <- rnorm(length(pn), 0, spread)
  mu <- mean(gamma * Dpi)
  Y <- matrix(rbinom(length(pn) * M, 1, 1 / (1 + exp(-(th - gamma * Dpi + mu)))), length(pn), M)
  storage.mode(Y) <- "double"
  Y
}

## General paper fit. The current engine fixes gamma=1, estimates sigma_z and
## sigma_w under half-Normal hyperpriors, and uses wrapped-normal positions.
## Legacy fixscale, wrapped-t and fixed-position-scale branches remain only in
## and are intentionally not exposed in this submission workflow.
fit_embed <- function(Y, hyp, mdim = 2L, niter = CFG$NITER, param = CFG$PARAM,
                      dof = 0, sd_fix = NULL, estimate_beta_sd = FALSE) {
  if (!identical(param, "fixgamma")) stop("the paper engine supports param='fixgamma' only")
  if (length(dof) != 1L || !is.finite(dof) || dof != 0 || !is.null(sd_fix)) {
    stop("wrapped-t and fixed-position-scale modes are not included in this submission workflow")
  }
  if (length(niter) != 1L || !is.finite(niter) || niter != as.integer(niter) || niter < 60L) {
    stop("niter must be an integer >= 60")
  }
  niter <- as.integer(niter)
  nburn <- min(CFG$NBURN, niter %/% 2L)
  first_saved <- ceiling(nburn / CFG$NTHIN) * CFG$NTHIN
  if (first_saved >= niter) stop("the requested MCMC schedule retains no draws")
  ndim <- if (hyp) as.integer(mdim + 1L) else as.integer(mdim)
  engine_args <- list(
    data = Y, ndim = ndim, niter = niter,
    nburn = nburn, nthin = CFG$NTHIN,
    jump_beta = .3, jump_theta = .5, jump_z = .15, jump_w = .15,
    pr_sd_beta = CFG$PR_SD_BETA, pr_a_beta = .001, pr_b_beta = .001,
    pr_a_theta = .001, pr_b_theta = .001, pr_sd_theta = 1,
    pr_sd_z = CFG$PR_SD_Z, pr_sd_w = CFG$PR_SD_Z,
    hyperbolic = hyp, estimate_beta_sd = estimate_beta_sd,
    adaptation = TRUE, adapt_interval = 25L, target_accept_1d = .44, adapt_exponent = .6,
    stable_intrinsic_proposal = CFG$STABLE_INTRINSIC_PROPOSAL,
    max_hyperbolic_radius = CFG$MAX_HYPERBOLIC_RADIUS
  )
  engine_formals <- names(formals(lsirm1pl_cpp))
  radial_control <- if (hyp) {
    list(enabled = CFG$RADIAL_ASIS, jump = CFG$RADIAL_ASIS_JUMP)
  } else {
    NULL
  }
  if ("radial_block_control" %in% engine_formals) {
    engine_args$radial_block_control <- radial_control
    engine_args$polar_split <- hyp && isTRUE(CFG$POLAR_SPLIT)
  } else if ("radial_asis_control" %in% engine_formals) {
    if (hyp && isTRUE(CFG$POLAR_SPLIT)) {
      stop("POLAR_SPLIT=true requires the current src/lsirmhyper.cpp engine")
    }
    engine_args$radial_asis_control <- radial_control
  } else {
    stop("unsupported sampler interface: no radial control argument found")
  }
  invisible(utils::capture.output(f <- do.call(lsirm1pl_cpp, engine_args),
    file = nullfile()
  ))

  ## Normalize renamed diagnostics so the simulation readers remain compatible
  ## with the supplied C++ engine.
  if (is.null(f$use_radial_asis)) f$use_radial_asis <- isTRUE(f$use_radial_block)
  if (is.null(f$use_radial_block)) f$use_radial_block <- isTRUE(f$use_radial_asis)
  if (is.null(f$polar_split)) f$polar_split <- FALSE
  if (is.null(f$use_gromov)) f$use_gromov <- FALSE
  f$geom <- if (hyp) "hyp" else "euc"
  f$mdim <- mdim
  f$param <- param
  f$dof <- 0
  f$sd_fix <- CFG$PR_SD_Z
  f$estimate_beta_sd <- estimate_beta_sd
  f$run_settings <- list(
    stable_intrinsic_proposal = CFG$STABLE_INTRINSIC_PROPOSAL,
    max_hyperbolic_radius = CFG$MAX_HYPERBOLIC_RADIUS,
    radial_asis = CFG$RADIAL_ASIS,
    radial_asis_jump = CFG$RADIAL_ASIS_JUMP,
    polar_split = hyp && isTRUE(CFG$POLAR_SPLIT),
    niter = niter, nburn = nburn, nthin = CFG$NTHIN,
    expected_retained = length(seq.int(first_saved, niter - 1L,
      by = CFG$NTHIN
    ))
  )
  f
}
## Back-compatible H^2 / R^2 helper returning only item draws W.
fit_geom <- function(Y, hyp, niter = 8000L) fit_embed(Y, hyp, mdim = 2L, niter = niter, param = "fixgamma")$w
## item-draws accessor: full-fit object -> $w; a bare W array passes through (readers stay format-robust).
getW <- function(x) if (is.list(x) && !is.null(x$w)) x$w else x

## multiplicative distortion of fitted distances vs the true distances (optimal global scale c*; lower=better).
distortion <- function(Dhat, Dtrue) {
  ut <- upper.tri(Dtrue)
  d <- Dhat[ut]
  g <- Dtrue[ut]
  o <- g > 0 & is.finite(g) & is.finite(d) & d >= 0
  den <- sum(d[o]^2)
  if (!any(o) || den <= 0) stop("distortion is undefined for these distances")
  c <- sum(d[o] * g[o]) / den
  mean(abs(c * d[o] - g[o]) / g[o])
}

## Paired per-seed H2/R2 baseline scores produced once by table_distortion.R.
## Reported sensitivity scripts reuse these values so every table refers to the
## same baseline fit and no identical 50,000-iteration chain is recomputed.
read_baseline_pair <- function() {
  p <- file.path(PRES, "reproduce", "results", "results_distortion_raw.csv")
  if (!file.exists(p)) stop("baseline scores are missing; run fit_sim.R and table_distortion.R first")
  x <- read.csv(p, stringsAsFactors = FALSE)
  h <- x[x$embed == "H2", c("seed", "distortion")]
  r <- x[x$embed == "R2", c("seed", "distortion")]
  names(h)[2] <- "H2"
  names(r)[2] <- "R2"
  z <- merge(h, r, by = "seed", all = TRUE, sort = TRUE)
  if (nrow(z) != length(CFG$SEEDS) || anyNA(z) || !identical(as.integer(z$seed), as.integer(CFG$SEEDS))) {
    stop("baseline H2/R2 scores do not match CFG$SEEDS")
  }
  z$gap <- z$R2 - z$H2
  z
}

## Exact normalized four-point hyperbolicity for the small (M=30) truths used
## by the negative control.
delta_over_diam_exact <- function(D) {
  if (nrow(D) < 4L || max(D) <= 0) {
    return(0)
  }
  idx <- utils::combn(nrow(D), 4L)
  best <- 0
  for (q in seq_len(ncol(idx))) {
    d <- D[idx[, q], idx[, q]]
    s <- sort(c(d[1, 2] + d[3, 4], d[1, 3] + d[2, 4], d[1, 4] + d[2, 3]), decreasing = TRUE)
    best <- max(best, (s[1] - s[2]) / 2)
  }
  best / max(D)
}

## score a decoded parent vector (as returned by decode_lex / decode_baseline: index1=virtual root, item j
## at decode index j+1) against the true parent structure t$tp. Shared by eval_recovery and any other
## decoder (e.g. lib_baseline_trees.R's decode_baseline) so all methods are scored identically.
score_pm <- function(pm, t) {
  M <- length(t$items)
  prd <- vapply(1:M, function(v) {
    p <- pm[v + 1L]
    if (is.na(p) || p == 1L) 0L else p - 1L
  }, integer(1))
  ok <- t$tp > 0
  per <- sapply(sort(unique(t$depi[ok])), function(L) {
    i <- which(ok & t$depi == L)
    mean(prd[i] == t$tp[i])
  })
  names(per) <- paste0("L", sort(unique(t$depi[ok])))
  list(prd = prd, overall = mean(prd[ok] == t$tp[ok]), per_level = per)
}

## decode + per-level parent-recovery accuracy. tp aug: index1=virtual root; item j -> decode index j+1.
## NOTE: decode_lex is NOT fork-safe -> call eval_recovery from a SERIAL loop (lapply), never mclapply.
eval_recovery <- function(W, t, metric, gap_floor = 0) {
  pm <- decode_lex(W, metric, sub = CFG$SUB, k = 1, root = "virtual", gap_floor = gap_floor)
  score_pm(pm, t)
}


## mclapply returns a worker's error as a try-error OBJECT instead of raising it, so a plain
## do.call(cbind, mclapply(...)) builds a CHARACTER matrix and the failure only surfaces much later
## as "rowMeans: 'x' must be numeric" -- with the worker's real message already lost.  Every seed
## loop in this pipeline goes through mc_cols() so the original message is printed where it happens
## and the nominal replicate count cannot silently change.
mc_cols <- function(seeds, fn, cores = sim_cores(length(seeds)), lab = "", nrow_expect = NULL) {
  L <- parallel::mclapply(seeds, fn, mc.cores = cores, mc.preschedule = FALSE)
  ok <- vapply(L, function(x) {
    is.numeric(x) && !anyNA(x) &&
      (is.null(nrow_expect) || length(x) == nrow_expect)
  }, logical(1))
  if (any(!ok)) {
    for (i in which(!ok)) {
      message(sprintf(
        "  [%s] seed %s failed: %s", lab, seeds[i],
        paste(trimws(as.character(L[[i]])), collapse = " ")
      ))
    }
    stop(sprintf(
      "%d of %d seeds failed%s -- canonical summaries require all seeds",
      sum(!ok), length(seeds), if (nzchar(lab)) paste0(" for '", lab, "'") else ""
    ))
  }
  m <- do.call(cbind, L[ok])
  attr(m, "nseed_ok") <- sum(ok)
  attr(m, "nseed") <- length(seeds)
  m
}

## Strict list-valued parallel map. A failed fork must fail the Slurm job rather
## than silently reduce the replicate count under a nominal 100-seed label.
mc_list <- function(X, fn, cores = sim_cores(length(X)), lab = "") {
  L <- parallel::mclapply(X, fn, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(L, inherits, logical(1), "try-error") | vapply(L, is.null, logical(1))
  if (any(bad)) {
    for (i in which(bad)) {
      message(sprintf(
        "  [%s] task %s failed: %s", lab, X[i],
        paste(trimws(as.character(L[[i]])), collapse = " ")
      ))
    }
    stop(sprintf(
      "%d of %d tasks failed%s", sum(bad), length(X),
      if (nzchar(lab)) paste0(" for '", lab, "'") else ""
    ))
  }
  L
}
