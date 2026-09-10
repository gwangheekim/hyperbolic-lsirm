# Appendix C: fit the additional HPO subnetworks and rebuild every artefact.
# One command does fit + analysis. SUPP_SKIP_FIT=1 reuses the stored chains;
# APPC_RESUME=1 skips subnetworks whose analysis is already current.

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
setwd(normalizePath(file.path(dirname(script_file), "..")))
source("R/util.R")
suppressPackageStartupMessages(library(posterior))

reported  <- "result/subnetworks/reported"
cases_dir <- file.path(reported, "cases")
dir.create(cases_dir, recursive = TRUE, showWarnings = FALSE)

top_k        <- 10L
support_gate <- 0.70

# One prespecified record per subnetwork, matching the Table C.2-C.5 headings.
case_ids <- c(
  cardiomyopathy_f3             = "OMIM:608758",   # familial hypertrophic cardiomyopathy 10
  leukodystrophy_f3             = "OMIM:601539",   # peroxisome biogenesis disorder 1B
  ciliopathy_f3                 = "OMIM:618161",   # Joubert syndrome 35
  syndromic_craniosynostosis_f3 = "Orphanet:1515" # cranioectodermal dysplasia
)


# 1. Fit -------------------------------------------------------------------

if (!nzchar(Sys.getenv("SUPP_SKIP_FIT"))) {
  suppressPackageStartupMessages(library(Rcpp))
  sourceCpp("cpp/lsirmhyper.cpp", rebuild = TRUE, verbose = FALSE)

  niter <- as.integer(Sys.getenv("NITER", "100000"))
  nburn <- as.integer(Sys.getenv("NBURN", "10000"))
  nthin <- as.integer(Sys.getenv("NTHIN", "10"))
  seeds <- c(20270101L, 20270201L, 20270301L)

  for (tag in setdiff(scope_tags, "dementia_f3")) {
    data_file <- scope_meta_path(tag)
    X <- readRDS(data_file)$X
    storage.mode(X) <- "double"
    dir.create(scope_fit_dir(tag), recursive = TRUE, showWarnings = FALSE)

    for (chain in 1:3) {
      message(sprintf("Fitting %s chain %d", tag, chain))
      set.seed(seeds[chain])
      fit <- do.call(lsirm1pl_cpp, c(list(data = X, niter = niter, nburn = nburn, nthin = nthin), h3_sampler_args))
      fit$pilot_scope <- tag
      fit$pilot_chain <- chain
      fit$seed        <- seeds[chain]
      fit$data_md5    <- unname(tools::md5sum(data_file))
      saveRDS(fit, file.path(scope_fit_dir(tag), sprintf("H3_chain%d.rds", chain)), compress = FALSE)
    }
  }
}


# 2. Per-subnetwork diagnostics and rank stability -----------------------

scalar_diagnostics <- function(fits) {
  quantities <- c("log_posterior", "log_likelihood", "sigma_theta", "sigma_z", "sigma_w")
  do.call(rbind, lapply(quantities, function(nm) {
    draws <- as_draws_matrix(do.call(cbind, lapply(fits, `[[`, nm)))
    chain_means <- vapply(fits, function(f) mean(f[[nm]]), numeric(1))
    data.frame(
      quantity  = nm,
      mean      = mean(unlist(lapply(fits, `[[`, nm))),
      chain_min = min(chain_means),
      chain_max = max(chain_means),
      rhat      = rhat(draws),
      ess_bulk  = ess_bulk(draws),
      ess_tail  = ess_tail(draws),
      method    = "H3"
    )
  }))
}

retention_summary <- function(values, ordering) {
  per_disease <- rowMeans(values)
  data.frame(
    ordering         = ordering,
    median           = median(per_disease),
    p10              = unname(quantile(per_disease, .1)),
    p25              = unname(quantile(per_disease, .25)),
    mean             = mean(per_disease),
    perfect_fraction = mean(per_disease == 1),
    minimum          = min(per_disease)
  )
}

resume <- nzchar(Sys.getenv("APPC_RESUME"))
for (tag in scope_tags) {
  chains  <- file.path(scope_fit_dir(tag), sprintf("H3_chain%d.rds", 1:3))
  outputs <- c(
    file.path(reported, sprintf("%s_scalar_diagnostics.csv", tag)),
    file.path(reported, sprintf("%s_rank_stability.csv", tag)),
    scope_rank_path(tag)
  )
  if (resume && all(file.exists(outputs)) && min(file.mtime(outputs)) > max(file.mtime(chains))) {
    message(sprintf("  %-32s current, skipped", tag))
    next
  }

  fits  <- lapply(chains, read_rds_retry)
  ranks <- rank_scope_fits(tag, top_k = top_k)

  write.csv(scalar_diagnostics(fits),
            file.path(reported, sprintf("%s_scalar_diagnostics.csv", tag)), row.names = FALSE)
  write.csv(rbind(retention_summary(ranks$distance_retention, "distance"),
                  retention_summary(ranks$gromov_retention,   "gromov")),
            file.path(reported, sprintf("%s_rank_stability.csv", tag)), row.names = FALSE)
  message(sprintf("  %-32s rhat=%.4f  retention=%.3f/%.3f",
                  tag, rhat(as_draws_matrix(do.call(cbind, lapply(fits, `[[`, "log_posterior")))),
                  ranks$summary$median[1], ranks$summary$median[2]))
  rm(fits); gc()
}


# 3. Annotation-internal case tables -----------------------------------

case_table <- function(tag, disease_id) {
  meta <- read_rds_retry(scope_meta_path(tag))
  X <- meta$X
  storage.mode(X) <- "double"

  keys       <- apply(X, 2L, paste0, collapse = "")
  profiles   <- unname(split(seq_len(ncol(X)), keys))
  profile_of <- integer(ncol(X))
  for (g in seq_along(profiles)) profile_of[profiles[[g]]] <- g

  record   <- match(disease_id, rownames(X))
  groups   <- unique(profile_of[which(X[record, ] == 1L)])
  columns  <- unlist(profiles[groups], use.names = FALSE)
  col_group <- rep(seq_along(groups), lengths(profiles[groups]))
  carriers <- vapply(groups, function(g) as.integer(sum(X[, profiles[[g]][1L]])), integer(1))
  labels   <- vapply(groups, function(g) paste(unname(meta$lab[profiles[[g]]]), collapse = "; "), character(1))

  per_chain <- lapply(1:3, function(chain) {
    fit <- read_rds_retry(file.path(scope_fit_dir(tag), sprintf("H3_chain%d.rds", chain)))
    n <- nrow(fit$theta)
    d <- dim(fit$z)[3]
    p <- length(columns)
    z <- fit$z[, record, , drop = FALSE]; dim(z) <- c(n, d)
    w <- fit$w[, columns, , drop = FALSE]; dim(w) <- c(n, p, d)

    inner <- matrix(w[, , 1], n, p) * z[, 1]
    spatial_w <- matrix(0, n, p)
    for (k in 2:d) {
      inner <- inner - matrix(w[, , k], n, p) * z[, k]
      spatial_w <- spatial_w + matrix(w[, , k], n, p)^2
    }
    dist <- acosh(pmax(inner, 1))
    grom <- asinh(sqrt(rowSums(z[, -1, drop = FALSE]^2))) + asinh(sqrt(spatial_w)) - dist

    collapse <- function(m) {
      out <- vapply(seq_along(groups), function(g) rowMeans(m[, col_group == g, drop = FALSE]), numeric(n))
      if (length(groups) == 1L) dim(out) <- c(n, 1L)
      out
    }
    list(distance = collapse(dist), gromov = collapse(grom))
  })

  cross_chain_retention <- function(ordering) {
    tops <- lapply(per_chain, function(ch)
      order(colMeans(ch[[ordering]]), decreasing = ordering == "gromov")[seq_len(min(top_k, length(groups)))])
    mean(combn(3, 2, function(p) length(intersect(tops[[p[1]]], tops[[p[2]]])) / top_k))
  }

  ordering_rows <- function(ordering) {
    means  <- lapply(per_chain, function(ch) colMeans(ch[[ordering]]))
    tops   <- lapply(means, function(v) order(v, decreasing = ordering == "gromov")[seq_len(min(top_k, length(v)))])
    common <- Reduce(intersect, tops)
    rows <- do.call(rbind, lapply(common, function(g) {
      chain_ranks <- vapply(means, function(v)
        rank(if (ordering == "distance") v else -v, ties.method = "min")[g], numeric(1))
      support <- min(vapply(per_chain, function(ch) {
        v <- ch[[ordering]]
        mean(1L + rowSums(if (ordering == "distance") v < v[, g] else v > v[, g]) <= top_k)
      }, numeric(1)))
      data.frame(disease = unname(meta$dis_lab[disease_id]), ordering = ordering,
                 profile = labels[g], carriers = carriers[g],
                 chain_ranks = paste(chain_ranks, collapse = "/"),
                 mean_rank = mean(chain_ranks), support = support)
    }))
    rows <- rows[rows$support >= support_gate, ]
    head(rows[order(rows$mean_rank), ], 3L)
  }

  tbl <- rbind(ordering_rows("distance"), ordering_rows("gromov"))
  rownames(tbl) <- NULL
  attr(tbl, "retention")  <- c(distance = cross_chain_retention("distance"),
                               gromov   = cross_chain_retention("gromov"))
  attr(tbl, "annotations") <- sum(X[record, ])
  attr(tbl, "profiles")    <- length(groups)
  tbl
}

case_retention <- list()
case_distance  <- list()
case_gromov    <- list()
for (tag in names(case_ids)) {
  tbl <- case_table(tag, case_ids[[tag]])
  ret <- attr(tbl, "retention")
  write.csv(tbl, file.path(cases_dir, sprintf("%s.csv", tag)), row.names = FALSE)
  case_retention[[tag]] <- data.frame(
    tag = tag, disease_id = case_ids[[tag]],
    annotations = attr(tbl, "annotations"), profiles = attr(tbl, "profiles"),
    retention_distance = ret["distance"], retention_gromov = ret["gromov"]
  )
  case_distance[[tag]] <- tbl$carriers[tbl$ordering == "distance"]
  case_gromov[[tag]]   <- tbl$carriers[tbl$ordering == "gromov"]
  message(sprintf("  %-30s %s  annot/prof=%d/%d  retention=%.4f/%.4f",
                  tag, case_ids[[tag]], attr(tbl, "annotations"), attr(tbl, "profiles"),
                  ret["distance"], ret["gromov"]))
}
write.csv(do.call(rbind, case_retention),
          file.path(reported, "appendix_c_case_retention.csv"), row.names = FALSE)
cat(sprintf("\nCase-table median records: distance %.1f, Gromov %.1f\n",
            median(unlist(case_distance)), median(unlist(case_gromov))))


# 4. Summary table ---------------------------------------------------

summary <- appendix_c_summary_table()
write.csv(summary, file.path(reported, "appendix_c_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
message("Appendix C complete.")
