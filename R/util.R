# Shared configuration and helpers for R/main.R and R/supplementary.R.

# Appendix C subnetworks. dementia_f3 is the Section 6 subnetwork; its row is
# built from the R/main.R fit.
scope_tags <- c(
  "movement_f3", "psychiatric_f3", "dementia_f3", "motor_neuron_f3",
  "intellectual_disability_f3", "cerebellar_f3", "basal_ganglia_f3",
  "ciliopathy_f3", "syndromic_craniosynostosis_f3", "leukodystrophy_f3",
  "hereditary_spastic_paraplegia_f3", "cardiomyopathy_f3",
  "lysosomal_storage_f3", "neurodegenerative_f3", "neurodevelopmental_f3"
)

scope_names <- c(
  movement_f3                      = "Movement disorders",
  psychiatric_f3                   = "Psychiatric disorders",
  dementia_f3                      = "Dementia",
  motor_neuron_f3                  = "Motor neuron disorders",
  intellectual_disability_f3       = "Intellectual disability",
  cerebellar_f3                    = "Cerebellar disorders",
  basal_ganglia_f3                 = "Basal ganglia disorders",
  ciliopathy_f3                    = "Ciliopathy",
  syndromic_craniosynostosis_f3    = "Syndromic craniosynostosis",
  leukodystrophy_f3                = "Leukodystrophy",
  hereditary_spastic_paraplegia_f3 = "Hereditary spastic paraplegia",
  cardiomyopathy_f3                = "Cardiomyopathy",
  lysosomal_storage_f3             = "Lysosomal storage diseases",
  neurodegenerative_f3             = "Neurodegenerative disease",
  neurodevelopmental_f3            = "Neurodevelopmental disorder"
)

# Sampler arguments shared by every H3 fit. ndim = 4 is H3.
h3_sampler_args <- list(
  ndim = 4L,
  jump_beta = .3, jump_theta = .5, jump_z = .15, jump_w = .15,
  pr_sd_beta = 1, pr_a_beta = .001, pr_b_beta = .001,
  pr_a_theta = .001, pr_b_theta = .001, pr_sd_theta = 1,
  pr_sd_z = 2.5, pr_sd_w = 2.5,
  hyperbolic = TRUE, estimate_beta_sd = FALSE,
  adaptation = TRUE, adapt_interval = 25L,
  target_accept_1d = .44, adapt_exponent = .6,
  radial_block_control = list(enabled = TRUE, jump = .15),
  stable_intrinsic_proposal = TRUE,
  max_hyperbolic_radius = 13,
  polar_split = TRUE
)

scope_meta_path <- function(tag) {
  if (tag == "dementia_f3") {
    "data/analysis_data/dementia/hpo_dementia_f3_meta.rds"
  } else {
    sprintf("data/analysis_data/subnetworks/hpo_%s_meta.rds", tag)
  }
}

scope_fit_dir <- function(tag) {
  if (tag == "dementia_f3") "result/model_fits/dementia"
  else file.path("result/model_fits/subnetworks", tag)
}

scope_rank_path <- function(tag) file.path("result/subnetworks/ranking_summaries", paste0(tag, ".rds"))

# The fits sit on a network volume that can drop a connection mid-read; retry a
# few times before giving up.
read_rds_retry <- function(path, attempts = 5L, wait = 10) {
  for (attempt in seq_len(attempts)) {
    value <- tryCatch(readRDS(path), error = function(e) e)
    if (!inherits(value, "error")) return(value)
    if (attempt == attempts) stop("readRDS failed for ", path, ": ", conditionMessage(value))
    message("  readRDS retry ", attempt, " for ", basename(path))
    Sys.sleep(wait)
  }
}

# Posterior-mean pairwise hyperbolic distance and Gromov product for one chain,
# over every retained draw.
mean_scope_scores <- function(fit) {
  distance <- matrix(0, ncol(fit$theta), ncol(fit$beta))
  gromov   <- matrix(0, ncol(fit$theta), ncol(fit$beta))
  for (draw in seq_len(nrow(fit$theta))) {
    z <- fit$z[draw, , ]
    w <- fit$w[draw, , ]
    d <- acosh(pmax(outer(z[, 1], w[, 1]) - z[, -1, drop = FALSE] %*% t(w[, -1, drop = FALSE]), 1))
    rz <- asinh(sqrt(rowSums(z[, -1, drop = FALSE]^2)))
    rw <- asinh(sqrt(rowSums(w[, -1, drop = FALSE]^2)))
    distance <- distance + d
    gromov   <- gromov + (outer(rz, rw, "+") - d) / 2
  }
  n <- nrow(fit$theta)
  list(distance = distance / n, gromov = gromov / n)
}

# Cross-chain leading-ten retention for one subnetwork; also cached to disk.
rank_scope_fits <- function(tag, top_k = 10L) {
  scores <- lapply(
    lapply(1:3, function(chain) read_rds_retry(file.path(scope_fit_dir(tag), sprintf("H3_chain%d.rds", chain)))),
    mean_scope_scores
  )
  leading <- function(field, decreasing) lapply(scores, function(s)
    t(apply(s[[field]], 1, function(x) head(order(x, decreasing = decreasing), top_k))))
  top_distance <- leading("distance", FALSE)
  top_gromov   <- leading("gromov", TRUE)

  overlap <- function(tops) sapply(seq_len(ncol(combn(3, 2))), function(k) {
    pair <- combn(3, 2)[, k]
    vapply(seq_len(nrow(tops[[1]])), function(i)
      length(intersect(tops[[pair[1]]][i, ], tops[[pair[2]]][i, ])) / top_k, numeric(1))
  })
  distance_retention <- overlap(top_distance)
  gromov_retention   <- overlap(top_gromov)

  result <- list(
    summary = data.frame(
      ordering = c("distance", "gromov"),
      median   = c(median(rowMeans(distance_retention)), median(rowMeans(gromov_retention))),
      mean     = c(mean(rowMeans(distance_retention)),   mean(rowMeans(gromov_retention)))
    ),
    distance_retention = distance_retention,
    gromov_retention   = gromov_retention,
    scores = scores, top_distance = top_distance, top_gromov = top_gromov, top_k = top_k
  )
  dir.create(dirname(scope_rank_path(tag)), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, scope_rank_path(tag))
  result
}

# One row per subnetwork: matrix size, split-Rhat for the log posterior and log
# likelihood, and median leading-ten retention.
appendix_c_summary_table <- function() {
  reported <- "result/subnetworks/reported"
  do.call(rbind, lapply(names(scope_names), function(tag) {
    meta   <- read_rds_retry(scope_meta_path(tag))
    rank   <- readRDS(scope_rank_path(tag))
    scalar <- read.csv(file.path(reported, sprintf("%s_scalar_diagnostics.csv", tag)))
    data.frame(
      tag                 = tag,
      subnetwork          = unname(scope_names[tag]),
      diseases            = nrow(meta$X),
      phenotypes          = ncol(meta$X),
      annotations         = sum(meta$X),
      rhat_log_posterior  = scalar$rhat[scalar$quantity == "log_posterior"],
      rhat_log_likelihood = scalar$rhat[scalar$quantity == "log_likelihood"],
      distance_retention  = rank$summary$median[rank$summary$ordering == "distance"],
      gromov_retention    = rank$summary$median[rank$summary$ordering == "gromov"]
    )
  }))
}
