# Section 6: build the dementia matrix, fit three H3 chains, and rebuild every
# Section 6 artefact. One command does everything. S6_SKIP=1 fits only;
# SKIP_MATRIX_BUILD=1 reuses the stored matrix.

suppressPackageStartupMessages(library(Rcpp))
have_igraph <- requireNamespace("igraph", quietly = TRUE)

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
setwd(normalizePath(file.path(dirname(script_file), "..")))
source("R/util.R")

data_file    <- scope_meta_path("dementia_f3")
reported_dir <- Sys.getenv("S6_REPORTED_DIR", "result/dementia/reported")
top_k        <- 10L


# 1. Dementia response matrix ------------------------------------------
# Mondo descendants of the dementia class; direct affirmative PCS/TAS
# phenotypic-abnormality annotations; phenotypes in at least three records.

raw_files <- c("data/raw/mondo/mondo.obo", "data/raw/hpo/phenotype.hpoa", "data/raw/hpo/hp.obo")
if (!nzchar(Sys.getenv("SKIP_MATRIX_BUILD")) && all(file.exists(raw_files))) {
  mondo <- readLines("data/raw/mondo/mondo.obo", warn = FALSE)
  term_start <- grep("^\\[Term\\]$", mondo)
  term_end   <- c(term_start[-1] - 1L, length(mondo))
  ids <- character(length(term_start))
  parents <- xrefs <- vector("list", length(term_start))
  obsolete <- logical(length(term_start))
  for (i in seq_along(term_start)) {
    block <- mondo[term_start[i]:term_end[i]]
    ids[i]      <- sub("^id: ", "", grep("^id: MONDO:", block, value = TRUE)[1])
    parents[[i]] <- sub(" .*", "", sub("^is_a: ", "", grep("^is_a: MONDO:", block, value = TRUE)))
    xrefs[[i]]   <- sub(" .*", "", sub("^xref: ", "", grep("^xref: (OMIM|Orphanet|DECIPHER):", block, value = TRUE)))
    obsolete[i] <- any(block == "is_obsolete: true")
  }
  keep <- !is.na(ids) & !obsolete
  ids <- ids[keep]; parents <- parents[keep]; xrefs <- xrefs[keep]
  names(parents) <- names(xrefs) <- ids

  children <- list()
  for (i in seq_along(ids)) for (parent in parents[[i]]) children[[parent]] <- c(children[[parent]], ids[i])
  descendants <- function(node) {
    found <- character(); frontier <- node
    while (length(frontier)) {
      frontier <- setdiff(unique(unlist(children[frontier], use.names = FALSE)), found)
      found <- c(found, frontier)
    }
    unique(c(node, found))
  }

  scope       <- descendants("MONDO:0001627")
  disease_ids <- unique(unlist(xrefs[intersect(scope, names(xrefs))], use.names = FALSE))
  annot <- read.delim("data/raw/hpo/phenotype.hpoa", comment.char = "#", quote = "", stringsAsFactors = FALSE)
  annot$database_id <- sub("^ORPHA:", "Orphanet:", annot$database_id)
  annot <- annot[annot$database_id %in% disease_ids & annot$aspect == "P" &
                   annot$evidence %in% c("PCS", "TAS") &
                   (is.na(annot$qualifier) | annot$qualifier != "NOT"),
                 c("database_id", "disease_name", "hpo_id")]
  dis_lab <- annot$disease_name[!duplicated(annot$database_id)]
  names(dis_lab) <- annot$database_id[!duplicated(annot$database_id)]

  annot     <- unique(annot[c("database_id", "hpo_id")])
  diseases   <- sort(unique(annot$database_id))
  phenotypes <- sort(unique(annot$hpo_id))
  X <- matrix(0L, length(diseases), length(phenotypes), dimnames = list(diseases, phenotypes))
  X[cbind(match(annot$database_id, diseases), match(annot$hpo_id, phenotypes))] <- 1L
  X <- X[, colSums(X) >= 3L]
  X <- X[rowSums(X) >= 1L, ]

  hpo <- readLines("data/raw/hpo/hp.obo", warn = FALSE)
  hpo_id   <- sub("^id: ", "", grep("^id: HP:", hpo, value = TRUE))
  hpo_name <- sub("^name: ", "", grep("^name: ", hpo, value = TRUE))
  term_labels <- setNames(hpo_name[match(grep("^id: HP:", hpo) + 1L, grep("^name: ", hpo))], hpo_id)

  meta <- list(X = X, id = colnames(X), lab = term_labels[colnames(X)], dis_lab = dis_lab[rownames(X)],
               mondo_release = "2026-08-04", hpo_release = "2026-06-06", evidence_filter = "PCS+TAS only")
  dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(meta, data_file)
  message(sprintf("Wrote %s: %d x %d, %d annotations", data_file, nrow(X), ncol(X), sum(X)))
}

meta <- readRDS(data_file)
X <- meta$X
storage.mode(X) <- "double"
disease_ids <- rownames(X)
Nd <- nrow(X)


# 2. Fit three H3 chains ---------------------------------------------

sourceCpp("cpp/lsirmhyper.cpp", rebuild = TRUE, verbose = FALSE)

niter <- as.integer(Sys.getenv("NITER", "100000"))
nburn <- as.integer(Sys.getenv("NBURN", "10000"))
nthin <- as.integer(Sys.getenv("NTHIN", "10"))
seeds <- c(20270101L, 20270201L, 20270301L)
fit_dir <- Sys.getenv("FIT_OUTPUT_DIR", "result/model_fits/dementia")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

fit_chain <- function(chain) {
  set.seed(seeds[chain])
  fit <- do.call(lsirm1pl_cpp, c(list(data = X, niter = niter, nburn = nburn, nthin = nthin), h3_sampler_args))
  fit$pilot_scope <- "dementia_f3"
  fit$pilot_chain <- chain
  fit$seed        <- seeds[chain]
  fit$data_md5    <- unname(tools::md5sum(data_file))
  fit
}

fits <- parallel::mclapply(1:3, fit_chain, mc.cores = 3, mc.preschedule = FALSE, mc.set.seed = FALSE)
for (chain in 1:3) {
  saveRDS(fits[[chain]], file.path(fit_dir, sprintf("H3_chain%d.rds", chain)), compress = FALSE)
  message("Saved chain ", chain)
}

if (nzchar(Sys.getenv("S6_SKIP"))) quit(save = "no")
dir.create(reported_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Disease geometry ----------------------------------------------
# Pooled posterior-mean pairwise hyperbolic distance and radius over every draw.

geo_distance <- matrix(0, Nd, Nd)
geo_radius   <- numeric(Nd)
draws <- 0L
for (fit in fits) {
  for (t in seq_len(dim(fit$z)[1])) {
    z <- fit$z[t, , ]
    geo_distance <- geo_distance + acosh(pmax(outer(z[, 1], z[, 1]) - z[, -1, drop = FALSE] %*% t(z[, -1, drop = FALSE]), 1))
    geo_radius   <- geo_radius + asinh(sqrt(rowSums(z[, -1, drop = FALSE]^2)))
  }
  draws <- draws + dim(fit$z)[1]
}
geo_distance <- geo_distance / draws
geo_radius   <- geo_radius / draws
dimnames(geo_distance) <- list(disease_ids, disease_ids)
names(geo_radius) <- disease_ids

community <- setNames(rep(NA_integer_, Nd), disease_ids)
if (have_igraph) {
  projection <- X %*% t(X); diag(projection) <- 0
  set.seed(1L)
  community <- as.integer(igraph::membership(igraph::cluster_louvain(
    igraph::graph_from_adjacency_matrix(projection, mode = "undirected", weighted = TRUE))))
  names(community) <- disease_ids
}

saveRDS(list(
  distance = geo_distance, radius = geo_radius,
  gromov = (outer(geo_radius, geo_radius, "+") - geo_distance) / 2,
  draws = draws, community = community,
  labels = setNames(unname(meta$dis_lab[disease_ids]), disease_ids)
), file.path(reported_dir, "dementia_f3_geometry.rds"))
message("  geometry: ", draws, " draws")


# 4. Full-vocabulary rank stability ------------------------------

chain_scores <- function(fit) {
  distance <- matrix(0, ncol(fit$theta), ncol(fit$beta))
  gromov   <- matrix(0, ncol(fit$theta), ncol(fit$beta))
  for (t in seq_len(nrow(fit$theta))) {
    z <- fit$z[t, , ]; w <- fit$w[t, , ]
    d <- acosh(pmax(outer(z[, 1], w[, 1]) - z[, -1, drop = FALSE] %*% t(w[, -1, drop = FALSE]), 1))
    rz <- asinh(sqrt(rowSums(z[, -1, drop = FALSE]^2)))
    rw <- asinh(sqrt(rowSums(w[, -1, drop = FALSE]^2)))
    distance <- distance + d
    gromov   <- gromov + (outer(rz, rw, "+") - d) / 2
  }
  list(distance = distance / nrow(fit$theta), gromov = gromov / nrow(fit$theta))
}

scores <- lapply(fits, chain_scores)
leading <- function(field, decreasing) lapply(scores, function(s)
  t(apply(s[[field]], 1, function(x) head(order(x, decreasing = decreasing), top_k))))
top_distance <- leading("distance", FALSE)
top_gromov   <- leading("gromov", TRUE)

pair_overlap <- function(tops) sapply(1:3, function(k) {
  pair <- combn(3, 2)[, k]
  vapply(seq_len(nrow(tops[[1]])), function(i)
    length(intersect(tops[[pair[1]]][i, ], tops[[pair[2]]][i, ])) / top_k, numeric(1))
})
distance_retention <- pair_overlap(top_distance)
gromov_retention   <- pair_overlap(top_gromov)

retention_row <- function(values, ordering) data.frame(
  ordering = ordering, median = median(rowMeans(values)),
  p10 = unname(quantile(rowMeans(values), .1)), p25 = unname(quantile(rowMeans(values), .25)),
  mean = mean(rowMeans(values)), perfect_fraction = mean(rowMeans(values) == 1),
  minimum = min(rowMeans(values))
)
rank_summary <- rbind(retention_row(distance_retention, "distance"),
                      retention_row(gromov_retention, "gromov"))
saveRDS(list(summary = rank_summary, distance_retention = distance_retention,
             gromov_retention = gromov_retention, top_distance = top_distance,
             top_gromov = top_gromov, scores = scores, top_k = top_k),
        file.path(reported_dir, "rank_stability.rds"))
write.csv(rank_summary, file.path(reported_dir, "rank_stability_summary.csv"), row.names = FALSE)


# 5. Exact profiles and annotation-internal summaries --------------

keys       <- apply(X, 2L, paste0, collapse = "")
profiles   <- unname(split(seq_len(ncol(X)), keys))
profile_of <- integer(ncol(X))
for (g in seq_along(profiles)) profile_of[profiles[[g]]] <- g
row_profiles <- lapply(seq_len(nrow(X)), function(r) unique(profile_of[which(X[r, ] == 1L)]))
carriers <- vapply(profiles, function(cols) as.integer(sum(X[, cols[1L]])), integer(1))
labels   <- vapply(profiles, function(cols) paste(unname(meta$lab[cols]), collapse = "; "), character(1))
eligible <- which(lengths(row_profiles) >= top_k)

set_overlap <- function(sets) mean(combn(3, 2, function(p)
  length(intersect(sets[[p[1]]], sets[[p[2]]])) / top_k))

full_by_chain <- do.call(rbind, lapply(1:3, function(chain) {
  td <- top_distance[[chain]][eligible, ]
  tg <- top_gromov[[chain]][eligible, ]
  data.frame(
    chain = chain,
    median_carriers_distance = median(carriers[profile_of[td]]),
    median_carriers_gromov   = median(carriers[profile_of[tg]]),
    share_le10_distance = mean(carriers[profile_of[td]] <= 10),
    share_le10_gromov   = mean(carriers[profile_of[tg]] <= 10),
    mean_overlap = mean(vapply(seq_len(nrow(td)), function(i)
      length(intersect(td[i, ], tg[i, ])), integer(1)))
  )
}))

disease_rows   <- list()
consensus_rows <- list()
for (record in eligible) {
  ids <- row_profiles[[record]]
  d_scores <- lapply(scores, function(s) vapply(ids, function(g) mean(s$distance[record, profiles[[g]]]), numeric(1)))
  g_scores <- lapply(scores, function(s) vapply(ids, function(g) mean(s$gromov[record, profiles[[g]]]), numeric(1)))
  d_top <- lapply(d_scores, function(v) ids[order(v)[1:top_k]])
  g_top <- lapply(g_scores, function(v) ids[order(v, decreasing = TRUE)[1:top_k]])
  d_common <- Reduce(intersect, d_top)
  g_common <- Reduce(intersect, g_top)

  disease_rows[[length(disease_rows) + 1L]] <- data.frame(
    row = record, disease_id = rownames(X)[record], disease = unname(meta$dis_lab[record]),
    annotations = sum(X[record, ]), profiles = length(ids),
    retention_distance = set_overlap(d_top), retention_gromov = set_overlap(g_top),
    common_distance = length(d_common), common_gromov = length(g_common),
    mean_overlap = mean(vapply(1:3, function(c) length(intersect(d_top[[c]], g_top[[c]])), integer(1)))
  )
  for (info in list(list("distance", d_common, d_scores), list("gromov", g_common, g_scores))) {
    for (g in info[[2]]) {
      ranks <- vapply(info[[3]], function(v) {
        local <- match(g, ids)
        rank(if (info[[1]] == "distance") v else -v, ties.method = "min")[local]
      }, numeric(1))
      consensus_rows[[length(consensus_rows) + 1L]] <- data.frame(
        row = record, disease_id = rownames(X)[record], disease = unname(meta$dis_lab[record]),
        ordering = info[[1]], profile_id = g, labels = labels[g], carriers = carriers[g],
        rank1 = ranks[1], rank2 = ranks[2], rank3 = ranks[3], mean_rank = mean(ranks)
      )
    }
  }
}
disease_summary <- do.call(rbind, disease_rows)
consensus       <- do.call(rbind, consensus_rows)

saveRDS(list(profiles = profiles, profile_id = profile_of, carriers = carriers, labels = labels,
             row_profiles = row_profiles, eligible = eligible, dis_sum = disease_summary),
        file.path(reported_dir, "s6_intermediate.rds"))
write.csv(full_by_chain,   file.path(reported_dir, "full_vocab_by_chain.csv"), row.names = FALSE)
write.csv(disease_summary, file.path(reported_dir, "annotation_internal_disease_summary.csv"), row.names = FALSE)
write.csv(consensus,       file.path(reported_dir, "annotation_internal_consensus.csv"), row.names = FALSE)

range_tag <- function(x, fmt) paste0(sprintf(fmt, min(x)), "-", sprintf(fmt, max(x)))
application_summary <- data.frame(
  diseases = nrow(X), phenotypes = ncol(X), annotations = sum(X),
  distinct_profiles = length(profiles), eligible_diseases = length(eligible),
  fullvocab_median_carriers_distance = range_tag(full_by_chain$median_carriers_distance, "%d"),
  fullvocab_median_carriers_gromov   = range_tag(full_by_chain$median_carriers_gromov, "%d"),
  fullvocab_share_le10_distance = range_tag(full_by_chain$share_le10_distance, "%.3f"),
  fullvocab_share_le10_gromov   = range_tag(full_by_chain$share_le10_gromov, "%.3f"),
  fullvocab_mean_shared = range_tag(full_by_chain$mean_overlap, "%.2f"),
  internal_median_retention_distance = median(disease_summary$retention_distance),
  internal_median_retention_gromov   = median(disease_summary$retention_gromov)
)


# 6. Choreoacanthocytosis case table (Table 3) ------------------
# Prespecified record; each ordering keeps up to three leading exact profiles
# whose cross-chain draw-wise top-ten support clears 0.70.

case_id  <- "OMIM:200150"
record   <- match(case_id, rownames(X))
screen   <- disease_summary[disease_summary$disease_id == case_id, ]
groups   <- row_profiles[[record]]
columns  <- unlist(profiles[groups], use.names = FALSE)
col_group <- rep(seq_along(groups), lengths(profiles[groups]))

per_chain <- lapply(1:3, function(chain) {
  fit <- fits[[chain]]
  n <- nrow(fit$theta); d <- dim(fit$z)[3]; p <- length(columns)
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

ordering_rows <- function(ordering) {
  means  <- lapply(per_chain, function(ch) colMeans(ch[[ordering]]))
  tops   <- lapply(means, function(v) order(v, decreasing = ordering == "gromov")[1:top_k])
  common <- Reduce(intersect, tops)
  rows <- do.call(rbind, lapply(common, function(g) {
    chain_ranks <- vapply(means, function(v)
      rank(if (ordering == "distance") v else -v, ties.method = "min")[g], numeric(1))
    support <- min(vapply(per_chain, function(ch) {
      v <- ch[[ordering]]
      mean(1L + rowSums(if (ordering == "distance") v < v[, g] else v > v[, g]) <= top_k)
    }, numeric(1)))
    data.frame(ordering = ordering, profile = labels[groups[g]], annotated_records = carriers[groups[g]],
               chain_ranks = paste(chain_ranks, collapse = "/"), mean_rank = mean(chain_ranks), support = support)
  }))
  rows <- rows[rows$support >= 0.70, ]
  head(rows[order(rows$mean_rank), ], 3L)
}

case <- rbind(ordering_rows("distance"), ordering_rows("gromov"))
rownames(case) <- NULL
write.csv(case, file.path(reported_dir, "selected_case_table.csv"), row.names = FALSE)

application_summary$selected_case               <- unname(meta$dis_lab[case_id])
application_summary$selected_case_id            <- case_id
application_summary$selected_profiles           <- screen$profiles
application_summary$selected_annotations        <- screen$annotations
application_summary$selected_retention_distance <- screen$retention_distance
application_summary$selected_retention_gromov   <- screen$retention_gromov
write.csv(application_summary, file.path(reported_dir, "application_summary.csv"), row.names = FALSE)


# 7. Table 3 and Figure 2 --------------------------------------

if (!nzchar(Sys.getenv("S6_SKIP_RENDER"))) {
  for (renderer in c("R/figures/table_dementia_case.R", "R/figures/figure_dementia_radial_tree.R")) {
    message("  render ", basename(renderer))
    if (system2("Rscript", shQuote(renderer), stdout = "", stderr = "") != 0L) stop("failed: ", renderer)
  }
}
message("Section 6 complete.")
