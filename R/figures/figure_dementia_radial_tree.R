# Writes figure/fig_dementia_radial_tree.pdf (manuscript Figure 2) and
# result/dementia/reported/dementia_groups.csv (Web Appendix D membership) from
# the dementia geometry and the three H3 chains. The dendrogram is Ward.D2 on
# the posterior-mean pairwise hyperbolic distance, cut at the silhouette-optimal
# number of groups over 2..12; each group is labelled by its signature phenotype
# (largest mean posterior Gromov product to the members). Run by
# R/main.R.

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
})
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "figure")
reported_dir <- Sys.getenv("S6_REPORTED_DIR", unset = file.path(ROOT, "result/dementia/reported"))
fit_dir <- Sys.getenv("FIT_OUTPUT_DIR", unset = file.path(ROOT, "result/model_fits/dementia"))
meta <- readRDS(file.path(ROOT, "data/analysis_data/dementia/hpo_dementia_f3_meta.rds"))
geo <- readRDS(file.path(reported_dir, "dementia_f3_geometry.rds"))
X <- meta$X
lab <- meta$lab
dn <- meta$dis_lab
dis <- rownames(X)
Nd <- nrow(X)
ann <- rowSums(X)
# Phenotype nodes whose label repeats a disease name are excluded from signature-phenotype selection; without this C1 and C2 would self-label.
selfref <- grepl("^(Alzheimer disease|Amyotrophic lateral sclerosis|Frontotemporal dementia)([.][0-9]+)?$", lab)
hc <- hclust(as.dist(geo$distance), "ward.D2")
ph <- ladderize(as.phylo(hc))
# Seven groups: the cut in 2..12 that maximizes the mean silhouette width. The cut and the groups are descriptive only.
sil <- function(k) mean(silhouette(cutree(hc, k), as.dist(geo$distance))[, "sil_width"])
K <- which.max(vapply(2:12, sil, numeric(1))) + 1L
cl <- cutree(hc, K)

SG <- matrix(0, Nd, ncol(X))
ndraw <- 0L
for (chain in 1:3) {
  fit <- readRDS(file.path(fit_dir, sprintf("H3_chain%d.rds", chain)))
  nd <- dim(fit$z)[1]
  for (t in seq_len(nd)) {
    Z <- fit$z[t, , ]
    W <- fit$w[t, , ]
    rz <- asinh(sqrt(rowSums(Z[, -1, drop = FALSE]^2)))
    rw <- asinh(sqrt(rowSums(W[, -1, drop = FALSE]^2)))
    SG <- SG + (outer(rz, rw, "+") - acosh(pmax(Z[, 1] %o% W[, 1] - Z[, -1, drop = FALSE] %*% t(W[, -1, drop = FALSE]), 1))) / 2
  }
  ndraw <- ndraw + nd
  rm(fit)
  gc()
}
G2 <- 2 * SG / ndraw

descr <- function(i) {
  nm <- tolower(paste(dn[i], collapse = " | "))
  kw <- c(
    "spinocerebellar ataxia" = "Spinocerebellar ataxia",
    "frontotemporal|amyotrophic lateral|pick disease|semantic dementia|primary progressive aphasia" = "Frontotemporal / ALS spectrum",
    "alzheimer" = "Alzheimer-type cortical dementia",
    "parkinson|parkinsonism|kufor|pallidopyramidal" = "Parkinsonism",
    "leukodystrophy|adrenoleuko|leukoenceph" = "Leukodystrophy",
    "acanthocyt|neuroferritin|pantothenate|brain iron|aceruloplasm|neurodegeneration with brain iron" = "Neuroacanthocytosis / brain-iron",
    "choreo|huntington|dentatorubral|chorea" = "Choreatic"
  )
  hit <- names(kw)[sapply(names(kw), function(p) grepl(p, nm))]
  frac <- sapply(hit, function(p) mean(grepl(strsplit(p, "\\|")[[1]][1], tolower(dn[i]))))
  if (length(hit) && max(frac) >= 0.45) {
    unname(kw[hit[which.max(frac)]])
  } else {
    e <- dn[i][which.max(ann[i])]
    if (nchar(e) > 30) paste0(substr(e, 1, 29), "..") else e
  }
}
# Integrating clinical labels, assigned by reading each group's disease
# membership (Web Appendix D) and keyed to the group's signature phenotype so
# the script is deterministic. Each label covers the group's dominant category
# together with the related disorders it also contains.
nm_by_hall <- c(
  "neurofibrillary tangles" = "Alzheimer disease",
  "mutism" = "Frontotemporal degeneration and ALS",
  "involuntary movements" = "Brain-iron and basal-ganglia degeneration",
  "chorea" = "Brain-iron and basal-ganglia degeneration",
  "compulsive behaviors" = "Brain-iron and basal-ganglia degeneration",
  "cerebellar hypoplasia" = "Early-onset cerebellar ataxia",
  "decreased nerve conduction velocity" = "Leukodystrophy and neuropathy",
  "decerebrate rigidity" = "Leukodystrophy and neuropathy",
  "dysmetria" = "Early-onset cerebellar ataxia",
  "limb ataxia" = "Dominant spinocerebellar ataxia",
  "gait ataxia" = "Dominant spinocerebellar ataxia",
  "gaze-evoked horizontal nystagmus" = "Dominant spinocerebellar ataxia",
  "hypomimic face" = "Parkinsonism and dystonia-parkinsonism",
  "abnormality of temperature regulation" = "Machado-Joseph and related ataxia",
  "progressive external ophthalmoplegia" = "Machado-Joseph and related ataxia"
)
grp <- lapply(1:K, function(g) {
  i <- which(cl == g)
  v <- colMeans(G2[i, , drop = FALSE])
  v[selfref] <- -Inf
  h <- tolower(lab[which.max(v)])
  list(g = g, n = length(i), name = if (h %in% names(nm_by_hall)) nm_by_hall[[h]] else descr(i), hall = h)
})


expected_hall <- c(
  "neurofibrillary tangles", "mutism", "compulsive behaviors", "gait ataxia",
  "cerebellar hypoplasia", "hypomimic face", "progressive external ophthalmoplegia",
  "decerebrate rigidity"
)
realized_hall <- vapply(grp, function(x) x$hall, character(1))
if (!setequal(realized_hall, expected_hall)) {
  stop(
    "signature phenotypes differ from the reported set; update nm_by_hall, ",
    "expected_hall, and the figure caption:\n",
    paste(sprintf("  C%d -> %s", seq_len(K), realized_hall), collapse = "\n")
  )
}

pal <- c("#6a4fb3", "#c0397f", "#2f8f78", "#c46210", "#4c7abf", "#5b6770", "#b58a00", "#8f2d56")[1:K]
tcl <- cl[ph$tip.label]
tcol <- pal[tcl]
desc <- prop.part(ph)
ecol <- rep("#bbb", nrow(ph$edge))
ew <- rep(0.7, nrow(ph$edge))
for (e in seq_len(nrow(ph$edge))) {
  nd2 <- ph$edge[e, 2]
  ti <- if (nd2 <= Ntip(ph)) nd2 else desc[[nd2 - Ntip(ph)]]
  gg <- unique(tcl[ph$tip.label[ti]])
  if (length(gg) == 1) {
    ecol[e] <- pal[gg]
    ew[e] <- if (nd2 <= Ntip(ph)) 1.3 else 1.1
  }
}

out_pdf <- file.path(OUT, "fig_dementia_radial_tree.pdf")
pdf(out_pdf, width = 11.0417, height = 7.25, useDingbats = FALSE)
par(mar = c(1.5, 11.5, 1.5, 11.5), xpd = NA)
plot(ph, type = "fan", show.tip.label = FALSE, edge.color = ecol, edge.width = ew, open.angle = 10)
pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
RM <- max(sqrt(pp$xx[1:Nd]^2 + pp$yy[1:Nd]^2))
points(pp$xx[1:Nd], pp$yy[1:Nd], pch = 19, cex = 0.55, col = tcol)
ang <- atan2(pp$yy[1:Nd], pp$xx[1:Nd])
for (g in 1:K) {
  i <- which(tcl[ph$tip.label] == g)
  aa <- ang[i]
  ma <- atan2(mean(sin(aa)), mean(cos(aa)))
  aord <- order((aa - ma + pi) %% (2 * pi))
  lines(RM * 1.05 * cos(aa[aord]), RM * 1.05 * sin(aa[aord]), col = pal[g], lwd = 3)
  gi <- grp[[g]]
  text(RM * 1.11 * cos(ma), RM * 1.11 * sin(ma),
    sprintf("%s  (n=%d)\n-> %s", gi$name, gi$n, gi$hall),
    cex = 0.82, col = pal[g], font = 2, adj = c(ifelse(cos(ma) >= 0, 0, 1), .5)
  )
}
dev.off()

pdfcrop <- Sys.which("pdfcrop")
if (nzchar(pdfcrop)) {
  status <- system2(pdfcrop, c("--margins", "6", shQuote(out_pdf), shQuote(out_pdf)),
    stdout = FALSE, stderr = FALSE)
  if (status != 0L) warning("pdfcrop failed; figure kept at full device size")
} else {
  warning("pdfcrop not found; figure kept at full device size with uneven margins")
}

cat(sprintf("K=%d\n", K))
for (g in 1:K) cat(sprintf("  C%d  n=%3d  %-34s  -> %s\n", g, grp[[g]]$n, grp[[g]]$name, grp[[g]]$hall))

# Record-level group membership, written next to the other reported artefacts
# as the replication file for Web Appendix D.
REP <- reported_dir
members <- data.frame(
  disease_id = dis,
  disease_label = dn,
  group = cl,
  group_name = vapply(cl, function(g) grp[[g]]$name, character(1)),
  annotations = ann,
  stringsAsFactors = FALSE
)
members <- members[order(members$group, -members$annotations), ]
write.csv(members, file.path(REP, "dementia_groups.csv"), row.names = FALSE)
