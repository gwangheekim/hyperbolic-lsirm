args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

sim <- read.csv(file.path(ROOT, "result/simulation/reported/results_distortion.csv"), stringsAsFactors = FALSE)
expected_distortion <- data.frame(
  embed = c("H2", "H3", "R2", "R3", "R4", "R5"),
  distortion = c(.1526, .1752, .2571, .2229, .2007, .1923),
  spearman = c(.9567, .9486, .8937, .9236, .9310, .9329),
  deepband = c(.1241, .1282, .2757, .2283, .2067, .1986)
)
stopifnot(identical(sim$embed, expected_distortion$embed))
for (name in c("distortion", "spearman", "deepband")) {
  stopifnot(identical(round(sim[[name]], 4L), expected_distortion[[name]]))
}

robustness <- read.csv(file.path(ROOT, "result/simulation/reported/results_robustness.csv"), stringsAsFactors = FALSE)
expected_robustness <- data.frame(
  H2 = c(.1526, .1725, .2082, .4734, .1451, .1732, .1428, .2117, .1554),
  R2 = c(.2571, .3060, .3199, .5352, .2812, .2580, .2484, .2974, .2264),
  gap = c(.1045, .1335, .1118, .0618, .1361, .0848, .1056, .0857, .0710)
)
stopifnot(nrow(robustness) == nrow(expected_robustness))
for (name in names(expected_robustness)) {
  stopifnot(identical(round(robustness[[name]], 4L), expected_robustness[[name]]))
}

n_sweep <- read.csv(file.path(ROOT, "result/simulation/reported/results_n_sweep.csv"), stringsAsFactors = FALSE)
expected_n_sweep <- data.frame(
  dist_h = c(.1736, .1627, .1526, .1438, .1289),
  dist_r = c(.2363, .2488, .2571, .2606, .2684),
  gap = c(.0626, .0861, .1045, .1168, .1395)
)
for (name in names(expected_n_sweep)) {
  stopifnot(identical(round(n_sweep[[name]], 4L), expected_n_sweep[[name]]))
}

negative_control <- read.csv(file.path(ROOT, "result/simulation/reported/results_negative_control.csv"), stringsAsFactors = FALSE)
expected_negative_control <- data.frame(
  H2 = c(.2115, .8015), R2 = c(.1007, .4572), gap = c(-.1108, -.3443)
)
for (name in names(expected_negative_control)) {
  stopifnot(identical(round(negative_control[[name]], 4L), expected_negative_control[[name]]))
}

meta_path <- file.path(ROOT, "data/analysis_data/dementia/hpo_dementia_f3_meta.rds")
meta <- readRDS(meta_path)
stopifnot(
  identical(unname(tools::md5sum(meta_path)), "30388b671cc7137cb5a656e76b2633b4"),
  identical(dim(meta$X), c(170L, 313L)), sum(meta$X) == 2884L
)
profiles <- apply(meta$X, 2L, paste0, collapse = "")
qualifying <- vapply(seq_len(nrow(meta$X)), function(i) {
  length(unique(profiles[meta$X[i, ] == 1L]))
}, integer(1))
stopifnot(sum(qualifying >= 10L) == 134L)

fit_dir <- file.path(ROOT, "result/model_fits/dementia")
for (chain in 1:3) {
  fit <- readRDS(file.path(fit_dir, sprintf("H3_chain%d.rds", chain)))
  stopifnot(
    identical(fit$data_md5, "30388b671cc7137cb5a656e76b2633b4"),
    identical(dim(fit$theta), c(9000L, 170L)),
    identical(dim(fit$w), c(9000L, 313L, 4L)),
    isTRUE(fit$polar_split), isTRUE(fit$stable_intrinsic_proposal)
  )
}

case <- read.csv(file.path(ROOT, "result/dementia/reported/selected_case_table.csv"), stringsAsFactors = FALSE)
stopifnot(
  identical(case$ordering, c("distance", "distance", "distance", "gromov", "gromov", "gromov")),
  identical(case$profile, c(
    "Depression", "Anxiety", "Chorea",
    "Increased circulating lactate dehydrogenase concentration", "Acanthocytosis", "Anxiety"
  )),
  identical(case$annotated_records, c(29L, 12L, 20L, 4L, 5L, 12L)),
  identical(case$chain_ranks, c("1/1/1", "2/2/2", "3/3/5", "1/1/1", "2/2/2", "3/3/3")),
  identical(round(case$support, 3L), c(.949, .834, .722, .933, .937, .931)),
  min(case$support) >= .7
)

groups <- read.csv(file.path(ROOT, "result/dementia/reported/dementia_groups.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(groups) == 170L, length(unique(groups$group)) == 8L)

required <- file.path(ROOT, c(
  "figure/fig_sim_negcontrol_truths.pdf",
  "figure/fig_dementia_radial_tree.pdf",
  "table/table_distortion.tex",
  "table/table_robustness.tex",
  "table/table_dementia_case.tex"
))
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0L))
cat("BIOMETRICS SUBMISSION VALIDATION PASS\n")
