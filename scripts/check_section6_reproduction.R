# Compares two sets of Section 6 artefacts and checks the manuscript numbers.
#
#   Rscript scripts/check_section6_reproduction.R <fresh_dir> [<reference_dir>]
#
# <fresh_dir>      directory holding a freshly regenerated set
# <reference_dir>  optional; the committed set to diff against
#
# With one argument only the manuscript-number checks run, against <fresh_dir>.

args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- normalizePath(file.path(dirname(gsub("~+~", " ", script_path, fixed = TRUE)), ".."), mustWork = TRUE)
stopifnot(length(args) >= 1L)
FRESH <- normalizePath(args[1], mustWork = TRUE)
REF <- if (length(args) >= 2L) normalizePath(args[2], mustWork = TRUE) else NA_character_

pass <- TRUE
note <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "OK" else "XX", msg))
  if (!isTRUE(ok)) pass <<- FALSE
}
have <- function(dir, name) file.exists(file.path(dir, name))

compare_csv <- function(name, tol = 1e-6) {
  if (!have(FRESH, name) || !have(REF, name)) {
    cat(sprintf("  [--] %s absent on one side, skipped\n", name)); return(invisible())
  }
  a <- read.csv(file.path(FRESH, name), stringsAsFactors = FALSE, check.names = FALSE)
  b <- read.csv(file.path(REF, name), stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
    note(FALSE, sprintf("%s shape/columns differ (%s vs %s)", name,
                        paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x"))); return(invisible())
  }
  worst <- 0; textbad <- character()
  for (col in names(a)) {
    if (is.numeric(a[[col]]) && is.numeric(b[[col]])) {
      worst <- max(worst, max(abs(a[[col]] - b[[col]]), na.rm = TRUE))
    } else if (!identical(as.character(a[[col]]), as.character(b[[col]]))) {
      textbad <- c(textbad, col)
    }
  }
  if (length(textbad)) note(FALSE, sprintf("%s text columns differ: %s", name, paste(textbad, collapse = ", ")))
  else note(worst <= tol, sprintf("%s  max|Δ| = %.2e", name, worst))
}

compare_rds_field <- function(name, field, tol = 1e-6) {
  if (!have(FRESH, name) || !have(REF, name)) {
    cat(sprintf("  [--] %s absent on one side, skipped\n", name)); return(invisible())
  }
  a <- readRDS(file.path(FRESH, name)); b <- readRDS(file.path(REF, name))
  av <- a[[field]]; bv <- b[[field]]
  d <- max(abs(as.numeric(av) - as.numeric(bv)))
  note(d <= tol, sprintf("%s$%s  max|Δ| = %.2e", name, field, d))
}

if (!is.na(REF)) {
  cat("== fresh vs reference ==\n")
  compare_rds_field("dementia_f3_geometry.rds", "distance")
  compare_rds_field("dementia_f3_geometry.rds", "radius")
  compare_rds_field("dementia_f3_geometry.rds", "gromov")
  for (f in c("full_vocab_by_chain.csv", "annotation_internal_disease_summary.csv",
              "annotation_internal_consensus.csv", "selected_case_table.csv",
              "rank_stability_summary.csv", "dementia_groups.csv", "application_summary.csv")) {
    compare_csv(f)
  }
  cat("\n")
}

cat("== manuscript Section 6 numbers (fresh set) ==\n")
meta <- readRDS(file.path(ROOT, "data/analysis_data/dementia/hpo_dementia_f3_meta.rds"))
note(identical(dim(meta$X), c(170L, 313L)) && sum(meta$X) == 2884L, "170 x 313 matrix, 2884 annotations")

summ <- read.csv(file.path(FRESH, "application_summary.csv"), check.names = FALSE, colClasses = "character")
note(summ$distinct_profiles == "303", sprintf("303 distinct exact response profiles (got %s)", summ$distinct_profiles))
note(summ$eligible_diseases == "134", sprintf("134 records with >= 10 profiles (got %s)", summ$eligible_diseases))
note(all(strsplit(summ$fullvocab_median_carriers_distance, "-")[[1]] == "18"),
     sprintf("median leading-ten carriers 18 under distance (got %s)", summ$fullvocab_median_carriers_distance))
note(all(strsplit(summ$fullvocab_median_carriers_gromov, "-")[[1]] == "6"),
     sprintf("median leading-ten carriers 6 under Gromov (got %s)", summ$fullvocab_median_carriers_gromov))
note(summ$selected_case_id == "OMIM:200150" && summ$selected_profiles == "21" && summ$selected_annotations == "21",
     sprintf("choreoacanthocytosis: 21 phenotypes / 21 profiles (got %s / %s)",
             summ$selected_profiles, summ$selected_annotations))

case <- read.csv(file.path(FRESH, "selected_case_table.csv"), stringsAsFactors = FALSE)
expect_case <- data.frame(
  ordering = c("distance", "distance", "distance", "gromov", "gromov", "gromov"),
  profile = c("Depression", "Anxiety", "Chorea",
              "Increased circulating lactate dehydrogenase concentration", "Acanthocytosis", "Anxiety"),
  annotated_records = c(29L, 12L, 20L, 4L, 5L, 12L),
  chain_ranks = c("1/1/1", "2/2/2", "3/3/5", "1/1/1", "2/2/2", "3/3/3"),
  support = c(.949, .834, .722, .933, .937, .931),
  stringsAsFactors = FALSE
)
ok_rows <- nrow(case) == 6L && identical(case$ordering, expect_case$ordering) &&
  identical(case$profile, expect_case$profile) &&
  identical(as.integer(case$annotated_records), expect_case$annotated_records)
note(ok_rows, "Table 3 profiles and carrier counts")
if (ok_rows) {
  note(identical(case$chain_ranks, expect_case$chain_ranks), "Table 3 chain ranks")
  note(identical(round(case$support, 3), expect_case$support),
       sprintf("Table 3 support values (got %s)", paste(sprintf("%.3f", case$support), collapse = " ")))
}

if (have(FRESH, "dementia_groups.csv")) {
  groups <- read.csv(file.path(FRESH, "dementia_groups.csv"), stringsAsFactors = FALSE)
  note(nrow(groups) == 170L && length(unique(groups$group)) == 8L,
       sprintf("radial tree: 170 records in %d groups", length(unique(groups$group))))
}

cat(sprintf("\n%s\n", if (pass) "SECTION 6 REPRODUCTION PASS" else "SECTION 6 REPRODUCTION FAIL"))
if (!pass) quit(status = 1L)
