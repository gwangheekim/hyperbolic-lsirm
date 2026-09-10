# Writes table/table_dementia_case.tex (manuscript Table 3) from
# result/dementia/reported/selected_case_table.csv. Run by R/main.R.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
reported_dir <- Sys.getenv("S6_REPORTED_DIR", unset = file.path(ROOT, "result/dementia/reported"))
case <- read.csv(file.path(reported_dir, "selected_case_table.csv"), stringsAsFactors = FALSE)

orderings <- c("distance", "gromov")
stopifnot(nrow(case) >= 2L, all(case$ordering %in% orderings),
          identical(case$ordering, case$ordering[order(match(case$ordering, orderings))]))

pretty <- c(distance = "Distance", gromov = "Gromov product")
label <- ifelse(!duplicated(case$ordering), pretty[case$ordering], "")

lines <- c(
  "\\begin{tabular}{p{0.13\\textwidth} p{0.38\\textwidth} c c c}",
  "\\toprule",
  "Ordering & Exact response profile & Annotated records & Chain ranks & Support \\\\",
  "\\midrule"
)
first_gromov <- match("gromov", case$ordering)
for (i in seq_len(nrow(case))) {
  if (!is.na(first_gromov) && i == first_gromov) lines <- c(lines, "\\addlinespace")
  lines <- c(lines, sprintf(
    "%s & %s & $%d$ & $%s$ & $%.3f$ \\\\",
    label[i], case$profile[i], case$annotated_records[i], case$chain_ranks[i], case$support[i]
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
dir.create(file.path(ROOT, "table"), showWarnings = FALSE)
writeLines(lines, file.path(ROOT, "table", "table_dementia_case.tex"))
