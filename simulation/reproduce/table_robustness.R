## Build the paper's combined robustness table from the reported simulations.
## Numerical values are never typed into this script.
PRES <- Sys.getenv("LSIRM_PRES")
if (!nzchar(PRES)) stop("LSIRM_PRES is not set")
RDIR <- file.path(PRES, "reproduce", "results")
A <- read.csv(file.path(RDIR, "results_robustness.csv"), stringsAsFactors = FALSE)
B <- read.csv(file.path(RDIR, "results_n_sweep.csv"), stringsAsFactors = FALSE)
C <- read.csv(file.path(RDIR, "results_negative_control.csv"), stringsAsFactors = FALSE)

stopifnot(nrow(A) == 9L, identical(sort(B$N), c(500L, 1000L, 1500L, 2000L, 4000L)), nrow(C) == 2L)
labA <- c(
  "balanced BF2 depth4 (baseline)" = "balanced, depth $4$, BF $2$ (baseline, $M=30$)",
  "balanced BF3 depth3" = "balanced, depth $3$, BF $3$ ($M=39$)",
  "balanced BF5 depth2 (bushy)" = "balanced, depth $2$, BF $5$, bushy ($M=30$)",
  "unbalanced shallow 5/23" = "unbalanced shallow, two branches with $5/23$ leaves ($M=30$)",
  "unbalanced mixed depth 2/3/4/6" = "unbalanced depth, branch depths $2/3/4/6$ ($M=30$)",
  "balanced BF2 row SD 1" = "depth $4$, BF $2$, $\\sigma_{\\theta_0}=1$",
  "balanced BF2 row SD 4" = "depth $4$, BF $2$, $\\sigma_{\\theta_0}=4$",
  "unbalanced random" = "unbalanced random growth, depth ceiling $5$",
  "unbalanced random deeper" = "unbalanced random growth, depth ceiling $6$"
)
labC <- c(
  "flat 2D plane" = "flat plane ($\\mathbb{R}^2$)",
  "flat 1D line" = "1D continuum (ability-test regime)"
)
if (anyNA(labA[A$family]) || anyNA(labC[C$truth])) stop("unexpected condition label")

eol <- "\\\\"
rowtex <- function(label, delta, h, r, gap, se) {
  paste0(sprintf(
    "\\quad %s & %.4f & %.4f & %.4f & $%+.4f{\\scriptstyle\\,\\pm%.4f}$ ",
    label, delta, h, r, gap, se
  ), eol)
}

lines <- c(
  "\\begin{tabular}{lrrrr}", "\\toprule",
  paste0("Setting & $\\delta/\\mathrm{diam}$ & $\\mathbb H^2$ $\\downarrow$ & $\\mathbb R^2$ $\\downarrow$ & gap ($\\mathbb R^2-\\mathbb H^2$) ", eol),
  "\\midrule", paste0("\\multicolumn{5}{l}{A. Generating tree family, $N=1500$} ", eol)
)
for (i in seq_len(nrow(A))) lines <- c(lines, rowtex(labA[A$family[i]], A$truth_delta[i], A$H2[i], A$R2[i], A$gap[i], A$gap_se[i]))
lines <- c(lines, "\\addlinespace", paste0("\\multicolumn{5}{l}{B. Number of respondents, baseline tree ($M=30$)} ", eol))
for (i in seq_len(nrow(B))) {
  lbl <- sprintf("$N=%s$%s", format(B$N[i], big.mark = ",", scientific = FALSE), if (B$N[i] == 1500L) " (default)" else "")
  lines <- c(lines, rowtex(lbl, B$truth_delta[i], B$dist_h[i], B$dist_r[i], B$gap[i], B$gap_se[i]))
}
lines <- c(lines, "\\addlinespace", paste0("\\multicolumn{5}{l}{C. Truths that do not branch} ", eol))
for (i in seq_len(nrow(C))) lines <- c(lines, rowtex(labC[C$truth[i]], C$truth_delta[i], C$H2[i], C$R2[i], C$gap[i], C$gap_se[i]))
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(CFG$TABLE, "table_robustness.tex"))
cat("saved figure/table_robustness.tex\n")
