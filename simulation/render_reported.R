## Render the two manuscript simulation tables and Figure 1 from reported results/specification.

args <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]), fixed = TRUE), mustWork = TRUE)
SIM <- normalizePath(dirname(script_file), mustWork = TRUE)
ROOT <- normalizePath(file.path(SIM, ".."), mustWork = TRUE)
setwd(ROOT)
dir.create("figure", showWarnings = FALSE)
dir.create("table", showWarnings = FALSE)
RDIR <- "result/simulation/reported"

D <- read.csv(file.path(RDIR, "results_distortion.csv"), stringsAsFactors = FALSE)
eol <- "\\\\"
texlab <- function(tag, d) sprintf("$\\mathbb{%s}^%d$", if (d$hyp[tag]) "H" else "R", d$mdim[tag])
dd <- D
rownames(dd) <- dd$embed
best_dist <- which.min(D$distortion)
best_spear <- which.max(D$spearman)
best_deep <- which.min(D$deepband)
lines <- c(
  "\\begin{tabular}{lccc}", "\\toprule",
  "embedding & distortion $\\downarrow$ & Spearman $\\uparrow$ & leaf--leaf distortion $\\downarrow$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(D))) {
  lab <- sprintf("$\\mathbb{%s}^%d$", if (D$hyp[i]) "H" else "R", D$mdim[i])
  fmt <- function(x, bold) if (bold) sprintf("\\textbf{%.4f}", x) else sprintf("%.4f", x)
  lines <- c(lines, sprintf(
    "%s & %s & %s & %s %s", lab,
    fmt(D$distortion[i], i == best_dist), fmt(D$spearman[i], i == best_spear),
    fmt(D$deepband[i], i == best_deep), eol
  ))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), "table/table_distortion.tex")

A <- read.csv(file.path(RDIR, "results_robustness.csv"), stringsAsFactors = FALSE)
B <- read.csv(file.path(RDIR, "results_n_sweep.csv"), stringsAsFactors = FALSE)
C <- read.csv(file.path(RDIR, "results_negative_control.csv"), stringsAsFactors = FALSE)
labA <- c(
  "balanced BF2 depth4 (baseline)" = "balanced, depth $4$, BF $2$ (baseline, $M=30$)",
  "balanced BF3 depth3" = "balanced, depth $3$, BF $3$ ($M=39$)",
  "balanced BF5 depth2 (bushy)" = "balanced, depth $2$, BF $5$, bushy ($M=30$)",
  "unbalanced shallow 5/23" = "unbalanced shallow, leaf split $5/23$ ($M=30$)",
  "unbalanced mixed depth 2/3/4/6" = "unbalanced, terminal depths $2/3/4/6$ ($M=30$)",
  "balanced BF2 row SD 1" = "depth $4$, BF $2$, $\\sigma_{\\theta_0}=1$",
  "balanced BF2 row SD 4" = "depth $4$, BF $2$, $\\sigma_{\\theta_0}=4$",
  "unbalanced random" = "unbalanced random growth, depth ceiling $5$",
  "unbalanced random deeper" = "unbalanced random growth, depth ceiling $6$"
)
labC <- c(
  "flat 2D plane" = "flat plane ($\\mathbb{R}^2$)",
  "flat 1D line" = "1D continuum (ability-test regime)"
)
rowtex <- function(label, delta, h, r, gap, se) {
  sprintf(
    "\\quad %s & %.4f & %.4f & %.4f & $%+.4f{\\scriptstyle\\,\\pm%.4f}$ %s",
    label, delta, h, r, gap, se, eol
  )
}
lines <- c(
  "\\begin{tabular}{lrrrr}", "\\toprule",
  paste0("Setting & $\\delta/\\mathrm{diam}$ & $\\mathbb H^2$ $\\downarrow$ & $\\mathbb R^2$ $\\downarrow$ & gap ($\\mathbb R^2-\\mathbb H^2$) ", eol),
  "\\midrule", paste0("\\multicolumn{5}{l}{A. Generating tree family, $N=1500$} ", eol)
)
for (i in seq_len(nrow(A))) lines <- c(lines, rowtex(labA[A$family[i]], A$truth_delta[i], A$H2[i], A$R2[i], A$gap[i], A$gap_se[i]))
lines <- c(lines, "\\addlinespace", paste0("\\multicolumn{5}{l}{B. Number of respondents, baseline tree ($M=30$)} ", eol))
for (i in seq_len(nrow(B))) {
  label <- sprintf("$N=%s$%s", format(B$N[i], big.mark = ",", scientific = FALSE), if (B$N[i] == 1500L) " (default)" else "")
  lines <- c(lines, rowtex(label, B$truth_delta[i], B$dist_h[i], B$dist_r[i], B$gap[i], B$gap_se[i]))
}
lines <- c(lines, "\\addlinespace", paste0("\\multicolumn{5}{l}{C. Truths that do not branch} ", eol))
for (i in seq_len(nrow(C))) lines <- c(lines, rowtex(labC[C$truth[i]], C$truth_delta[i], C$H2[i], C$R2[i], C$gap[i], C$gap_se[i]))
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), "table/table_robustness.tex")

set.seed(4)
M <- 30L
pal <- colorRampPalette(c("#9ecae1", "#08306b"))(7)
edge <- "#2b6a4d"
finish_tree <- function(par) {
  dep <- integer(length(par))
  for (i in 2:length(par)) dep[i] <- dep[par[i]] + 1L
  list(par = par, dep = dep)
}
balanced_tree <- function(depth, bf) {
  par <- 0L
  add <- function(node, level) {
    if (level >= depth) {
      return()
    }
    for (b in seq_len(bf)) {
      par <<- c(par, node)
      add(length(par), level + 1L)
    }
  }
  add(1L, 0L)
  finish_tree(par)
}
unbalanced_5_23 <- function() finish_tree(c(0L, 1L, 1L, rep(2L, 5L), rep(3L, 23L)))
unbalanced_mixed <- function() {
  par <- 0L
  spines <- vector("list", 4L)
  depths <- c(2L, 3L, 4L, 6L)
  for (a in seq_along(depths)) {
    parent <- 1L
    for (d in seq_len(depths[a])) {
      par <- c(par, parent)
      parent <- length(par)
      spines[[a]] <- c(spines[[a]], parent)
    }
  }
  for (k in seq_len(4L)) par <- c(par, spines[[1]][1])
  for (k in seq_len(4L)) par <- c(par, spines[[2]][2])
  for (k in seq_len(4L)) par <- c(par, spines[[3]][3])
  for (k in seq_len(3L)) par <- c(par, spines[[4]][5])
  stopifnot(length(par) - 1L == M)
  finish_tree(par)
}
tree_layout <- function(tree) {
  kids <- lapply(seq_along(tree$par), function(i) which(tree$par == i))
  x <- rep(NA_real_, length(tree$par))
  leaf <- 0L
  place <- function(v) {
    if (!length(kids[[v]])) {
      leaf <<- leaf + 1L
      x[v] <<- leaf
    } else {
      for (child in kids[[v]]) place(child)
      x[v] <<- mean(x[kids[[v]]])
    }
  }
  place(1L)
  list(x = x, y = max(tree$dep) - tree$dep)
}
panel_tree <- function(tree, main) {
  xy <- tree_layout(tree)
  plot(NA, xlim = range(xy$x) + c(-.7, .7), ylim = range(xy$y) + c(-.45, .45), axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 1.15, font.main = 1)
  for (k in 2:length(tree$par)) segments(xy$x[tree$par[k]], xy$y[tree$par[k]], xy$x[k], xy$y[k], col = edge, lwd = 1)
  points(xy$x[-1], xy$y[-1], pch = 21, cex = 1.05, lwd = .5, bg = pal[pmin(tree$dep[-1] + 1L, length(pal))], col = "white")
  points(xy$x[1], xy$y[1], pch = 21, cex = 1.05, lwd = .5, bg = "grey75", col = "white")
}
pdf("figure/fig_sim_negcontrol_truths.pdf", width = 12, height = 7.1)
par(mfrow = c(2, 3), mar = c(1.6, 1.2, 2.8, 1.2), xpd = NA)
panel_tree(balanced_tree(4L, 2L), "(A) balanced binary tree, depth 4")
panel_tree(balanced_tree(2L, 5L), "(B) balanced bushy tree, depth 2")
panel_tree(unbalanced_5_23(), "(C) unbalanced shallow tree\n5/23 split")
panel_tree(unbalanced_mixed(), "(D) unbalanced mixed-depth tree\nterminal depths 2/3/4/6")
xy <- matrix(runif(2 * M), ncol = 2)
plot(NA, xlim = c(-.06, 1.06), ylim = c(-.06, 1.06), asp = 1, axes = FALSE, xlab = "", ylab = "", main = "(E) flat plane", cex.main = 1.15, font.main = 1)
rect(0, 0, 1, 1, border = "grey55", lwd = .9)
points(xy, pch = 21, cex = 1.15, lwd = .5, bg = "grey55", col = "white")
u <- sort(runif(M))
plot(NA, xlim = c(-.06, 1.06), ylim = c(-.35, .35), axes = FALSE, xlab = "", ylab = "", main = "(F) one-dimensional continuum", cex.main = 1.15, font.main = 1)
segments(0, 0, 1, 0, col = "grey55", lwd = .9)
points(u, rep(0, M), pch = 21, cex = 1.15, lwd = .5, bg = colorRampPalette(c("#9ecae1", "#08306b"))(M), col = "white")
text(0, -.20, "easiest", cex = .85, col = "grey30")
text(1, -.20, "hardest", cex = .85, col = "grey30")
arrows(.08, -.13, .92, -.13, length = .07, col = "grey45", lwd = 1)
dev.off()
cat("rendered reported simulation tables and Figure 1\n")
