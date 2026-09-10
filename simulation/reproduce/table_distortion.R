## table_distortion.R — REPRODUCIBLE distortion table for the simulation slide.
## Reads the saved fits (fits/sim_s{seed}_{tag}.rds) and computes, per embedding (H^2, R^2..R^5):
##   distortion (multiplicative, optimal scale) and Spearman rank-fidelity to the tree distance,
##   plus the deep-band distortion (leaf-leaf pairs) for H^2 vs R^2.
## Outputs: reproduce/results/results_distortion.csv  and  figure/table_distortion.tex (\input-ready).
PRES <- Sys.getenv("LSIRM_PRES") # override: export LSIRM_PRES=/path/to/presentation
source(file.path(PRES, "reproduce/lib_sim.R"))
## A downloaded server bundle keeps completed fits at PRES/fits.  Older and
## live-run layouts use reproduce/fits.
downloaded_fits <- file.path(PRES, "fits")
if (dir.exists(downloaded_fits)) CFG$FITS <- normalizePath(downloaded_fits)
S <- readRDS(file.path(CFG$FITS, "sim_truth.rds"))
t <- build_tree()
tg <- t$G[S$items, S$items]
ut <- upper.tri(tg)
gg <- tg[ut]
leaf <- S$depi == max(S$depi)
deep <- outer(leaf, leaf, "&")[ut] # leaf-leaf pairs
dband <- function(D) {
  d <- D[ut][deep]
  g <- gg[deep]
  o <- g > 0 & is.finite(g) & is.finite(d) & d >= 0
  den <- sum(d[o]^2)
  if (!any(o) || den <= 0) stop("deep-band distortion is undefined")
  c <- sum(d[o] * g[o]) / den
  mean(abs(c * d[o] - g[o]) / g[o])
}

raw <- do.call(rbind, lapply(names(CFG$EMBED), function(tag) {
  e <- CFG$EMBED[[tag]]
  do.call(rbind, lapply(CFG$SEEDS, function(s) {
    W <- getW(readRDS(file.path(CFG$FITS, sprintf("sim_s%d_%s.rds", s, tag))))
    D <- per_draw_dist(W, e$hyp, CFG$SUB)
    data.frame(
      seed = s, embed = tag, mdim = e$mdim, hyp = e$hyp,
      distortion = distortion(D, tg), spearman = cor(D[ut], gg, method = "spearman"),
      deepband = dband(D), niter = CFG$NITER,
      stringsAsFactors = FALSE
    )
  }))
}))
dir.create(file.path(PRES, "reproduce/results"), recursive = TRUE, showWarnings = FALSE)
write.csv(raw, file.path(PRES, "reproduce/results/results_distortion_raw.csv"), row.names = FALSE)

rows <- lapply(names(CFG$EMBED), function(tag) {
  x <- raw[raw$embed == tag, ]
  data.frame(
    embed = tag, mdim = x$mdim[1], hyp = x$hyp[1],
    distortion = mean(x$distortion), distortion_se = sd(x$distortion) / sqrt(nrow(x)),
    spearman = mean(x$spearman), deepband = mean(x$deepband), nseed = nrow(x),
    niter = CFG$NITER, stringsAsFactors = FALSE
  )
})
df <- do.call(rbind, rows)
write.csv(df, file.path(PRES, "reproduce/results/results_distortion.csv"), row.names = FALSE)

cat(sprintf(
  "=== distortion table (depth%d BF%d, N=%d, NITER=%d, PARAM=%s, %d seeds) ===\n",
  CFG$DEPTH, CFG$BF, CFG$N, CFG$NITER, CFG$PARAM, length(CFG$SEEDS)
))
cat("embedding  distortion(SE)   Spearman   deepband\n")
for (i in 1:nrow(df)) {
  cat(sprintf(
    "%-9s  %.4f (%.4f)    %.4f      %.4f\n",
    df$embed[i], df$distortion[i], df$distortion_se[i], df$spearman[i], df$deepband[i]
  ))
}
cat(sprintf(
  "deep-band (leaf-leaf): H^2 %.2f  vs  R^2 %.2f\n",
  df$deepband[df$embed == "H2"], df$deepband[df$embed == "R2"]
))

## LaTeX snippet (matches the slide's table shape; \input{figure/table_distortion.tex})
texlab <- function(tag) {
  e <- CFG$EMBED[[tag]]
  sprintf("$\\mathbb{%s}^%d$", if (e$hyp) "H" else "R", e$mdim)
}
con <- file(file.path(CFG$TABLE, "table_distortion.tex"), "w")
writeLines(c(
  "\\begin{tabular}{lccc}", "\\toprule",
  "embedding & distortion $\\downarrow$ & Spearman $\\uparrow$ & leaf--leaf distortion $\\downarrow$ \\\\", "\\midrule"
), con)
best_dist <- which.min(df$distortion)
best_spear <- which.max(df$spearman)
best_deep <- which.min(df$deepband)
for (i in 1:nrow(df)) {
  writeLines(sprintf(
    "%s & %s & %s & %s \\\\", texlab(df$embed[i]),
    if (i == best_dist) sprintf("\\textbf{%.4f}", df$distortion[i]) else sprintf("%.4f", df$distortion[i]),
    if (i == best_spear) sprintf("\\textbf{%.4f}", df$spearman[i]) else sprintf("%.4f", df$spearman[i]),
    if (i == best_deep) sprintf("\\textbf{%.4f}", df$deepband[i]) else sprintf("%.4f", df$deepband[i])
  ), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), con)
close(con)
cat("saved reproduce/results/results_distortion_raw.csv + results_distortion.csv + figure/table_distortion.tex\n")
