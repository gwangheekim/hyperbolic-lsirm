args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
PKG <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
set.seed(4)
M <- 30L
pal <- colorRampPalette(c("#9ecae1", "#08306b"))(7)
EDGE <- "#2b6a4d"

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

unbalanced_5_23 <- function() {
  par <- c(0L, 1L, 1L)
  par <- c(par, rep(2L, 5L), rep(3L, 23L))
  finish_tree(par)
}

unbalanced_mixed <- function() {
  par <- 0L
  spines <- vector("list", 4L)
  for (a in seq_along(c(2L, 3L, 4L, 6L))) {
    parent <- 1L
    spines[[a]] <- integer(0)
    for (d in seq_len(c(2L, 3L, 4L, 6L)[a])) {
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
  par <- tree$par
  dep <- tree$dep
  n <- length(par)
  kids <- lapply(seq_len(n), function(i) which(par == i))
  x <- rep(NA_real_, n)
  leaf <- 0L
  place <- function(v) {
    if (!length(kids[[v]])) {
      leaf <<- leaf + 1L
      x[v] <<- leaf
    } else {
      for (ch in kids[[v]]) place(ch)
      x[v] <<- mean(x[kids[[v]]])
    }
  }
  place(1L)
  list(x = x, y = max(dep) - dep)
}

panel_tree <- function(tree, main) {
  xy <- tree_layout(tree)
  par <- tree$par
  dep <- tree$dep
  plot(NA,
    xlim = range(xy$x) + c(-.7, .7), ylim = range(xy$y) + c(-.45, .45),
    axes = FALSE, xlab = "", ylab = "", main = main, cex.main = 1.28, font.main = 1
  )
  for (k in 2:length(par)) {
    segments(xy$x[par[k]], xy$y[par[k]], xy$x[k], xy$y[k], col = EDGE, lwd = 1.0)
  }
  points(xy$x[-1], xy$y[-1],
    pch = 21, cex = 1.05, lwd = .5,
    bg = pal[pmin(dep[-1] + 1L, length(pal))], col = "white"
  )
  points(xy$x[1], xy$y[1],
    pch = 21, cex = 1.05, lwd = .5,
    bg = "grey75", col = "white"
  )
}

out <- file.path(PKG, "figure/fig_sim_negcontrol_truths.pdf")
pdf(out, width = 12.0, height = 7.1)
par(mfrow = c(2, 3), mar = c(1.6, 1.2, 2.8, 1.2), xpd = NA)

panel_tree(balanced_tree(4L, 2L), "(A) balanced binary tree, depth 4")
panel_tree(balanced_tree(2L, 5L), "(B) balanced bushy tree, depth 2")
panel_tree(unbalanced_5_23(), "(C) unbalanced shallow tree\n5/23 split")
panel_tree(unbalanced_mixed(), "(D) unbalanced mixed-depth tree\nterminal depths 2/3/4/6")

xy <- matrix(runif(2 * M), ncol = 2)
plot(NA,
  xlim = c(-.06, 1.06), ylim = c(-.06, 1.06), asp = 1, axes = FALSE,
  xlab = "", ylab = "", main = "(E) flat plane", cex.main = 1.28, font.main = 1
)
rect(0, 0, 1, 1, border = "grey55", lwd = .9)
points(xy, pch = 21, cex = 1.15, lwd = .5, bg = "grey55", col = "white")

u <- sort(runif(M))
plot(NA,
  xlim = c(-.06, 1.06), ylim = c(-.35, .35), axes = FALSE,
  xlab = "", ylab = "", main = "(F) one-dimensional continuum",
  cex.main = 1.28, font.main = 1
)
segments(0, 0, 1, 0, col = "grey55", lwd = .9)
points(u, rep(0, M),
  pch = 21, cex = 1.15, lwd = .5,
  bg = colorRampPalette(c("#9ecae1", "#08306b"))(M), col = "white"
)
text(0, -0.20, "easiest", cex = .85, col = "grey30")
text(1, -0.20, "hardest", cex = .85, col = "grey30")
arrows(0.08, -0.13, 0.92, -0.13, length = .07, col = "grey45", lwd = 1)
dev.off()

cat("wrote ", out, "\n", sep = "")
