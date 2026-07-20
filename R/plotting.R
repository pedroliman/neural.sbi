#' Visualize posterior samples
#'
#' A dependency-free (base graphics) pair plot: 1-D marginal densities on the
#' diagonal and 2-D scatter/contours off-diagonal, with optional markers for a
#' reference (e.g. true) parameter value. Analogous to `sbi`'s `pairplot`.
#'
#' @param samples A matrix of posterior draws (rows = draws), or an
#'   `nsbi_samples` object.
#' @param truth Optional reference parameter vector to overlay.
#' @param labels Optional parameter labels.
#' @param limits Optional list/matrix of per-parameter c(lo, hi) axis limits.
#' @param col Point colour.
#' @param ... Passed to plotting calls.
#' @return Invisibly, the samples.
#' @export
pairplot <- function(samples, truth = NULL, labels = NULL, limits = NULL,
                     col = grDevices::adjustcolor("steelblue", 0.4), ...) {
  X <- as_theta_matrix(samples)
  d <- ncol(X)
  if (is.null(labels)) {
    labels <- colnames(X) %||% paste0("theta[", seq_len(d), "]")
  }
  op <- graphics::par(mfrow = c(d, d), mar = c(2.2, 2.2, 0.6, 0.6),
                      mgp = c(1.2, 0.4, 0), tcl = -0.3)
  on.exit(graphics::par(op), add = TRUE)
  lim <- function(j) {
    if (!is.null(limits)) {
      if (is.list(limits)) return(limits[[j]])
      return(limits[j, ])
    }
    range(X[, j])
  }
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      if (i == j) {
        dens <- stats::density(X[, i])
        graphics::plot(dens, main = "", xlab = labels[i], ylab = "",
                       xlim = lim(i), col = "steelblue", ...)
        if (!is.null(truth)) graphics::abline(v = truth[i], col = "firebrick",
                                              lwd = 2)
      } else {
        graphics::plot(X[, j], X[, i], pch = 16, cex = 0.3, col = col,
                       xlab = labels[j], ylab = labels[i],
                       xlim = lim(j), ylim = lim(i), ...)
        if (!is.null(truth)) {
          graphics::points(truth[j], truth[i], col = "firebrick", pch = 4,
                           lwd = 2, cex = 1.5)
        }
      }
    }
  }
  invisible(samples)
}

#' Plot an SBC rank histogram
#'
#' Uniform bars indicate calibration; a U shape means the posterior is too
#' narrow (overconfident); an inverted-U means it is too wide.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param param Which parameter index to plot (default 1).
#' @param bins Number of histogram bins.
#' @return Invisibly, the rank vector.
#' @export
plot_sbc <- function(sbc_result, param = 1L, bins = 20L) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  r <- sbc_result$ranks[, param]
  L <- sbc_result$n_posterior_samples
  h <- graphics::hist(r, breaks = seq(0, L, length.out = bins + 1L),
                      plot = FALSE)
  expected <- sbc_result$n_sbc / bins
  graphics::plot(h, col = "grey80", border = "white",
                 main = sprintf("SBC ranks: parameter %d", param),
                 xlab = "rank of true value")
  graphics::abline(h = expected, col = "firebrick", lwd = 2, lty = 2)
  ci <- stats::qbinom(c(0.005, 0.995), sbc_result$n_sbc, 1 / bins)
  graphics::abline(h = ci, col = "firebrick", lwd = 1, lty = 3)
  invisible(r)
}

#' Plot nominal vs. empirical credible-interval coverage
#'
#' Well-calibrated posteriors lie on the diagonal. Curves above the diagonal
#' mean the posterior is too wide (conservative); below means overconfident.
#' A shaded band shows the Monte-Carlo uncertainty from the finite number of
#' SBC trials.
#'
#' @param sbc_result An `nsbi_sbc` object from [sbc()].
#' @param levels Nominal credibility levels to evaluate.
#' @return Invisibly, the coverage data frame from [expected_coverage()].
#' @export
plot_coverage <- function(sbc_result, levels = seq(0.05, 0.95, by = 0.05)) {
  stopifnot(inherits(sbc_result, "nsbi_sbc"))
  cov <- expected_coverage(sbc_result, levels = levels)
  d <- ncol(cov) - 1L
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
                 xlab = "nominal credibility level",
                 ylab = "empirical coverage",
                 main = "Expected coverage")
  # Monte-Carlo band around the diagonal (99% binomial interval)
  n <- sbc_result$n_sbc
  lo <- stats::qbinom(0.005, n, cov$nominal) / n
  hi <- stats::qbinom(0.995, n, cov$nominal) / n
  graphics::polygon(c(cov$nominal, rev(cov$nominal)), c(lo, rev(hi)),
                    col = grDevices::adjustcolor("grey60", 0.3), border = NA)
  graphics::abline(0, 1, col = "grey40", lty = 2)
  cols <- grDevices::hcl.colors(max(d, 2L), "Dark 2")
  for (j in seq_len(d)) {
    graphics::lines(cov$nominal, cov[[j + 1L]], col = cols[j], lwd = 2)
  }
  graphics::legend("topleft", legend = colnames(cov)[-1L], col = cols[seq_len(d)],
                   lwd = 2, bty = "n", cex = 0.8)
  invisible(cov)
}

#' @export
print.nsbi_samples <- function(x, ...) {
  acc <- attr(x, "acceptance_rate")
  cat(sprintf("<nsbi_samples> %d draws x %d parameters\n", nrow(x), ncol(x)))
  if (!is.null(acc)) cat(sprintf("  support acceptance rate: %.3f\n", acc))
  cls <- setdiff(class(x), "nsbi_samples")
  print(utils::head(`class<-`(unclass(x), cls)))
  invisible(x)
}
