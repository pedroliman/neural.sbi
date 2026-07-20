#' Summaries and tidy accessors
#'
#' `summary()` methods for fits, posteriors, and samples, plus
#' `as.data.frame()` for posterior draws so results drop straight into
#' data-frame workflows (dplyr, ggplot2, ...).
#'
#' @name summaries
NULL

#' @export
as.data.frame.nsbi_samples <- function(x, row.names = NULL, optional = FALSE, ...) {
  m <- unclass(x)
  attr(m, "acceptance_rate") <- NULL
  colnames(m) <- colnames(m) %||% paste0("theta", seq_len(ncol(m)))
  as.data.frame(m, row.names = row.names, optional = optional, ...)
}

#' Summarize posterior samples
#'
#' @param object An `nsbi_samples` matrix from [sample()].
#' @param probs Quantiles to report.
#' @param ... Unused.
#' @return A data frame with one row per parameter: mean, sd, and quantiles.
#' @export
summary.nsbi_samples <- function(object,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                 ...) {
  m <- unclass(object)
  q <- t(apply(m, 2, stats::quantile, probs = probs))
  colnames(q) <- paste0("q", 100 * probs)
  out <- data.frame(
    parameter = colnames(m) %||% paste0("theta", seq_len(ncol(m))),
    mean = colMeans(m),
    sd = apply(m, 2, stats::sd),
    q, row.names = NULL, check.names = FALSE
  )
  out
}

#' Summarize a posterior by drawing samples
#'
#' @param object An `nsbi_posterior`.
#' @param n Number of draws used for the summary.
#' @param x Observation to condition on (defaults to the posterior's `x_obs`).
#' @param ... Passed to [summary.nsbi_samples()].
#' @return A data frame as in [summary.nsbi_samples()].
#' @export
summary.nsbi_posterior <- function(object, n = 1000L, x = NULL, ...) {
  draws <- sample.nsbi_posterior(object, n = n, obs = x)
  summary(draws, ...)
}

#' @export
summary.nsbi_npe <- function(object, ...) {
  info <- list(
    density_estimator = object$density_estimator,
    dim_theta = object$dim_theta,
    dim_x = object$dim_x,
    n_simulations = object$n_simulations,
    best_val_loss = object$de$best_val_loss %||% NA_real_,
    epochs_trained = if (!is.null(object$de$history)) nrow(object$de$history)
                     else NA_integer_
  )
  print(object)
  if (!is.na(info$epochs_trained)) {
    cat(sprintf("  epochs trained    : %d\n", info$epochs_trained))
  }
  invisible(info)
}
