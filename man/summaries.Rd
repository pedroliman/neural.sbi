% Hand-written (mirrors roxygen tags); regenerate with devtools::document()
\name{summaries}
\alias{summaries}
\alias{summary.nsbi_samples}
\alias{summary.nsbi_posterior}
\alias{summary.nsbi_npe}
\alias{as.data.frame.nsbi_samples}
\title{Summaries and tidy accessors}
\usage{
\method{summary}{nsbi_samples}(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)

\method{summary}{nsbi_posterior}(object, n = 1000L, x = NULL, ...)

\method{summary}{nsbi_npe}(object, ...)

\method{as.data.frame}{nsbi_samples}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{object}{An \code{nsbi_samples} matrix, \code{nsbi_posterior}, or
\code{nsbi_npe} fit.}

\item{probs}{Quantiles to report.}

\item{n}{Number of draws used to summarize a posterior.}

\item{x}{Observation to condition on (defaults to the posterior's
\code{x_obs}); for \code{as.data.frame}, the samples object.}

\item{row.names, optional, ...}{Standard \code{as.data.frame}/method
arguments.}
}
\value{
For samples/posteriors, a data frame with one row per parameter (mean, sd,
quantiles). For fits, an invisible list of training metadata.
}
\description{
\code{summary()} methods for fits, posteriors, and samples, plus
\code{as.data.frame()} for posterior draws so results drop straight into
data-frame workflows.
}
