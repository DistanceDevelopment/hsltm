% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsltm.R
\name{hmmlt.gof.y}
\alias{hmmlt.gof.y}
\title{Goodness-of-fit in forward dimension.}
\usage{
hmmlt.gof.y(hmltm, ks.plot = TRUE, seplots = FALSE, smult = 5,
  ymax = hmmlt$fitpars$survey.pars$ymax)
}
\arguments{
\item{hmltm}{fitted model, as output by \code{\link{est.hmltm}}}

\item{ks.plot}{If TRUE, does Q-Q plot. Point corresponding to largest difference between
empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.}

\item{seplots}{if TRUE does additional diagnostic plots}

\item{smult}{multiplier to size circles in third plot.}

\item{ymax}{forward distance at which detection probability is assumed to be zero.}
}
\description{
Calculates goodness-of-fit in forward dimension, plots fit, and returns p-value and other stuff.
Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
based on only the largest difference between emprical and theoretical cdfs), and Cramer-von Mises
p-value (which is based on all cdf values).
}
