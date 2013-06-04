\name{fitted.esw}
\alias{fitted.esw}
\title{Calculates esw from model.}
\usage{
  fitted.esw(hmmlt, obs = 1:dim(hmmlt$xy)[1], nx = 100,
    to.x = FALSE, all = FALSE)
}
\arguments{
  \item{hmmlt}{output from \code{\link{fit.hmltm}}}

  \item{obs}{observations (row numbers of \code{hmmlt$xy})
  for which to calculate esw}

  \item{b}{parameters (on link scale) at which to evaluate
  mu. If NULL, uses parameter estimates in hmmlt}

  \item{nx}{number of x values to use to implement
  Simpson's rule in perp dist dimension;}

  \item{to.x}{If TRUE integrates only to observed x, else
  integrates to W}

  \item{all}{If TRUE then returns esw for every
  observation, else returns only that for first obs if
  there are no covariates; always returns esw for every
  observation if there are covariates.}
}
\description{
  Calculates esw from model.
}
\details{
  Calls \code{\link{hmltm.esw}} to calclate effective
  stript width (esw) for fitted object \code{hmmlt}.
}

