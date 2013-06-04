\name{fitted.invesw}
\alias{fitted.invesw}
\title{Calculates E[1/esw] from model.}
\usage{
  fitted.invesw(hmmlt, obs = 1:dim(hmmlt$xy)[1], nx = 100)
}
\arguments{
  \item{hmmlt}{output from fit.hmltm()}

  \item{obs}{observations (row numbers of \code{hmmlt$xy})
  for which to calculate esw}

  \item{nx}{number of x values to use to implement
  Simpson's rule in perp dist dimension;}
}
\description{
  Calculates E[1/esw] from model.
}

