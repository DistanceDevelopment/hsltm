\name{hmmlt.gof.y}
\alias{hmmlt.gof.y}
\title{Goodness-of-fit in forward dimension.}
\usage{
  hmmlt.gof.y(hmltm, ks.plot = TRUE, seplots = FALSE,
    smult = 5)
}
\arguments{
  \item{hmltm}{fitted model, as output by
  \code{\link{est.hmltm}}}

  \item{ks.plot}{If TRUE, does Q-Q plot.}

  \item{seplots}{if TRUE does additional diagnostic plots}

  \item{smult}{multiplier to size circles in third plot.}
}
\description{
  Calculates goodness-of-fit in forward dimension, plots
  fit, and returns p-value and other stuff.
}

