\name{fitted.invp1}
\alias{fitted.invp1}
\title{Calculates E[1/p(x)] from model.}
\usage{
  fitted.invp1(hmmlt, obs = 1:dim(hmmlt$xy)[1], nx = 100,
    W = NULL)
}
\arguments{
  \item{hmmlt.}{}

  \item{obs.}{}

  \item{nx.}{}

  \item{W.}{}
}
\description{
  Calculates E[1/p(x)] from model.
}
\details{
  Identical to fitted.invp but returns data frame instead
  of numerical scalar or vector. This is to allow it to be
  used in NDest for estimating density and abundance. (Also
  has extra parameter: W)
}

