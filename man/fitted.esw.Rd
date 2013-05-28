\name{fitted.esw}
\alias{fitted.esw}
\title{Calculates esw from model.}
\usage{
  fitted.esw(hmmlt, obs = 1:dim(hmmlt$xy)[1], nx = 100,
    to.x = FALSE, all = FALSE)
}
\arguments{
  \item{hmmlt.}{}

  \item{obs.}{}

  \item{nx.}{}

  \item{to.x.}{}

  \item{all.}{}
}
\description{
  Calculates esw from model.
}
\details{
  Calls hmltm.esw to calclate effective stript width (esw)
  for fitted object hmmlt. Input: hmmlt: output from
  fit.hmltm() obs : observatoins (row numbers of hmmlt$xy)
  for which to calculate esw b : parameters (on link scale)
  at which to evaluate mu. If NULL, uses parameter
  estimates in hmmlt nx : number of x values to use to
  implement Simpson's rule in perp dist dimension; to.x :
  If TRUE integrates only to observed x, else integrates to
  W all : If TRUE then returns esw for every observation,
  else returns only that for first obs if there are no
  covariates; always returns esw for every observation if
  there are covariates.
}

