\name{hmltm.esw}
\alias{hmltm.esw}
\title{Calculates esw from model.}
\usage{
  hmltm.esw(pars, hfun, models, cov, survey.pars, hmm.pars,
    nx = 100, type = "response", to.x = FALSE)
}
\arguments{
  \item{pars.}{}

  \item{hfun.}{}

  \item{models.}{}

  \item{cov.}{}

  \item{survey.pars.}{}

  \item{hmm.pars.}{}

  \item{nx.}{}

  \item{type.}{}

  \item{to.x.}{}
}
\description{
  Calculates esw from model.
}
\details{
  Returns effective stript width (esw) for fitted object
  hmmlt, integrating using Simpson's rule.

  Input: hmmlt: output from fit.hmltm obs : observatoins
  (row numbers of hmmlt$xy) for which to calculate esw b :
  parameters (on link scale) at which to evaluate mu. If
  NULL, uses parameter estimates in hmmlt nx : number of x
  values to use to implement Simpson's rule; to.x : If TRUE
  integrates only to observed x, else integrates to W type
  : If "link", interprets pars as being on link scale, else
  natural scale
}

