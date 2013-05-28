\name{bs.hmltm}
\alias{bs.hmltm}
\title{Bootstrap for hmltm model.}
\usage{
  bs.hmltm(hmltm.est, B, hmm.pars.bs = NULL, bs.trace = 0,
    report.by = 10, fixed.avail = FALSE)
}
\arguments{
  \item{\code{hmltm.est}}{output from
  \code{\link{est.hmltm}}.}

  \item{B}{number of bootstrap replicates to do.}

  \item{hmm.pars.bs}{output \code{\link{hmmpars.boot}},
  containing sets of refitted HMM parameter values.}

  \item{bs.trace}{amount of reporting that
  \code{\link{optim}} should do while fitting models to
  bootstrapped datasets. (\code{bs.trace}=0 is no
  reporting, which is fasters. See component \code{trace}
  of parameter \code{control.opt} of function
  \code{\link{optim}} for more details.)}

  \item{report.by}{frequency with which to report count of
  number of bootstraps completed.}

  \item{fixed.avail}{whether to treat the hidden Markov
  availability model as fixed (\code{fixed.avail}=TRUE) or
  random (\code{fixed.avail}=FALSE). If
  \code{fixed.avail}=FALSE, the availability model is also
  bootstrapped (see Details below).}
}
\description{
  Stratified nonparameteric bootstrap of hidden Markov line
  transect model (hmltm) with transects as the sampling
  units.
}
\details{
  If
}
\seealso{
  \code{\link{bootsum}} summarises output from this
  function.
}

