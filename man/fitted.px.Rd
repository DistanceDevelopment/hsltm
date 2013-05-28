\name{fitted.px}
\alias{fitted.px}
\title{Calculates fitted values for p(x) from model.}
\usage{
  fitted.px(hmmlt, obs = 1:dim(hmmlt$xy)[1], at.x = NULL,
    ny = 100)
}
\arguments{
  \item{hmmlt.}{}

  \item{obs.}{}

  \item{at.x.}{}

  \item{ny.}{}
}
\description{
  Calculates fitted values, p(x) for given observations,
  from model (optionally at given x-value).
}
\details{
  If hmmlt$models is NULL and at.x is NULL : returns vector
  of values of p(x=hmmlt$xy$x[obs]) if obs given, or vector
  of values of p(x=hmmlt$xy$x) if obs not given; at.x
  specified: returns vector of values of p(x=at.x). If
  hmmlt$models is not NULL and at.x is NULL : returns
  vector of values of p_obs[i](x=hmmlt$xy$x[j]) if obs
  given, or vector of values of p_i(x=hmmlt$xy$x[j]) for
  all i if obs not given; at.x specified: as when at.x is
  NULL, but with x=at.x instead of hmmlt$xy$x.
}

