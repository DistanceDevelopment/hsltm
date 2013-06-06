\name{p.xy}
\alias{p.xy}
\title{Calculates probability of first observing animal at (x,y).}
\usage{
  p.xy(x, y, hfun, b, pm, Pi, delta, ymax, dy, theta.f,
    theta.b, ally = FALSE, cdf = FALSE)
}
\arguments{
  \item{x}{perpendicular distance.}

  \item{y}{forward distance.}

  \item{hfun}{detection hazard function name.}

  \item{b}{detection hazard function parameters.}

  \item{pm}{state-dependent Bernoulli parameters.}

  \item{Pi}{Markov model transition probability matrix.}

  \item{delta}{Markov model stationary distribution.}

  \item{ymax}{forward distance at or beyond which things
  are undetectable.}

  \item{dy}{distance step in forward distance direction
  corresponding to the time step size of the Markov model
  governed by \code{Pi}.}

  \item{theta.f}{REDUNDANT must = 0.}

  \item{theta.b}{REDUNDANT must = 90.}

  \item{ally}{If TRUE calculates detection probability at
  all forward distances, else at zero.}

  \item{cdf}{If TRUE returns cdf, not detection
  probability.}
}
\description{
  Just calls p.xy1 >= once, cycling through 3rd index of Pi
  and 2nd indices of pm and delta.
}

