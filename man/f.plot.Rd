\name{f.plot}
\alias{f.plot}
\title{Survivor function plotting.}
\usage{
  f.plot(hfun, pars, dat = NULL, models = NULL, pm, Pi,
    delta, ymax, dy, theta.f = 0, theta.b = 90,
    xrange = c(0, 50), yrange = xrange, nx = 50, ny = nx,
    xlab = "Perpendicular distance",
    ylab = "Forward distance", type = "contour",
    nlevels = 20, add = FALSE, col = "black",
    logscale = FALSE, theta = 90, phi = 35,
    zlab = "pdf of first detections", values = FALSE, ...)
}
\arguments{
  \item{hfun}{detection hazard function name}

  \item{dat}{data frame.}

  \item{models}{model list, as for \code{\link{est.hmltm}}
  for example.}

  \item{pm}{state-dependent Bernoulli parameter values.}

  \item{Pi}{Markov chain transition probability matrix.}

  \item{delta}{stationary distribution of \code{Pi}.}

  \item{ymax}{a forward distance beyond which things can't
  be detected.}

  \item{theta.f}{REDUNDANT must = 0.}

  \item{theta.b}{REDUNDANT must = 90.}

  \item{xrange}{range of x-axis.}

  \item{yrange}{range of y-axis.}

  \item{nx}{number of points on x-axis at which to evaluate
  detection hazard function.}

  \item{ny}{number of points on y-axis at which to evaluate
  detection hazard function.}

  \item{xlab}{x label.}

  \item{ylab}{y label.}

  \item{type}{"contour", "persp", "image" or "both" (for
  image and contour).}

  \item{nlevels}{number of contour levels.}

  \item{add}{If TRUE, adds to current plot, else makes new
  plot.}

  \item{col}{colour of plot.}

  \item{logscale}{If TRUE, plots hazard values on log
  scale.}

  \item{theta}{argument for \code{\link{image}}.}

  \item{phi}{argument for \code{\link{image}}.}

  \item{zlab}{label in z-dimension for persp plots}

  \item{values}{If TRUE, returns hazard function values on
  nx by by grid.}

  \item{...}{other arguments to image, contour or persp.}
}
\description{
  Plots contours of (1- the survivor function) with HMM
  availbility pars (Pi, pm) and detection hazard
  hfun(x,y,b,...).
}

