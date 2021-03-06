% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsltm.R
\name{h.plot}
\alias{h.plot}
\title{Detection hazard plotting.}
\usage{
h.plot(hfun, pars, dat = NULL, models = NULL, xrange = c(0, 50),
  yrange = xrange, nx = 50, ny = nx, type = "contour",
  nlevels = 20, add = FALSE, col = "black", logscale = FALSE,
  xlab = "Perpendicular distance", ylab = "Forward distance",
  theta = 90, phi = 35, ...)
}
\arguments{
\item{hfun}{detection hazard function name}

\item{pars}{detection hazard function parameter vector}

\item{dat}{data frame.}

\item{models}{model list, as for \code{\link{est.hmltm}} for example.}

\item{xrange}{range of x-axis.}

\item{yrange}{range of y-axis.}

\item{nx}{number of points on x-axis at which to evaluate detection hazard function.}

\item{ny}{number of points on y-axis at which to evaluate detection hazard function.}

\item{type}{"contour", "persp", "image" or "both" (for image and contour).}

\item{nlevels}{number of contour levels.}

\item{add}{if TRUE adds to existing plot, else creates new plot.}

\item{col}{colour of plot.}

\item{logscale}{If TRUE, plots hazard values on log scale.}

\item{xlab}{x label.}

\item{ylab}{y label.}

\item{theta}{argument for \code{\link{image}}.}

\item{phi}{argument for \code{\link{image}}.}

\item{...}{other arguments to image, contour or persp.}
}
\description{
Plots detection hazard contours.
}
