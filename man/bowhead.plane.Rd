\docType{data}
\name{bowhead.plane}
\alias{bowhead.plane}
\title{Bowhead whale aerial survey dataset from Greenland.}
\format{An mrds data frame (with 2 rows per detection) with 86 observations on the following 11
variables.
 \describe{
   \item{\code{stratum}:}{ stratum number (a numeric vector).}
   \item{\code{area}:}{ stratum surface area (a numeric vector).}
   \item{\code{transect}:}{ transect number (a numeric vector).}
   \item{\code{L}:}{ transect length (a numeric vector).}
   \item{\code{size}:}{ size of each detected group (NA if no detection).}
   \item{\code{object}:}{ unique identifier for each detection (a numeric vector; NA if no detection).}
   \item{\code{side}:}{ the side of the plane from which each detection was made (NA if no detection).}
   \item{\code{obs}:}{ the observer who made each detection (NA if no detection).}
   \item{\code{bf}:}{ Beaufort sea state at the time of each detection (NA if no detection).}
   \item{\code{x}:}{ perpendicular distances to detections (NA if no detection).}
   \item{\code{y}:}{ perpendicular distances to detections (NA if no detection).}
 }}
\source{
  Mads Peter Heide-Jorgessen.
}
\description{
  Data from the 2012/13 aerial survey of bowhead whales off
  Greenland. These data were used together with the dataset
  \code{\link{bowhead.hmm.pars}} in the analyses of (**
  reference to Mads Peter's paper **).
}
\examples{
data(bowhead.plane)
}
\references{
  None as yet
}
\keyword{datasets}

