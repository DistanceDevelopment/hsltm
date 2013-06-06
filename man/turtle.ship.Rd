\docType{data}
\name{turtle.ship}
\alias{turtle.ship}
\title{Shipboard turtle survey dataset from US east coast.}
\format{A data frame with 48 observations on the following 12 variables.
 \describe{
   \item{\code{x}:}{ perpendicular distances to detections (NA if no detection).}
   \item{\code{y}:}{ perpendicular distances to detections (NA if no detection).}
   \item{\code{bf}:}{ Beaufort sea state at the time of each detection (NA if no detection).}
   \item{\code{we}:}{ weather code at the time of each detection (NA if no detection).}
   \item{\code{cu}:}{ cue tpe of each detection (NA if no detection).}
   \item{\code{transect}:}{ transect number (a numeric vector). In this dataset this is a dummy
   value - just there to have something in this compulsory variable.}
   \item{\code{stratum}:}{ stratum number (a numeric vector). In this dataset this is a dummy
   value - just there to have something in this compulsory variable.}
   \item{\code{area}:}{ stratum surface area (a numeric vector). In this dataset this is a dummy
   value - just there to have something in this compulsory variable.}
   \item{\code{L}:}{ transect length (a numeric vector). In this dataset this is a dummy
   value - just there to have something in this compulsory variable.}
   \item{\code{object}:}{ unique identifier for each detection (a numeric vector; NA if no detection).}
   \item{\code{size}:}{ size of each detected group (NA if no detection).}
   \item{\code{seen}:}{ REDUNDANT: 1 for all detections (NAA if no detection).}
 }}
\source{
  Virginia Aquarium.
}
\description{
  Subset of a shipboard survey dataset from the US east
  coast.
}
\details{
  Temporary dataset with which to test package.
}
\examples{
data(turtle.ship)
}
\references{
  None as yet
}
\keyword{datasets}

