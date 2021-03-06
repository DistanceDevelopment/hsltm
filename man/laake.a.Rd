% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsltm.R
\name{laake.a}
\alias{laake.a}
\title{Laake's availability correction factor calculation for multiple availability models.}
\usage{
laake.a(hmm.pars, ymax, spd = NULL)
}
\arguments{
\item{hmm.pars}{is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
in element $Pi (where m is number of availability parameter sets).}

\item{ymax}{is maximum forward distance to consider (`w' in Laake et al. (1997), Eqn (4)).}

\item{spd}{speed of observer. Needed to convert time (units of \code{hmm.pars}) to distance (units
of \code{ymax}).}
}
\description{
Calculates the probability that an animal is available at least once while within detectable
region, using the method of Laake et al. (1997), Eqn (4). Does this for each of a set of m 
2-state Markov model availability parameters passed to it and returns probabilities and their mean.
}
\details{
See Laake et al. (1997), Eqn (4), or Borchers et al. (2013) for details of Laake's method.
}
\examples{
Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
laake.a(hmm.pars,ymax=200,spd=4) # because hmm.pars is in units of time, not distance, you need to specify spd.
}
\references{
Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.

Laake, J., Calambokidis, J., Osmek, S., and Rugh, D. 1997. Probability of detecting harbor 
porpoise from aerial surveys: estimating g(0). Journal of Wildlife Management 61, 63-75.
}
\seealso{
\code{\link{instant.a}}, \code{\link{mclaren.a}}, \code{\link{richard.a}}
}
