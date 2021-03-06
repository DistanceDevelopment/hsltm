% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hsltm.R
\name{instant.a}
\alias{instant.a}
\title{Instantaneous availability correction factor calculation for multiple availability models.}
\usage{
instant.a(hmm.pars, Et = NULL)
}
\arguments{
\item{hmm.pars}{is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
in element $Pi (where m is number of availability parameter sets).}

\item{Et}{is a 2xm matrix with expect times Unavailable (row 1) and Available (row 2).}
}
\description{
Calculates the probability that an animal is available at an instant, for each of a set of m 
2-state Markov model availability parameters passed to it and returns probabilities and their mean.
}
\details{
If Et is given, it is used and hmm.pars is ignored.
}
\examples{
Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3);pm=NULL
hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
instant.a(hmm.pars)
instant.a(NULL,Et=matrix(c(Eu,Ea),ncol=2,byrow=TRUE))
instant.a(NULL,c(Eu[1],Ea[1]))
}
\seealso{
\code{\link{mclaren.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
}
