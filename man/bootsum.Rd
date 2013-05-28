\name{bootsum}
\alias{bootsum}
\title{Summarise bootstrap results.}
\usage{
  bootsum(bests, cilevel = 0.95)
}
\arguments{
  \item{bests}{: output from bs.hmltm().}

  \item{cilevel:}{confidence level for confidence intervals
  by percentile method}
}
\value{
  Returns a list with elements as follows for each
  statistic in \code{bests}: \itemize{ \item{nbad} {number
  of bad estimates of the statistic in question that were
  excluded from the summary. (Bad estimates occur, for
  example, when a bootstrap sample involves no detections
  and the estimation function tries to calculate mean group
  size.)} \item{means} {Bootstrap means of the statistic in
  question.} \item{cv} {Percentage CVs of the statistic in
  question.} \item{se} {SEs of the statistic in question.}
  \item{lower} {Lower \code{cilevel} percentiles of the
  statistic in question.} \item{upper} {Upper
  \code{cilevel} percentiles of the statistic in question.}
  }
}
\description{
  Produces summarry of bootstrap results, in a list
}

