#                          --------------------------
#------------------------- Goodness of Fit functions -----------------------


#' @title Goodness-of-fit in forward dimension.
#'
#' @description
#' Calculates goodness-of-fit in forward dimension, plots fit, and returns p-value and other stuff.
#' Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs), and Cramer-von Mises
#' p-value (which is based on all cdf values).
# 
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param ks.plot If TRUE, does Q-Q plot. Point corresponding to largest difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
#' @param seplots if TRUE does additional diagnostic plots
#' @param smult multiplier to size circles in third plot.
#' @param ymax forward distance at which detection probability is assumed to be zero. 
#' 
#' @export
hmmlt.gof.y <- function(hmltm,ks.plot=TRUE,seplots=FALSE,smult=5,ymax=hmmlt$fitpars$survey.pars$ymax)
{
  hmmlt <- hmltm$hmltm.fit
  dat <- hmmlt$xy[!is.na(hmmlt$xy$y),]
  hfun <- hmmlt$h.fun
  b <- hmmlt$fit$b
  
  if(is.null(dat$y)) 
    stop("Can't do goodness-of-fit in y-dimension when don't have y observations!")
  
  models <- hmmlt$models
  dy <- hmmlt$fitpars$survey.pars$dy
  Pi <- hmmlt$fitpars$hmm.pars$Pi
  pm <- hmmlt$fitpars$hmm.pars$pm
  delta <- hmmlt$fitpars$hmm.pars$delta
  n <- length(dat$x)
  covb <- makeCovPar(b,hfun,models,dat) # put covariates into paramerters
  Fy <- p.xy(dat$x,dat$y,hfun,b=covb,pm,Pi,delta,ymax,dy,ally=FALSE,cdf=TRUE)
  F0 <- p.xy(dat$x,dat$y,hfun,b=covb,pm,Pi,delta,ymax,dy,ally=TRUE,cdf=FALSE)
  #  Fy=Fx.cox.Q(dat$x,dat$y,mu,ystart,Q,b) # area up to y
  #  F0=Fx.cox.Q(dat$x,rep(0,n),mu,ystart,Q,b) # area up to y=0
  Fy0 <- Fy/F0
  Fy0.order <- order(Fy0)
  yy <- dat$y[Fy0.order]
  xx <- dat$x[Fy0.order]
  cdf <- Fy0[Fy0.order]
  e.cdf <- order(cdf)/n
  
  # K-S statistic
  dF <- cdf-e.cdf
  worst <- which(abs(dF)==max(abs(dF))) # mark point on which Kolmogarov test hinges
  Dn <- max(abs(dF))*sqrt(n)
  p <- p.kolomogarov(Dn)
  p.cvm <- cvm.test(Fy0)$p.value
  
  # plots
  if(ks.plot) {
    plot(1-e.cdf,cdf,xlab="Empirical Distribution Function",ylab="Cumulative Distribution Function",
         main="Forward Dist. Q-Q Plot",xlim=c(0,1),ylim=c(0,1),pch="+")
    lines(c(0,1),c(1,0))
    points(1-e.cdf[worst],cdf[worst],col="red") # mark point on which Kolmogarov test hinges
    if(seplots){
      plot(xx,dF,xlab="Perpendicular distance",ylab="CDF-Empirical CDF")
      lines(c(0,max(xx)),c(0,0),lty=2)
      size <- smult*dF/max(dF)
      dFcol <- c("red","black","black")
      plot(dat$x[Fy0.order],dat$y[Fy0.order],xlab="Perpendicular distance",ylab="Forward distance",
           cex=abs(size),col=dFcol[sign(dF)+2],main="CDF-Empirical CDF (red=negative)")
    }
  }
  
  return(list(p.ks=1-p,p.cvm=p.cvm,qq.x=e.cdf,qq.y=cdf,y=yy))
}


#' @title Kolmogarov-Smirnov goodness-of-fit p-value.
#'
#' @description
#' Kolmogarov-Smirnov goodness-of-fit p-value calculation.
# 
#' @param x value of Kolmogarov-distributed random variable at which to evaluate.
#' @param inf approximation to infinity (a large number).
#' @param dp approximation convergence criterion.
#' 
#' @details
#' Calculates p-value for Kolmogarov distribution at x, approximating infite sum
#' \code{sqrt(2*pi)/x*sum{i=1}^infty exp(-(2i-1)^2*pi^2/(8x^2)))}
#' by a finite sum to inf (default 1000) if sum to inf and inf+1 differ by less
#' than dp (default 1e-4), else sum until difference is less than dp.
#' 
#' @export
p.kolomogarov <- function(x,inf=1000,dp=1e-4)
{
  infsum <- rep(0,inf)
  i <- 1:inf
  K <- sqrt(2*pi)/x
  p <- p1 <- K*sum(exp(-(2*i-1)^2*pi^2/(8*x^2)))
  dp <- 1
  
  while(dp>1e-4) {
    inf <- inf+1
    p <- p1+K*exp(-(2*inf-1)^2*pi^2/(8*x^2))
    dp <- abs(p-p1)
  }
  
  return(p=p)
}


#' @title Goodness-of-fit in perpendicular dimension.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular dimension, plots fit, and returns p-value and 
#' other stuff. Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs), and Cramer-von Mises
#' p-value (which is based on all cdf values).
# 
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param ks.plot If TRUE, does Q-Q plot. Point corresponding to largest difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
#' 
#' @importFrom goftest cvm.test
#' 
#' @export
hmmlt.gof.x <- function(hmltm,ks.plot=TRUE)
{
  hmmlt <- hmltm$hmltm.fit
  n <- length(hmmlt$xy$x)
  edf <- (1:n)/n
  cdf <- fitted_esw(hmmlt,to.x=TRUE,all=TRUE)/fitted_esw(hmmlt,all=TRUE)
  cdf.order <- order(cdf)
  cdf <- cdf[cdf.order]
  e.cdf <- cdf.order/n
  
  # K-S statistic
  dF <- cdf-edf
  worst <- which(abs(dF)==max(abs(dF)))
  Dn <- max(abs(dF))*sqrt(n)
  p.ks <- p.kolomogarov(Dn)
  p.cvm <- cvm.test(cdf)$p.value # Under model, cdf values are from uniform; default for cvm.test is "punif"
  
  if(ks.plot) {
    plot(edf,cdf,pch="+",xlim=c(0,1),ylim=c(0,1),xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",main="Perp. Dist. Q-Q Plot")
    lines(c(0,1),c(0,1))
    points(edf[worst],cdf[worst],col="red")
  }
  
  return(list(p.ks=1-p.ks,p.cvm=p.cvm,qq.x=edf,qq.y=cdf,x=hmmlt$xy$x[cdf.order]))
}
