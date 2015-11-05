#                          ------------------------------------
#------------------------- Derived Stats Calculation functions -----------------------

#' @title Calculates perp dist detection probability p(x).
#'
#' @description
#' Calculates perp dist detection probability p(x).
#'
#' @param x perpendicular distances at which to evaluate function.
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param cov covariate matrix with each row corresponding to an observation. Must contain columns 
#' with variables appearing in \code{models}, and named accordingly, as well as column of perpendicular
#' distance, named "x". (Perpendicular distances only used if \code{to.x} is TRUE.)
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars availability hmm parameters, as for \code{\link{est.hmltm}}. Must have elements
#' \code{$Et} and \code{$Sigma.Et}
#' @param type if "link", parameter vector \code{pars} is assumed to be on the link scale, else on 
#' the natural scale
#' @param ally If TRUE calculates detection probability at all forward distances, else at zero.
#' 
#' @export
hmltm.px <- function(x,pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,
                  type="response",ally=TRUE)
{
  ymax <- survey.pars$ymax
  dy <- survey.pars$dy
  pm <- hmm.pars$pm
  Pi <- hmm.pars$Pi
  delta <- hmm.pars$delta
  
  b <- pars
  if(type!="link") 
    b <- n2w_rcpp(pars,hfun)
  
  nx <- length(x)
  n <- 1
  covb <- b
  
  if(!is.null(cov) & !(is.nullmodel(models))) { # only use covars if have them and model uses them
    n <- dim(cov)[1]
    covb <- make.covb(b,hfun,models,cov) # put covariates into paramerters
  } 
  
  nb <- length(covb)/n
  px <- matrix(rep(0,nx*n),nrow=n)
  
  for(i in 1:n) {
    start <- (i-1)*nb+1
    bi <- c(rep(covb[start:(start+nb-1)],nx)) # nx replicates of covb for ith detection
    px[i,] <- p.xy(x=x,y=rep(0,nx),hfun=hfun,b=bi,pm=pm,Pi=Pi,delta=delta,ymax=ymax,dy=dy,ally=ally)
  }
  
  return(px)
}


#' @title Calculates a bunch of derived statistics from model.
#'
#' @description
#' Calculates a one of a variety of derived statistics from model (see below).
# 
#' @param stat the statistic name (character variable).
#' @param hmmlt the model (as ouput by \code{\link{est.hmltm}}.
#' @param obs indices of rows of \code{hmmlt$xy} (i.e. which observations) to use in calculating 
#' the statistics.
#' 
#' @details
#' The following are the options for argument \code{stat}:
#' \describe{
#' \item{esw:}{effective strip width estimate.}
#' \item{invesw:}{inverse effective strip width estimate.}
#' \item{p0:}{estimated probability of detection at perpendicular distance zero.}
#' \item{p:}{estimated mean probability of detection.}
#' \item{invp:}{estimated inverse mean probability of detection.}
#' }
#' 
#' @export
calc.derived <- function(stat,hmmlt,obs=1:dim(hmmlt$xy)[1])
{
  if(stat=="esw") 
    return(fitted_esw(hmmlt,obs))
  else if(stat=="invesw")
    return(fitted_invesw(hmmlt,obs))
  else if(stat=="p0") 
    return(fitted_px(hmmlt,obs,at.x=0))
  else if(stat=="p") 
    return(fitted_p(hmmlt,obs))
  else if(stat=="invp") 
    return(fitted_invp(hmmlt,obs))
  else 
    stop(paste(stat," is an invalid stat type"))
}


#' @title Calculates fitted values for p(x) from model.
#'
#' @description
#' Calculates fitted values, p(x) for given observations, from model (optionally at given x-value).
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param at.x values at which to evaluate p(x) - see below.
#' 
#' @details 
#' If \code{hmmlt$models} is NULL and
#' \describe{
#'   \item{at.x is NULL}{returns vector of values of p(x=hmmlt$xy$x[obs]) if obs given, or vector 
#'   of values of p(x=hmmlt$xy$x) if obs not given;}
#'   \item{at.x is specified}{returns vector of values of p(x=at.x).}
#' }
#' If \code{hmmlt$models} is not NULL and
#' \describe{
#'   \item{at.x is NULL}{returns vector of values of p_obs[i](x=hmmlt$xy$x[j]) if obs given, or 
#'   vector of values of p_i(x=hmmlt$xy$x[j]) for all i if obs not given;}
#'   \item{at.x is specified}{as when at.x is NULL, but with x=at.x instead of hmmlt$xy$x.}
#'}
#'
#'@export
fitted_px <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],at.x=NULL)
{
  cov <- hmmlt$xy
  
  if(max(obs)>dim(hmmlt$xy)[1]) 
    stop("obs greater than number observations in hmmlt$xy")
  if(min(obs)<1) 
    stop("obs < 1")
  
  cov <- cov[obs,]
  
  if(!is.null(at.x)) {
    if(length(at.x)!=1 & length(at.x)!=dim(cov)[1]) 
      stop("Length of at.x inconsistent with covariate data frame.")
    
    if(length(at.x)==1) {
      cov$x <- rep(at.x,length(cov$x))
      at.x <- cov$x
    }
    else
      cov$x <- at.x
  } else {
    at.x <- cov$x
  }
  
  if(is.nullmodel(hmmlt$models))
    cov <- cov[1,]

  pars <- hmmlt$fit$par
  hfun <- hmmlt$h.fun
  models <- hmmlt$models
  survey.pars <- hmmlt$fitpars$survey.pars
  hmm.pars <- hmmlt$fitpars$hmm.pars
  px <- hmltm.px(at.x,pars,hfun,models,cov,survey.pars,hmm.pars,ally=TRUE)
  
  if(!is.nullmodel(hmmlt$models))
    px <- diag(px)
  
    return(px)
}

#' @title Calculates E[p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of fitted values, p(x) for given observations, 
#' from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
#' 
#' @export
fitted_p <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100)
{
  W <- hmmlt$fitpars$survey.pars$W
  esw <- fitted_esw(hmmlt,obs,nx)
  return(esw/W)
}

#' @title Calculates E[1/p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of inverse of fitted values, 1/p(x) for given 
#' observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
#' 
#' @export
fitted_invp <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100)
{
  W <- hmmlt$fitpars$survey.pars$W
  esw <- fitted_esw(hmmlt,obs,nx)
  return(W/esw)
}

#' @title Calculates 1/esw from model.
#'
#' @description
#' Calculates inverse of effective strip half-width, 1/(W*E[p(x)]) (where W is actual half-width) for 
#' given observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
#' 
#' @export
fitted_invesw <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100)
{
  return(1/fitted_esw(hmmlt,obs,nx))
}

#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip half-width, W*E[p(x)], where W is actual half-width, from model.
# 
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x values to use to implement Simpson's rule in perp dist dimension;
#' @param to.x If TRUE integrates only to observed x, else integrates to W
#' @param all If TRUE then returns esw for every observation, else returns only that for first obs 
#' if there are no covariates; always returns esw for every observation if there are covariates.
#' 
#' @details 
#' Calls \code{\link{hmltm.esw}} to calclate effective stript width (esw) for fitted object \code{hmmlt}.
#' 
#' @export
fitted_esw <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100,to.x=FALSE,all=FALSE)
{
  if(!is.null(obs)){
    if(max(obs)>dim(hmmlt$xy)[1]) 
      stop("obs greater than number observations in hmmlt$xy")
    if(min(obs)<1) 
      stop("obs < 1")
    
    cov <- hmmlt$xy[obs,]
  } 
  if(is.nullmodel(hmmlt$models) & !all)
    cov <- hmmlt$xy[1,]
  else
    cov <- hmmlt$xy[obs,]
  
  pars <- hmmlt$fit$par
  hfun <- hmmlt$h.fun
  models <- hmmlt$models
  survey.pars <- hmmlt$fitpars$survey.pars
  hmm.pars <- hmmlt$fitpars$hmm.pars
  
  esw <- hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type="response",to.x)
  return(esw)
}


#' @title Calculates various statistics from model.
#'
#' @description
#' Calculates various statistics from model.
# 
#' @param stat name of statistic to calculate. Valid statistics are "p0" for estimated probability 
#' at perpendicular distance zero, "p" for mean estimated detection probability over all perpendicular
#' distances, "invp" for the inverse of mean estimated detection probability over all perpendicular
#' distances, "esw" for estimated effective strip width, and "invesw" for estimated inverse of effective 
#' strip width.
#' @param b detection hazard parameter vector.
#' @param hfun detection hazard name (character).
#' @param models covariate models (see \code{\link{est.hmltm}} for details).
#' @param cov covariate values.
#' @param survey.pars survey parameter specification  (see \code{\link{est.hmltm}} for details).
#' @param hmm.pars hidden Markov model parameter specification  (see \code{\link{est.hmltm}} for 
#' details).
#' @param nx number of points at which to evaluate detection function in perpendicular distance
#' dimension.
#' @param type if "link", assumes that parameter vector \code{b} is on link scale, else assumes 
#' it is on natural scale.
#' 
#' @export
hmltm.stat <- function(stat,b,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                    type="link")
{
  if(type!="link") 
    b <- n2w_rcpp(b,hfun)
  
  if(stat=="p0")
    return(hmltm.px(x=0,b,hfun,models,cov,survey.pars,hmm.pars,type))
  else if(stat=="p") 
    return(hmltm.p(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))
  else if(stat=="invp") 
    return(1/hmltm.p(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))
  else if(stat=="esw") 
    return(hmltm.esw(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))
  else if(stat=="invesw") 
    return(1/hmltm.esw(b,hfun,models,cov,survey.pars,hmm.pars,nx,type))
  else 
    stop("Invalid stat")
}


#' @title Calculates E[p(x)] from model.
#'
#' @description
#' Calculates expected p(x) from model.
# 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' 
#' @export
hmltm.p <- function(pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                 type="response")
{
  esw <- hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type=type)
  return(esw/survey.pars$W)
}


#' @title Calculates 1/esw from model.
#'
#' @description
#' Calculates inverse of effective strip half-width from model.
# 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' 
#' @export
hmltm.invp <- function(pars,hfun,models=list(y=NULL,x=NULL),cov=NULL,survey.pars,hmm.pars,nx=100,
                    type="response")
{
  esw <- hmltm.esw(pars,hfun,models,cov,survey.pars,hmm.pars,nx,type=type)
  return(survey.pars$W/esw)
}

#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip half-width from model.
#'
#' @param pars starting parameter values, as for \code{\link{est.hmltm}}.
#' @param hfun detection hazard function name; same as argument \code{FUN} of \code{\link{est.hmltm}}.
#' @param models detection hazard covariate models, as for \code{\link{est.hmltm}}.
#' @param cov covariate matrix with each row corresponding to an observation. Must contain columns 
#' with variables appearing in \code{models}, and named accordingly, as well as column of perpendicular
#' distance, named "x". (Perpendicular distances only used if \code{to.x} is TRUE.)
#' @param survey.pars survey parameters, as for \code{\link{est.hmltm}}.
#' @param hmm.pars availability hmm parameters, as for \code{\link{est.hmltm}}. Must have elements
#' \code{$Et} and \code{$Sigma.Et}
#' @param nx number of x-values (perpendicular distances) to use in evaluating esw.
#' @param type if "link", parameter vector \code{pars} is assumed to be on the link scale, else on 
#' the natural scale
#' @param to.x if TRUE integrates only out to \code{cov$x[i]} for observation i (else integrates to 
#' \code{survey.pars$W}).
#' 
#' @details
#' Returns effective strip half-width (esw) for fitted object hmmlt, integrating 
#' using Simpson's rule.
#'
#' @export
hmltm.esw <- function(pars,hfun,models,cov,survey.pars,hmm.pars,nx=100,type="response",to.x=FALSE)
{
  n <- dim(cov)[1]
  
  if(to.x) 
    maxx <- cov$x
  else 
    maxx <- rep(survey.pars$W,n)
  
  p <- rep(0,n)
  for(i in 1:n) {
    if(maxx[i]>0){
      xs <- seq(0,maxx[i],length=nx) # set of poiints on which to evaluate p(see|x)
      px <- hmltm.px(xs,pars,hfun,models,cov[i,],survey.pars,hmm.pars,type)
      p[i] <- sintegral(px,xs)
    }
  }
  
  return(p)
}
