#                     -------------------------------------------
#-------------------- Density and Abundance Estimation Functions --------------------

#' @title Calculates E[1/p(x)] from model.
#'
#' @description
#' Calculates mean over perpendicular distance (x) of inverse of fitted values, 1/p(x) for given 
#' observations, from model.
#'  
#' @param hmmlt output from \code{\link{fit.hmltm}}
#' @param obs observations (row numbers of \code{hmmlt$xy}) for which to calculate esw 
#' @param nx number of x-values (perpendicular distance values) to use in calculation.
#' @param W actual half-width over which to integrate (overrides W in \code{hmmlt}).
#' 
#' @details
#' Identical to fitted_invp but returns data frame instead of numerical scalar or vector.
#' This is to allow it to be used in NDest for estimating density and abundance. (Also has extra 
#' parameter: W)
#' 
#' @export
fitted_invp1 <- function(hmmlt,obs=1:dim(hmmlt$xy)[1],nx=100,W=NULL){
  if(is.null(W))
    W <- hmmlt$fitpars$survey.pars$W
  esw <- fitted_esw1(hmmlt,obs,nx,W=W)
  
  return(data.frame(stratum=esw$stratum,transect=esw$transect,object=esw$object,invp=W/esw$esw))
}


#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip width from fitted model.
# 
#' @param hmmlt output from fit.hmltm()
#' @param obs observations (row numbers of hmmlt$xy) for which to calculate esw.
#' @param nx number of x values to use to implement Simpson's rule in perp dist dimension.
#' @param to.x If TRUE integrates only to observed x, else integrates to W.
#' @param all If TRUE then returns esw for every observation, else returns only that for
#'        first obs if there are no covariates; always returns esw for every observation
#'        if there are covariates.
#' @param W limit of perp. dist integration. If NULL, uses survey.pars$W.
#' 
#' @details
#' Designed to be called by \code{\link{fitted_invp1}}.
#' Identical to fitted_esw but returns list instead of numerical scalar or vector. This is to allow 
#' it to be used in \code{\link{NDest}} for estimating density and abundance.
#
#' Calls \code{\link{hmltm.esw}} to calclate effective stript width (esw) for 1 observer for fitted 
#' object \code{hmmlt}.
#' 
#' @export
fitted_esw1 <- function(hmmlt,obs=1:nrow(hmmlt$xy),nx=100,to.x=FALSE,all=FALSE,W=NULL){
  if(!is.null(obs)) {
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
  
  # unique sighting identifier
  sID <- list(stratum=hmmlt$xy$stratum,transect=hmmlt$xy$transect,object=hmmlt$xy$object) 
  
  esw <- hmltm.esw1(pars,hfun,models,cov,survey.pars,hmm.pars,ID=sID,nx,to.x,type="response",W=W)
  return(esw)
}


#' @title Calculates esw from model.
#'
#' @description
#' Calculates effective strip width from fitted model.
#' 
#' @param pars parameters (e.g. \code{$fit$par} output from \code{\link{fit.hmltm}})
#' @param hfun detection hazard type (character) 
#' @param models list with elements \code{$y} and \code{$x} speficying covariate model (as for 
#' \code{\link{est.hmltm}}).
#' @param cov data frame with covariates for detections.
#' @param survey.pars survey pars list (as for \code{\link{est.hmltm}}).
#' @param hmm.pars hmm pars list (as for \code{\link{est.hmltm}})
#' @param ID list with elements \code{$stratum}, \code{$transect}, \code{$object} - to uniquely 
#' identify detections.
#' @param nx number of x values for Simpson's rule integration
#' @param type "response" (default) or "link". If "link", interprets pars as being on link scale, 
#' else natural scale
#' @param to.x If TRUE integrates only to observed x, else integrates to W
#' @param W limit of perpendicular dist integration. If NULL, uses \code{survey.pars$W}
#' 
#' @details
#' Designed to be called by \code{\link{fitted_esw1}}.
#' Identical to hmltm.esw but returns list instead of numerical scalar or vector, and allows limit 
#' of integration (W) to be specified explicitly, which overrides limit \code{survey.pars$W}.
#' This is to allow it to be used in \code{\link{NDest}} for estimating density and abundance.
#
#' Returns effective stript half-width (esw) for 1 observer for fitted object \code{hmmlt}, 
#' integrating using Simpson's rule.
#' 
#' @export
hmltm.esw1 <- function(pars,hfun,models,cov,survey.pars,hmm.pars,ID,nx=100,type="response",to.x=FALSE,
                       W=NULL){
  nmax <- nrow(cov)
  if(is.nullmodel(models)) 
    smax <- length(ID$object) 
  else 
    smax <- nmax # number detections
  
  if(to.x) 
    maxx=cov$x
  else {
    maxx <- rep(survey.pars$W,nmax)
    if(!is.null(W)) 
      maxx <- rep(W,nmax)
  }
  
  # vectors for unique stratum, transect, sighting numbers and esws
  ustrat <- utrans <- uobject <- esw <- rep(0,smax) 
  
  for(i in 1:nmax) {
    if(maxx[i]>0){
      xs <- seq(0,maxx[i],length=nx) # set of poiints on which to evaluate p(see|x)
      px <- hmltm.px(xs,pars,hfun,models,cov[i,],survey.pars,hmm.pars,type)
      ustrat[i] <- ID$stratum[i]
      utrans[i] <- ID$transect[i]
      uobject[i] <- ID$object[i] # record sighting ID
      esw[i] <- sintegral(px,xs)
    }
  }
  if(smax>nmax){ # repeat single esw for all smax detections:
    for(i in 2:smax) {
      ustrat[i] <- ID$stratum[i]
      utrans[i] <- ID$transect[i]
      uobject[i] <- ID$object[i] # record sighting ID
      esw[i] <- esw[1]      
    }
  }
  
  return(list(stratum=ustrat,transect=utrans,object=uobject,esw=esw))
}



#' @title Data truncation.
#'
#' @description
#' Left- and/or right-truncates any variable in data frame dat, inserting effort-only record if 
#' truncation removes all detections on a transect.
#' Subtracts left-trunction point off all variable values - after all truncation.

#' 
#' @param dat distance data frame. Must have columns for stratum stratum.area, transect, transect.length, 
#' <any-observer-level_variable>, and if \code{twosit}==TRUE then object as well. These can 
#' have any names, but the names must be specified, via argument \code{colnames}.
#' @param minval left-truncation point.
#' @param maxval right-truncation point.
#' @param twosit If TRUE, assumes this is an mrds-type dataset with two lines per detection.
#' @param colnames name of columns containing stratum, stratum area, transect, transect length, 
#' <any-observer-level_variable>, and if \code{twosit}==TRUE then object as well, IN THIS ORDER, 
#' in a character vector. The default value is 
#' colnames=c("stratum","area","transect","L","x","obs").
#' 
#' @export
truncdat <- function(dat,minval=0,maxval=NULL,twosit=FALSE,colnames=c("stratum","area","transect","L","x","obs")){
  tdat <- dat
  keepcols <- rep(NA,4)
  
  for(i in 1:4) {
    keepcols[i] <- which(names(dat)==colnames[i])
    if(is.null(keepcols[i])) 
      stop("No column in dat called ",colnames[i])
  }
  
  NAs <- dat[1,,drop=FALSE]
  NAs[,-keepcols] <- NA
  xcol <- which(names(dat)==colnames[5])
  
  if(is.null(xcol)) 
    stop("No column in dat called ",colnames[5])
  
  if(is.null(maxval)) 
    maxval <- max(na.omit(dat[,xcol]))
  
  if(twosit){
    if(length(colnames)<6) 
      stop("With double-observer data, need 6th colname for observer; only 5 colnames:",colnames)
    
    obscol <- which(names(dat)==colnames[6])
    
    if(is.null(obscol)) 
      stop("No column in dat called ",colnames[6])
    
    out1 <- which(dat[,obscol]==1 & (dat[,xcol]<minval | dat[,xcol]>maxval))
    out2 <- which(dat[,obscol]==2 & (dat[,xcol]<minval | dat[,xcol]>maxval))
    
    if(length(out1) != length(out2)) 
      stop("Different number of obs1 and obs2 detections for left-truncation")
    if(unique(out2-out1) != 1) 
      stop("Looks like non-consecutive obs1 and obs2 detections for left-truncation")
    
    nout <- length(out1)
    svalues <- tvalues <- rep(NA,nout)
    
    if(is.factor(dat[,keepcols[1]])) {
      svalues <- rep(levels(dat[,keepcols[1]])[1],nout)
      levels(svalues) <- levels(dat[,keepcols[1]])
      NAs[keepcols[1]] <- svalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
      levels(NAs) <- levels(svalues)
    }
    if(is.factor(dat[,keepcols[3]])) {
      tvalues <- rep(levels(dat[,keepcols[3]])[3],nout)
      levels(tvalues) <- levels(dat[,keepcols[3]])
      NAs[keepcols[3]] <- tvalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
      levels(NAs) <- levels(tvalues)
    }
    
    outeff <- data.frame(stratum=svalues,area=rep(NA,nout),transect=tvalues,L=rep(NA,nout))
    
    if(is.factor(outeff$stratum)) 
      levels(outeff$stratum) <- levels(svalues)
    if(is.factor(outeff$transect)) 
      levels(outeff$stratum) <- levels(tvalues)
    
    for(i in 1:nout) 
      outeff[i,] <- dat[out1[i],keepcols] # store effort info for truncated sightings
    
    tdat <- dat[-c(out1,out2),] # remove all truncated sightings
    
    for(i in 1:nout) {
      got <- which(tdat[,keepcols[1]]==outeff$stratum[i] & tdat[,keepcols[2]]==outeff$area[i] & 
                     tdat[,keepcols[3]]==outeff$transect[i] & tdat[,keepcols[4]]==outeff$L[i])
      
      if(length(got)==0) { # transect no longer in truncated dataset
        tdat <- rbind(tdat,NAs) # add row of NAs
        tdat[dim(tdat)[1],keepcols] <- outeff[i,] # add in missing transect info
      }
    }    
  } else {
    out <- which(dat[,xcol]<minval | dat[,xcol]>maxval)
    nout <- length(out)
    svalues <- tvalues <- rep(NA,nout)
    
    if(is.factor(dat[,keepcols[1]])) {
      svalues <- rep(levels(dat[,keepcols[1]])[1],nout)
      levels(svalues) <- levels(dat[,keepcols[1]])
      NAs[keepcols[1]] <- svalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
      levels(NAs) <- levels(svalues)
    }
    
    if(is.factor(dat[,keepcols[3]])) {
      tvalues <- rep(levels(dat[,keepcols[3]])[3],nout)
      levels(tvalues) <- levels(dat[,keepcols[3]])
      NAs[keepcols[3]] <- tvalues[1] # put a factor in the NA row (arbitrary level, as it will be overwritten)
      levels(NAs) <- levels(tvalues)
    }
    outeff <- data.frame(stratum=svalues,area=rep(NA,nout),transect=tvalues,L=rep(NA,nout))
    
    if(is.factor(outeff$stratum)) 
      levels(outeff$stratum) <- levels(svalues)
    
    if(is.factor(outeff$transect)) 
      levels(outeff$stratum) <- levels(tvalues)
    
    for(i in 1:nout) 
      outeff[i,] <- dat[out[i],keepcols] # store effort info for truncated sightings
    
    tdat <- dat[-out,] # remove all truncated sightings
    
    for(i in 1:nout) {
      got <- which(tdat[,keepcols[1]]==outeff$stratum[i] & tdat[,keepcols[2]]==outeff$area[i] & 
                     tdat[,keepcols[3]]==outeff$transect[i] & tdat[,keepcols[4]]==outeff$L[i])
      
      if(length(got)==0) { # transect no longer in truncated dataset
        tdat <- rbind(tdat,NAs) # add row of NAs
        tdat[dim(tdat)[1],keepcols] <- outeff[i,] # add in missing transect info
      }
    }
  }
  
  tdat <- tdat[order(tdat[keepcols[1]],tdat[keepcols[3]]),]
  tdat[!is.na(tdat[,xcol]),xcol] <- tdat[!is.na(tdat[,xcol]),xcol]-minval # shift all left
  
  return(tdat)
}


#' @title Reduces MRDS data frame to CDS data frame.
#'
#' @param dat distance data frame in mrds form. Must have columns 'object' (unique detection 
#' identifier), 'seen' (binary indicating detecteb by observer or not) and 'y' (forward) 
#' detection distance.
#' @param prefer which of the two observers' data to prefer when forward distances are 
#'   missing/equal must be 1 or 2.
#'   
#' @description
#' Reduces mark-recapture Distance sampling (MRDS) data frame dat, with two lines per detection, to a 
#' conventional distance sampling (CDS) data frame with a single line per detection. In the case of 
#' duplicates ("recaptures"), takes the information from the detection made farthest ahead (i.e. that
#' with larger \code{y}) With duplicates that have the same \code{y} or neither of which have a 
#' \code{y}, it chooses according to the parameter \code{prefer}. If only one of the duplicates 
#' has a \code{y}, it chooses that one.
#' 
#' @export
make.onesit <- function(dat,prefer=1) {
  if(prefer!=1 & prefer!=2) 
    stop("Argument 'prefer' must be 1 or 2.")
  
  n <- dim(dat)[1]
  out <- rep(FALSE,n)
  i <- 1
  
  while(i<=n){
    if(!is.na(dat$object[i])){
      if(dat$seen[i]==0) 
        out[i] <- TRUE
      
      if(dat$seen[i+1]==0) 
        out[i+1] <- TRUE
      
      if(dat$seen[i]==1 & dat$seen[i+1]==1) {
        if(is.na(dat$y[i]) & is.na(dat$y[i+1])) # no y's; remove not preferred detection
          out[i+(3-prefer)-1] <- TRUE
        else if(is.na(dat$y[i]) & !is.na(dat$y[i+1])) # keep only non-NA y
          out[i] <- TRUE
        else if(!is.na(dat$y[i]) & is.na(dat$y[i+1])) # keep only non-NA y
          out[i+1] <- TRUE
        else if(dat$y[i]==dat$y[i+1])# remove not preferred detection
          out[i+(3-prefer)-1] <- TRUE
        else if(dat$y[i]<dat$y[i+1]) # remove later (closer) detection
          out[i] <- TRUE
        else
          out[i+1] <- TRUE
      }
      
      i <- i+2
    } else 
      i <- i+1
  }
  
  dat <- dat[!out,] # remove worse observer's lines
  return(dat)
}



#' @title Estimates density and abundance.
#'
#' @description
#' Horvitz-Thompson like estimation of density and abundance of groups and of individuals, as well as
#' of group size (estimated as the ratio of individual density and group density estimates). Produces 
#' estimates by stratum and over all strata.
#' 
#' @param dat MRDS data frame. Must have cols "stratum","area","transect","L","size","object","x","y" 
#' (and possibly others).
#' @param hmltm.fit output from \code{\link{fit.hmltm}}.
#' @param W perpendicular truncation distance for estimation.
#' 
#' @export
NDest <- function(dat,hmltm.fit,W){
  maxx <- max(na.omit(dat$x))
  
  if(maxx>W) {
    cat("Maximum perp. dist=",maxx,"is greater than W=",W,"\n")
    stop("You need a bigger W.")
  }
  
  # Add 1/p column
  dat$invp <- rep(NA,dim(dat)[1])
  invp <- fitted_invp1(hmltm.fit,W=W)
  
  for(i in 1:length(invp$object)) {
    row <- which(dat$stratum==invp$stratum[i] & dat$transect==invp$transect[i] & dat$object==invp$object[i])
    
    if(length(row)>1) {
      cat("Target stratum:",invp$stratum[i],"\n")
      cat("Target transect:",invp$transect[i],"\n")
      cat("Target sighting:",invp$object[i],"\n")
      cat("Found >1: at rows",row,"\n")
      stop("")
    }
    
    dat$invp[row] <- invp$invp[i]
  }
  
  # Calculate density and abundance by stratum
  m2km <- 1/1000
  strat <- unique(dat$stratum)
  nstrat <- length(strat)
  n <- L <- a <- A <- Dg <- D <- Ng <- N <- sbar <- rep(0,nstrat+1)
  stratname <- rep("",nstrat+1)
  
  for(i in 1:nstrat){
    stratname[i] <- as.character(strat[i])
    vdat <- dat[dat$stratum==strat[i],]
    
    trans <- unique(vdat$transect)
    L.tr <- 0
    
    for(tr in 1:length(trans))
      L.tr <- L.tr+vdat$L[min(which(vdat$transect==trans[tr]))]
    
    L[i] <- L.tr
    a[i] <- L[i]*2*W*m2km
    A[i] <- vdat$area[1]
    svdat <- vdat[!is.na(vdat$object),]
    n[i] <- length(svdat$invp)
    Dg[i] <- sum(svdat$invp)/a[i]
    D[i] <- sum(svdat$size*svdat$invp)/a[i]
    sbar[i] <- D[i]/Dg[i]
    Ng[i] <- Dg[i]*A[i]
    N[i] <- D[i]*A[i]
  }
  
  stratname[nstrat+1] <- "Total"
  Ng[nstrat+1] <- sum(Ng[1:nstrat])
  N[nstrat+1] <- sum(N[1:nstrat])
  A[nstrat+1] <- sum(A[1:nstrat])
  Dg[nstrat+1] <- Ng[nstrat+1]/sum(A[1:nstrat])
  D[nstrat+1] <- N[nstrat+1]/sum(A[1:nstrat])
  n[nstrat+1] <- sum(n[1:nstrat])
  L[nstrat+1] <- sum(L[1:nstrat])
  a[nstrat+1] <- sum(a[1:nstrat])
  sbar[nstrat+1] <- D[nstrat+1]/Dg[nstrat+1]
  
  # add transect frequency:
  tfreq <- apply(table(dat$stratum,dat$transect)>0,1,sum)
  k <- c(tfreq,sum(tfreq))
  
  return(list(invp=invp,
              ests=data.frame(stratum=stratname,n=n,k=k,L=L,covered.area=a,stratum.Area=A,
                              Dgroups=signif(Dg,3),Ngroups=signif(Ng,3),mean.size=round(sbar,1),
                              D=signif(D,5),N=round(N,1))
  )
  )
}

#' @title Line transect estimation with a hidden Markov availability model.
#'
#' @description
#' \code{est.hmltm} estimates group and individual density and abundance, together with mean 
#' group size, by stratum, from (1) line transect data that includes forward detection distances and 
#' (2) estimated Markov model or hidden Markov model availability prameters. 
#'
#' @param dat data frame in distance-like format, but including forward distances of detections. The 
#' following are compulsory elements (field name in quotes, contents in brackets): 
#' "stratum" (survey stratum: must be numeric), "area" (stratum area), "transect" (transect number: must
#' be numeric), "L" (transect length), "size" (group size), "object" (unique detection identifier: must
#' be numeric), "x" (perpendicular distance), "y" (forward distance).
#' @param pars starting parameter values.
#' @param FUN detection hazard functional form name (character). Currently implemented forms are 
#' "h.IP.0", "h.EP1.0", "h.EP2.0", "h.EP1x.0", "h.EP2x.0". (See Vignette "Specifying models and parameter 
#' starting values" for details.)
#' @param models list of characters with elements \code{$y} and \code{$x} specifying models for the y- 
#' and x-dimension detection hazard scale parameters. Must be either \code{NULL} or regression model 
#' format (without response on left, e.g. "~size").
#' @param survey.pars a list containing the following elements (in any order):
#'  \itemize{
#'  \item {$spd} {speed of observer,}
#'  \item {$W} {perpendicular distance right-truncation point,}
#'  \item {$Wl} {perpendicular distance left-truncation point,}
#'  \item {$ymax} {forward distance by which detection probability is effectively zero,}
#'  \item {$dT} {availability process (Markov chain) time step size.}
#' }
#' @param hmm.pars a list containing the parameters of animals' availability processs hidden Markov 
#' model (HMM), as follows (in any order):
#' \itemize{
#'  \item {$Pi} {a 2x2xm HMM transition probability matrix, where m is the number of availability HMMs 
#'  being used to model animal availability. If m>1, each set of HMM parameters is treated as a random
#'  sample from the set of HMM parameters in the population.}
#'  \item {$pm} {a2xm matrix of HMM state-dependent Bernoulli distribution parameters (the probabilities
#'  of being available, given the animal's "behavioural" state - i.e. the state of the hidden Markov 
#'  chain)}
#'  \item {$delta} {a 2xm matrix of stationary distribution of a Markov chain, the ith of which has
#'  transition probability matrix Pi[,,i].}
#' }
#' And if the HMM was constructed from mean times animals are available and unavailable (by means of
#' function \code{\link{make.hmm.pars.from.Et}} for example), then also
#' \itemize{
#'  \item {$Et} {a 2xm matrix in which the first element is the mean time animals are UNavailable
#'  in a single available-unavailable cycle, and the second element is the corresponding mean time that
#'  they are available,}
#'  \item {Sigma.Et} {a 2x2xm matrix, in which Sigma.Et[,,i] is the variance-covariance matrix of 
#'  Et[,i.] (i.e. the variance-covariance matrix of Et for the ith availability model).}
#' }
#' @param control.fit list with elements
#' \itemize{
#'  \item{$hessian} {logical) - if TRUE Hessian is estimated and returned, else not,}
#'  \item{$nx} {(scalar) - the number of intervals to use with Simpson's rule integration over y. 
#'    \code{nx=64} seems safe; smaller number makes computing faster.}
#' }
#' @param control.opt as required by \code{\link{optim}} (and hence by \code{\link{fit.hmltm}}).
#' @param twosit TRUE if \code{dat} is in mrds format (with two lines per detection), else assumes
#'              that \code{dat} is in cds format (with one line per detection).
#' @param notrunc if TRUE, does not do any perp dist truncation, else uses \code{survey.pars$W}
#' and \code{$Wl} to do perp dist truncation.
#' @param W.est right truncation perpendicular distance for estimation. Can't be less than maximum 
#' perpendicular distance (x) in the line transect data frame \code{dat}, but can be less than 
#' the max perpendicular distance used for fitting (\code{survey.pars$W}).
#' @param groupfromy a forward distance (y) below which all y's are grouped into a single
#' interval in the likelihood function (i.e. exact y,s < groupfromy are combined into
#' an interval rather than passed as exact distances).
#'
#' @return A list with four elements: \code{hmltm.fit}, \code{point}, \code{dat}, \code{W.est}. 
#' Their contents are as follows:
#' 
#' \code{hmltm.fit} is the output from \code{fit.hmltm}, i.e. a list containing the following elements:
#' \itemize{
#'  \item{xy} {dat used in fitting (input reflection).}
#'  \item{phats} {estimated detection probabilities of all detections.}
#'  \item{phat} {1/mean(1/phat).}
#'  \item{pzero} {estimated detection probabilities at perpendicular distance.}
#'  \item{h.fun} {=FUN (input reflection).}
#'  \item{models} {=models (input reflection).}
#'  \item{fit} {output from \code{\link{fit.xy}}.}
#'  \item{Loglik} {log-likelihood function at MLE.}
#'  \item{AIC} {AIC.}
#'  \item{x} {vector of x-values for plotting perpendicular distance fit.}
#'  \item{p} {vector of detection function values for plotting perpendicular distance fit.}
#'  \item{fitpars} {a list containing all the given parameters controlling the fit (survey.pars,hmm.pars,
#'  control.fit,control.optim).}
#' }
#' \code{point} is a list containing two elements:
#' \itemize{
#'  \item{invp} {is a data frame containing one row for every observation, with the first three columns
#'   giving the stratum, transect and object identifier for the observation, and the final column 
#'   (invp) giving the estimate of the inverse of the probability of detection for the observation.}
#'  \item{ests} {is a data frame with one row per stratum and a final row for all strata combined,
#'  and columns giving the number of detections in the stratum (n), the line lingth in the stratum 
#'  (L), the covered area in the stratum (covered.area=2WL), the stratum area (stratum.Area), the 
#'  estimated group density in the stratum (Dgroups), the estimated group abunance in the stratum 
#'  (Ngroups), the estimated mean group size in the stratum (mean.size), the individual denstiy
#'  in the stratum (D), and the abundance in the stratum (N).}
#' }
#' \code{dat} is the data frame passed to \code{est.hmltm}.
#' 
#' \code{W.est} is the right perpendicular distance used for estimation (and passed to 
#' \code{est.hmltm}.)
#' 
#' @references Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 
#' 2013. Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' @export
#' @useDynLib hsltm
est.hmltm <- function(dat,pars,FUN,models=list(y=~NULL,x=~NULL),survey.pars,hmm.pars,
                      control.fit,control.opt,twosit=FALSE,notrunc=FALSE,W.est=NULL,
                      groupfromy=NULL)
{
  # If no covariate in either dimension, set to NULL formula
  if(is.null(models$y))
    models$y <- ~NULL
  if(is.null(models$x))
    models$x <- ~NULL
  
  # data truncation:
  if(!notrunc) {
    if(min(dat$x,na.rm=TRUE) < survey.pars$Wl | survey.pars$W < max(dat$x,na.rm=TRUE))
      dat <- truncdat(dat,minval=survey.pars$Wl,maxval=survey.pars$W,twosit=FALSE)
  }
  
  srows <- !is.na(dat$object)
  sdat <- dat[srows,]
  
  # fit the model:
  hmltm.fit <- fit.hmltm(sdat,pars,FUN,models,survey.pars,hmm.pars,control.fit,control.opt,
                         groupfromy=groupfromy)
  
  # estimate density, etc:
  if(is.null(W.est)) 
    W.est <- survey.pars$W
  
  if(!is.null(survey.pars$Wl)) 
    W.est <- survey.pars$W-survey.pars$Wl # since in this case data all shifted left by $Wl
  
  point <- NDest(dat,hmltm.fit,W.est)
  hmltm.obj <- list(hmltm.fit=hmltm.fit,point=point,dat=dat,W.est=W.est)
  class(hmltm.obj) <- c(class(hmltm.obj),"hmltm")
  
  return(hmltm.obj)
}
