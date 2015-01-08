# AIC and related summary functions
#==================================

summary.hmltm=function(object,...,sort=TRUE,k=2,dmax=10,criterion=c("AICc","AIC")) {
    allargs <- list(...)
    modelnames <- (c(as.character(match.call(expand.dots = FALSE)$object), 
                     as.character(match.call(expand.dots = FALSE)$...)))
    allargs <- hmltmlist(object, allargs)
    names(allargs) <- modelnames
    p0=sapply(allargs,getp0)
    Nhat=sapply(allargs,getNhat)
    Ngrp=sapply(allargs,getNgroups)
    grpsize=sapply(allargs,getsize)
    outp=AIC(allargs, sort = sort, k = k, dmax = dmax, criterion = criterion)
    outp$Nhat.grp=Ngrp
    outp$grp.size=grpsize
    outp$Nhat=Nhat
    outp$p0=p0
    return(outp)
}

getp0=function(object) {
  if(!inherits(belIP.wait,"hmltm")) stop("object must be of class `hmltm'")
  return(object$hmltm.fit$pzero)
}
getNhat=function(object) {
  if(!inherits(belIP.wait,"hmltm")) stop("object must be of class `hmltm'")
  return(object$point$ests[dim(object$point$ests)[1],dim(object$point$ests)[2]])
}
getNgroups=function(object) {
  if(!inherits(belIP.wait,"hmltm")) stop("object must be of class `hmltm'")
  return(object$point$ests[dim(object$point$ests)[1],7])
}
getsize=function(object) {
  if(!inherits(belIP.wait,"hmltm")) stop("object must be of class `hmltm'")
  return(object$point$ests[dim(object$point$ests)[1],8])
}

AIC.hmltm=function (object,...,sort=TRUE,k=2,dmax=10,criterion=c("AICc","AIC")) 
{
  allargs <- list(...)
  modelnames <- (c(as.character(match.call(expand.dots = FALSE)$object), 
                   as.character(match.call(expand.dots = FALSE)$...)))
  allargs <- hmltmlist(object, allargs)
  names(allargs) <- modelnames
  AIC(allargs, sort = sort, k = k, dmax = dmax, criterion = criterion)
}

AIC.hmltmlist=function (object,...,sort=TRUE,k=2,dmax=10,criterion= c("AICc","AIC")) 
{
  if (k != 2) 
    stop("AIC.hmltm defined only for k = 2")
  if (length(list(...)) > 0) 
    warning("... argument ignored in 'AIC.hmltmlist'")
  criterion <- match.arg(criterion)
  modelnames <- names(object)
  allargs <- object
  if(any(!sapply(allargs,inherits,"hmltm")))
    stop("components of 'object' must be 'hmltm' objects")
  output <- data.frame(t(sapply(allargs, oneline.hmltm)), stringsAsFactors = F)
  output$delta <- output[, criterion] - min(output[, criterion])
  OK <- abs(output$delta) < abs(dmax)
  sumdelta <- sum(exp(-output$delta[OK]/2))
  output$wt <- ifelse(OK, round(exp(-output$delta/2)/sumdelta,4), 0)
  row.names(output) <- modelnames
  if (sort) output <- output[order(output[, criterion]), ]
  names(output)[5] <- paste("d", criterion, sep = "")
  names(output)[6] <- paste(criterion, "wt", sep = "")
  if (nrow(output) == 1) {
    output[, 6] <- NULL
    output[, 5] <- NULL
  }
  output
}

oneline.hmltm=function(object) {
  n=length(obj$hmltm.fit$phats) # number detections
  loglik=logLik(object) # logLik structure
  nll=-as.numeric(loglik)
  p=as.integer(attributes(loglik)$df) # number of parameters
  AICval=2*(nll+p) # AIC
  AICcval=ifelse((n-p-1)> 0, 2*(nll+p)+2*p*(p+1)/(n-p-1), NA) # AICc
  c(npar=p,logLik= -nll, AIC=round(AICval, 3), AICc=round(AICcval,3))
}


hmltmlist=function (...) 
{
  dots <- match.call(expand.dots = FALSE)$...
  allargs <- list(...)
  allargs <- lapply(allargs, function(x) if (inherits(x, "hmltm")) 
    list(x)
    else x)
  temp <- do.call(c, allargs)
  if (is.null(names(temp))) 
    names(temp) <- paste("hmltm", 1:length(temp), sep = "")
  if (!all(sapply(temp, function(x) inherits(x, "hmltm")))) 
    stop("objects must be of class 'hmltm' or 'hmltmlist'")
  class(temp) <- "hmltmlist"
  temp
}

logLik.hmltm=function(object, ...) 
{
  npar <- length(object$hmltm.fit$fit$par)
  structure(object$hmltm.fit$Loglik, df = npar, class = "logLik")
}


