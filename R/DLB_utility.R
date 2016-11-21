#                      ------------------------
#--------------------- DLB's utility functions --------------------

#' @title Clone of \code{sample} without weird behaviour for length(x)==1.
#'
#' @usage bsample(x, size, replace = FALSE, prob = NULL)
#' 
#' @description
#'  Identical to \code{sample} but when x is an integer it just returns x rather than an
#'  integer in 1:x (which is what \code{sample} does).
#'  
#' @param x Either a vector of one or more elements from which to choose, or a positive integer. 
#' See \link{sample} for details.
#' @param size a positive number, the number of items to choose from.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the elements of the vector being sampled.
bsample <- function(x,size,replace=FALSE,prob=NULL) {
  if(length(x)==1) 
    return(x)
  else 
    return(sample(x,size,replace,prob))
}


#' @title Caclulates coefficient of variation.
#'
#' @description
#'  Utility function
#'  
#' @param x Random variable.
cv <- function(x) sd(x)/mean(x) # calculates coefficient of variation

#' @title Converts nm to metres.
#'
#' @description
#'  Utility function
#'  
#' @param x Distance in nautical miles.
nm2m <- function(x) return(x*1852) # converts nautical miles to metres


#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
#' 
#' @export
histline <- function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                     xlab="x",ylab="y",...)
{
  n <- length(height)
  if(length(breaks)!=(n+1)) 
    stop("breaks must be 1 longer than height")
  
  if(outline) {
    y <- c(0,rep(height,times=rep(2,n)),0)
    x <- rep(breaks,times=rep(2,(n+1)))
  } else {
    y <- rep(0,4*n)
    x <- rep(0,4*n+2)
    
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)] <- c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)] <- c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x <- x[1:(4*n)]
  }
  
  if(lineonly) {
    if(!fill) 
      lines(x,y,...)
    else 
      polygon(x,y,...)
  } else {
    if(!fill) 
      plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}


#' @title Integration using Simpson's method.
#'
#' @description
#'  Modified version of function sintegral from library Bolstad.
#'  
#' @param fx Values at \code{x} (see below) of the function to be integrated.
#' @param x Values at which function fx has been evaluated.
#' @param n.pts The number of points to be used in integration. If \code{x} contains more than
#'\code{n.pts} then \code{n.pts} will be set to \code{length(x)}.
#' @param type if equal to `\code{cdf}' a list comprising x-values and cdf values at each x-value 
#' is returned, else the integral over the range of \code{x} is returned.
sintegral <- function (fx, x, n.pts=16, type="int") 
{
  #  if (class(fx) == "function") 
  #    fx = fx(x)
  #  n.x = length(x)
  #  if (n.x != length(fx)) 
  #    stop("Unequal input vector lengths")
  #  if (n.pts < 64) 
  #    n.pts = 64
  ap <- approx(x, fx, n = 2 * n.pts + 1)
  h <- diff(ap$x)[1]
  integral <- h*(ap$y[2*(1:n.pts)-1]+4*ap$y[2*(1:n.pts)]+ap$y[2*(1:n.pts)+1])/3
  
  if(type!="cdf") 
    return(sum(integral)) 
  else 
    return(list(x=ap$x[2*(1:n.pts)],y=cumsum(integral)))
}


# ============== Other availability correction methods and realted ==========================

#' @title Laake's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the probability that an animal is available at least once while within detectable
#' region, using the method of Laake et al. (1997), Eqn (4). Does this for each of a set of m 
#' 2-state Markov model availability parameters passed to it and returns probabilities and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param ymax is maximum forward distance to consider (`w' in Laake et al. (1997), Eqn (4)).
#' @param spd speed of observer. Needed to convert time (units of \code{hmm.pars}) to distance (units
#'  of \code{ymax}).
#'
#' @details See Laake et al. (1997), Eqn (4), or Borchers et al. (2013) for details of Laake's method.
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' Laake, J., Calambokidis, J., Osmek, S., and Rugh, D. 1997. Probability of detecting harbor 
#' porpoise from aerial surveys: estimating g(0). Journal of Wildlife Management 61, 63-75.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{mclaren.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' 
#' # because hmm.pars is in units of time, not distance, you need to specify spd.
#' laake.a(hmm.pars,ymax=200,spd=4) 
#' 
#' @export

laake.a <- function(hmm.pars,ymax,spd=NULL){
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a <- rep(NA,nav)
  
  for(i in 1:nav)
    a[i] <- jeffa(hmm.pars$Pi[,,i],ymax,spd)
  
  return(list(mean=mean(a),a=a))
}


#' @title Laake's availability correction factor calculation for a single 2-state Markov model.
#'
#' @description
#' Calculates the probability that an animal is availabel at least once whil within detectable
#' region, using the method of Laake et al. (1997), Eqn (4).
#'
#' @param Pi is Markov model transition matrix (state 1 is UNavailable).
#' @param w is maximum forward distance things can be detected.
#' @param spd is observer speed. If NULL w is assumed to be in the units of the Markov chain.
#' @param E is mean times unavailable (E[1]) and available (E[2]) (in unist of Markov chain). 
#' If E is not NULL, it is used in preference to Pi for calculations, else it is ignored.
#'
#' @details See Laake et al. (1997), Eqn (4), or Borchers et al. (2013) for details of the method.
#' 
#' @references 
#' Borchers, D.L., Zucchini, W., Heide-Jorgenssen, M.P., Canadas, A. and Langrock, R. 2013. 
#' Using hidden Markov models to deal with availability bias on line transect surveys. Biometrics.
#' 
#' Laake, J., Calambokidis, J., Osmek, S., and Rugh, D. 1997. Probability of detecting harbor 
#' porpoise from aerial surveys: estimating g(0). Journal of Wildlife Management 61, 63-75.
#' 
#' @export
jeffa <- function(Pi,w,spd=NULL,E=NULL)
{
  if(!is.null(spd)) 
    w <- w/spd # convert distance to time if spd given
  
  if(!is.null(E))
    warning("Used E, not Pi for calculations")
  else
    E <- makeE(Pi) # expected time up, expected time down
  
  a <- (E[2] + E[1]*(1-exp(-w/E[1])))/sum(E)    
  return(a)
}


#' @title Instantaneous availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the probability that an animal is available at an instant, for each of a set of m 
#' 2-state Markov model availability parameters passed to it and returns probabilities and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param Et is a 2xm matrix with expect times Unavailable (row 1) and Available (row 2).
#'
#' @details If Et is given, it is used and hmm.pars is ignored.
#' 
#' @seealso \code{\link{mclaren.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3);pm=NULL
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' instant.a(hmm.pars)
#' instant.a(NULL,Et=matrix(c(Eu,Ea),ncol=2,byrow=TRUE))
#' instant.a(NULL,c(Eu[1],Ea[1]))
#' 
#' @export
instant.a <- function(hmm.pars,Et=NULL){
  if(!is.null(Et)) {
    if(!is.null(hmm.pars)) 
      warning("Used Et, not Pi for calculations")
    
    if(is.vector(Et)) 
      Et <- matrix(Et,ncol=1)
    
    if(dim(Et)[1]!=2) 
      stop("1st dimension of Et must be 2.")
    
    nav=dim(Et)[2]
  } else {
    if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) 
      stop("1st two dimensions of hmm.pars$Pi must be 2.")
    
    if(length(dim(hmm.pars$Pi))==2) 
      hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
    
    nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  }
  
  a <- rep(NA,nav)
  
  for(i in 1:nav) {
    if(is.null(Et)) 
      a[i] <- simplea(hmm.pars$Pi[,,i],Et)
    else 
      a[i] <- simplea(NULL,Et[,i])
  }
  
  return(list(mean=mean(a),a=a))
}


#' @title Simple availability correction factor calculation.
#'
#' @description
#' Calculates proportion of time an animal with a Markov availability process is available.
#'
#' @param Pi is Markov model transition matrix (state 1 is UNavailable).
#' @param E is expected time in each state (state 1 in UNavailable).
#'
#' @details If \code{E} is NULL, uses Pi to calculate proportion of time available, else uses \code{E}.
#' 
#' @export
simplea <- function(Pi,E=NULL)
{
  if(is.null(E)) 
    E <- makeE(Pi) # expected time up, expected time down
  a <- E[2]/sum(E)    
  
  return(a)
}

#' @title McLaren's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates McLaren's availability correction factor, for each of a set of m 2-state Markov model 
#' availability parameters passed to it and returns these and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element \code{\$Pi} (where m is number of availability parameter sets).
#' @param w is max forward distance things can be seen at (or max forward time). Must be scalar.
#' @param spd is observer speed; omit if w is max forward TIME.
#'
#' @references
#' McLaren, I.A. 1961. Methods of determining the numbers and availability of ringed seals in the 
#' eastern Canadian Arctic. Arctic 14:162--175.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{laake.a}}, \code{\link{richard.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#' mclaren.a(hmm.pars,w=10,spd=4)
#' mclaren.a(hmm.pars,w=100,spd=4) # can be greater than 1 (!)
#' 
#' @export
mclaren.a <- function(hmm.pars,w,spd=1){
  if(dim(hmm.pars$Pi)[1]!=2 | dim(hmm.pars$Pi)[2]!=2) 
    stop("1st two dimensions of hmm.pars$Pi must be 2.")
  
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  a <- rep(NA,nav)
  
  for(i in 1:nav){
    Eti <- makeE(hmm.pars$Pi[,,i])*spd # Pi is TIME; if w is DISTANCE need to multiply time by speed
    a[i] <- (Eti[2]+w)/sum(Eti)
  }
  
  return(list(mean=mean(a),a=a))
}

#' @title Richard's availability correction factor calculation for multiple availability models.
#'
#' @description
#' Calculates the availability correction factor "C_{ca}" of Richard et al. (2010). Does this for each 
#' of a set of m 2-state Markov model availability parameters passed to it and returns probabilities 
#' and their mean.
#'
#' @param hmm.pars is a list with 2x2xm Markov model transition matrices (in which state 1 is UNavailable) 
#' in element $Pi (where m is number of availability parameter sets).
#' @param w vector of forward distances.
#' @param spd observer speed: must be entered if y is not time, since hmm.pars always time.
#'
#' @details See Richard et al. (2010), equation on botto mof page 91 for details of method.
#' 
#' @references 
#' Richard, P.R., Laake, J.L., Hobbs, R.C., Heide-Jorgensen, M.P., Asselin, N.C. and Cleator, H. 2010.
#' Baffin Bay narwhal population distribution and numbers: aerial surveys in the Canadian High Arctic,
#' 2002-04. Arctic 63: 85-99.
#' 
#' @seealso \code{\link{instant.a}}, \code{\link{mclaren.a}}, \code{\link{laake.a}}
#' 
#' @examples
#' Ea=c(10,12);Eu=c(20,22);seEa=c(2,3);seEu=c(4,6);covEt=c(2,3) # mean avail and unavail times (& var)
#' hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt) # make hmm.pars object
#' richard.a(hmm.pars,w=10,spd=4)
#' richard.a(hmm.pars,w=100,spd=4) # can be greater than 1 (!)
#' richard.a(hmm.pars,w=rexp(20,1/100),spd=4)
#' 
#' @export
richard.a <- function(hmm.pars,w,spd=1){
  y <- na.omit(w)
  n <- length(y)
  
  if(length(dim(hmm.pars$Pi))==2) 
    hmm.pars$Pi <- array(hmm.pars$Pi,dim=c(2,2,1)) # need 3D array below
  
  nav <- dim(hmm.pars$Pi)[3] # number of HMM parameter sets
  ina <- matrix(rep(instant.a(hmm.pars)$a,n),nrow=n,byrow=TRUE)
  mca <- matrix(rep(NA,n*nav),nrow=n)
  
  for(i in 1:n) 
    mca[i,] <- mclaren.a(hmm.pars,y[i],spd)$a
  
  sumfb <- apply(ina/mca,2,sum)
  a <- 1/((1/ina[1,])*(sumfb/n))
  
  return(list(mean=mean(a),a=a))
}


#' @title Makes Markov transition matrices from mean times available and unavailable.
#'
#' @description
#' Makes Markov transition matrices from mean times available (Ea) and unavailable (Eu). 
#' If Ea and Eu are vectors of length m (they must be the same length), returns a 2x2xm array 
#' in which element [,,i] is the ith Markov transition matrix; else returns a single 2x2 matrix.
#'
#' @param Eu is the mean time UNavailable in one available-unavailable cycle.
#' @param Ea is the mean time available in one available-unavailable cycle.
#'
#' @examples
#' Ea=c(10,12);Eu=c(20,22)
#' makePi(Eu,Ea)
#' makePi(Eu[1],Ea[1])
#' 
#' @export
makePi <- function(Eu,Ea)
{
  nav <- length(Eu)
  if(length(Ea)!=nav) 
    stop("Lengths of Eu and Ea must be the same")
  #  if(nav==1) {
  #    Pi=matrix(rep(0,4),nrow=2,
  #              dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available")))
  #    Pi[1,2]=1/Eu
  #    Pi[2,1]=1/Ea
  #    Pi[1,1]=1-Pi[1,2]
  #    Pi[2,2]=1-Pi[2,1]
  #  } else {
  Pi <- array(rep(0,2*2*nav),dim=c(2,2,nav),
              dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                            Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Pi[1,2,i] <- 1/Eu[i]
    Pi[2,1,i] <- 1/Ea[i]
    Pi[1,1,i] <- 1-Pi[1,2,i]
    Pi[2,2,i] <- 1-Pi[2,1,i]    
  }
  #  }
  
  return(Pi)
}

#' @title Returns expected time in state for 2-state Markov transition matrices.
#'
#' @description
#' Returns expected time in state for the m 2-state Markov transition matrices Markov transition 
#' matrices contained in the matrix or array Pi. If m>1 returns a matrix with column i being the
#' expected time in state 1 ("unavailable") as its first element and expected time in state 2 
#' ("available") as its second; else returns a vector.
#'
#' @param Pi either a 2x2 Markov transition probability matrix or a 2x2xm array with element [,,i]
#' being a 2x2 Markov transition probability matrix.
#'
#' @examples
#' Ea=c(10,12);Eu=c(20,22)
#' Pi2=makePi(Eu,Ea) # make a set of transition matrices
#' makeE(Pi2) # recover c(Eu,Ea)
#' Pi1=makePi(Eu[1],Ea[1]) # make a single transition matrix
#' makeE(Pi1) # recover c(Eu,Ea)
#' 
#' @export
makeE <- function(Pi){
  #----------------------------------------------------------
  # Returns expected time in states 1 and 2 for the 2x2 
  # probability transition matrix Pi for a 2-state
  # Markov process.
  #----------------------------------------------------------
  if(length(dim(Pi))==2){
    E <- c(1/Pi[1,2],1/Pi[2,1])
    names(E) <- c("Unavailable","Available")
  } else {
    nav <- dim(Pi)[3]
    E <- matrix(rep(NA,2*nav),nrow=2,
                dimnames=list(State=c("Unavailable","Available"),Animal=as.character(1:nav)))
    
    for(i in 1:nav)
      E[,i] <- c(1/Pi[1,2,i],1/Pi[2,1,i])
  }
  
  return(E)
}


#' @title Makes 2-state hidden Markov model list suitable for passing to 
#' \code{\link{est.hmltm}} via argument \code{hmm.pars}.
#'
#' @description
#' \code{make.hmm.pars.from.Et} creates 2-state hidden Markov availability model specification from 
#' mean times unavailable (Eu) and available (Ea).
#'
#' @param Ea vector of length m>=1 specifying mean distance (=mean time * observer speed) animals are in 
#' the more available state in one cycle (e.g. dive cycle: surface-dive).
#' @param Eu vector of length m>=1 specifying mean distance (=mean time * observer speed) animals are in 
#' the more UNavailable state in one cycle.
#' @param seEa standard error of Ea.
#' @param seEu standard error of Eu.
#' @param covEt vector of length m>=1 containing covariance of each pair of Ea and Eu.
#' @param pm is 2xm matrix containing state-dependent Bernoulli distribution parameters for the m
#' pairs of Ea and Eu, with first being probability of being available when in state i (i=1,2), 
#' where i=1 is the more UNavailable state and i=2 is the more available state.
#'
#' @details Calculates 2-state hidden Markov model parameters such that the Markov process is in states 1 
#' (more unavailable) and 2 (more available) for mean times Ea and Eu, with Bernoulli state-dependent response
#' probability pm[1,] and pm[2,], respectively. Also constructs a covariance matrix for Ea and Eu. If 
#' pm[1,i]=0 and pm[2,i]=1 for animal i, the availability process is a Markov process and the states 
#' are actual unavailbile (state 1) and availabile (state 2).
#' 
#'  @examples
#'  # Some arbitrary numbers for illustration:
#'  Ea=c(10,12);Eu=c(100,120);seEa=c(2,3);seEu=c(20,30);covEt=c(10,12)
#'  make.hmm.pars.from.Et(Ea[1],Eu[1],seEa[1],seEu[1],covEt[1]) # single animal
#'  make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt) # two animals
#'  
#'  # Here's how the data porpoise.hmm.pars was created from numbers in Westgate et al. (1995):
#'  ppn=c(40,52,36,34,49,60,33)/100 # proportion of time available
#'  ET=c(76,44,52,64,70,46,103) # mean dive cycle duration
#'  seET=c(48,37,52,65,59,32,67) # SE of mean dive cycle duration
#'  cvET=seET/ET # CV of mean dive cycle duration
#'  Ea=ET*ppn  # mean time available
#'  Eu=ET*(1-ppn)  # mean time UNavailable
#'  # For lack of better info, assume independence of Ea and Eu, and that cv(Ea)=cv(Eu)=cv, 
#'  # which means that
#'  cv=sqrt((cvET*ET)^2/(Ea^2+Eu^2))
#'  # and hence:
#'  seEa=Ea*cv
#'  seEu=Eu*cv
#'  covEt=seET*0 # assume independence
#'  porpoise.hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu,covEt)
#'  
#'  # Here's how the dataset beaked.hmm.pars were created:
#'  Ea=121.824
#'  Eu=1580.256
#'  seEa=9.618659
#'  seEu=134.9212
#'  beaked.hmm.pars=make.hmm.pars.from.Et(Ea,Eu,seEa,seEu)
#'  
#' @references 
#' Westgate, A. J., Read, A. J., Berggren, P., Koopman, H. N., and Gaskin, D. E. 1995. Diving behaviour 
#' of harbour porpoises, Phocoena phocoena. Canadian Journal of Fisheries and Aquatic Sciences 52, 
#' 1064-1073.
#' 
#' @export
make.hmm.pars.from.Et <- function(Ea,Eu,seEa,seEu,covEt=0,pm=NULL) {
  nav <- length(Ea)
  
  if(length(Eu)!=nav |length(seEa)!=nav |length(seEu)!=nav |length(covEt)!=nav) 
    stop("Lengths of Ea, Eu, seEa, seEu, covEt must all be the same.")
  
  if(is.null(pm)) 
    pm <- matrix(c(rep(0,nav),rep(1,nav)),nrow=2,byrow=TRUE)
  
  if(is.vector(pm)) {
    if(length(pm)!=2) 
      stop("pm must either be a vector of length 2 or a matrix of dimension length(Ea)x2.")
    
    pm <- matrix(c(pm[1],pm[2]),ncol=2)
  }
  
  if(dim(pm)[2]!=nav) 
    stop("Inconsistent dimensions of Ea and pm.")
  
  Pi <- Sigma.Et <- array(rep(NA,2*2*nav),dim=c(2,2,nav),
                          dimnames=list(From=c("Unavailable","Available"),To=c("Unavailable","Available"),
                                        Animal=as.character(1:nav)))
  
  Et <- delta <- newpm <- matrix(rep(NA,2*nav),ncol=nav,dimnames=list(State=c("Unavailable","Available"),
                                                                      Animal=as.character(1:nav)))
  for(i in 1:nav) {
    Et[,i] <- c(Eu[i],Ea[i])
    Sigma.Et[,,i] <- diag(c(seEu[i],seEa[i])^2)
    #    cvEt=c(seEu[i]/Et[1,i],seEa[i]/Et[2,i])
    #    Sigma.Et[,,i]=diag((cvEt*Et)^2)
    Sigma.Et[1,2,i] <- Sigma.Et[2,1,i] <- covEt[i]
    Pi[,,i] <- makePi(Et[1,i],Et[2,i])
    delta[,i] <- compdelta(Pi[,,i])
    newpm[,i] <- pm[,i]
  }
  
  hmm.pars <- list(pm=newpm,Pi=Pi,delta=delta,Et=Et,Sigma.Et=Sigma.Et)
  return(hmm.pars)  
}


# ==============- utility functions specific to avail estimation ------------------

#' @title Makes survey parameter list (\code{survey.pars}) suitable for passing to 
#' \code{\link{est.hmltm}}.
#'
#' @description
#' \code{make.survey.pars} just puts a bunch of stuff in a list in a format that \code{\link{est.hmltm}}
#'  expects.
#'
#' @param spd is observer speed.
#' @param W is perpendicular right-truncation distance for estimation.
#' @param ymax is maximum forward distance for estimation - it should be far enough ahead that it is
#' not possible to detect anyting beyond \code{ymax}.
#' @param Wl is perpendicular left-truncation distance for estimation. (After truncating \code{Wl} is 
#' subtracted from all perpendicular distances.)
#' @param dT is the time step on which the availability hidden Markov model operates.
#'
#' @details Packs the above in a list suitable for passing as \code{survey.pars} to 
#' \code{\link{est.hmltm}}.
#' 
#' @export
make.survey.pars <- function(spd,W,ymax,Wl=0,dT=1){
  return(list(spd=spd,W=W,ymax=ymax,Wl=Wl,dT=dT,dy=spd*dT))
}

#' @title Decides if model is a null model.
#'
#' @description
#'  Logical function: true if \code{models} includes no covariates
#'  
#' @param models list of characters with elements \code{$y} and \code{$x} specifying y- and 
#' x-covariate models. Either \code{NULL} or regression model format (without response on left).
#' 
#' @export
is.nullmodel <- function(models){
  null <- TRUE
  
  for(i in 1:length(models)) 
    null <- null & (is.null(models[[i]]) | models[[i]]==~NULL)
  
  return(null)
}


#' @title Constructs linear predictor.
#'
#' @description
#' Returns parameter vector on linear predictor scale for hazard function \code{FUN}, using 
#' \code{model} and data frame \code{dat}. Makes a matrix with each row the parameter values for an 
#' observation, then concatenates these into a single vector (so can pass as vector to C++ code).
#'  
#' @param b parameter vector.
#' @param FUN detection hazard function name (character).
#' @param models list with two components (\code{$y} and \code{$x}) specifying models for detection 
#' hazard function scale parameters in forward (y) and perpendicular (x) dimensions. If detection 
#' hazard function form does not have separate scale parameters in y- and x- dimensions, the model
#' given in \code{$y} is take as the scale parameter model. \code{$y} and \code{$x} must be either 
#' \code{NULL} or "~<regspec>", where <regspec> is a regression model spefication (e.g. 
#' \code{height+weight} or \code{height:weight} or \code{height*weight}, etc.).
#' @param dat data frame, which must have columns corresponding to variable names in \code{$y} and
#' \code{$x}.
#' 
#' @export
makeCovPar <- function(b, FUN, models, dat)
{
  nbObs <- nrow(dat)
  
  cov_x <- model.matrix(as.formula(models$x),data=dat)
  nbCov_x <- ncol(cov_x) - 1
  cov_y <- model.matrix(as.formula(models$y),data=dat)
  nbCov_y <- ncol(cov_y) - 1
  
  if(hfun=="h.IP" | hfun=="h.EP1" | hfun=="h.EP1x")
    nbShapePar <- 1
  else
    nbShapePar <- 2
  
  if(hfun=="h.IP" | hfun=="h.EP1" | hfun=="h.EP2")
    nbScalePar <- 1
  else 
    nbScalePar <- 2
  
  if(nbShapePar==1) {
    shape <- b[1]
    covb <- matrix(shape,nrow=nbObs,ncol=1)
  } else {
    shape_x <- b[1]
    shape_y <- b[2]
    covb <- matrix(c(rep(shape_x,nbObs),rep(shape_y,nbObs)),ncol=2)
  }
  
  if(nbScalePar==1) {
    beta <- b[(nbShapePar+1):length(b)]
    scale <- cov_x%*%beta
    covb <- cbind(covb,exp(scale))
  } else {
    beta_x <- b[(nbShapePar+1):(nbShapePar+nbCov_x+1)]
    beta_y <- b[(nbShapePar+nbCov_x+2):length(b)]
    scale_x <- cov_x%*%beta_x
    scale_y <- cov_y%*%beta_y
    covb <- cbind(covb,exp(scale_x),exp(scale_y))
  }
  
  if(is.null(dat$id)) {
    # if single-observer case
    covb <- as.vector(t(covb))
  } else {
    # if double-observer case
    covb_obs1 <- covb[which(dat$id==1),]
    covb_obs2 <- covb[which(dat$id==2),]
    covb <- matrix(c(as.vector(t(covb_obs1)),as.vector(t(covb_obs2))),
                   ncol=2)
  }

  return(covb)
}


#' @title Constructs a HMM equivalent to a Poisson process.
#'
#' @description
#'  Returns list with $Pi, $pm and $delta of Poisson with same event rate as in the input HMM object 
#'  availhmm.
#'  
#' @param availhmm availability model list of the sort passed to \code{\link{est.hmltm}}.
#' @param zero REDUNDANT (I think - CHECK)
#' 
#' @export
poiss.equiv <- function(availhmm,zero=0){
  Pois.availhmm <- availhmm
  Pi <- availhmm$Pi
  pcu <- FALSE
  
  if(is.element("pcu",names(availhmm))) {
    names(availhmm)[which(names(availhmm)=="pcu")] <- "pm"
    pcu <- TRUE
  }
  
  pm <- availhmm$pm
  delta <- availhmm$delta
  if(is.vector(pm)&!is.matrix(Pi) | !is.vector(pm)&is.matrix(Pi)) 
    stop("Single animal: pcu/pm is not a vector or Pi is not a matrix")
  
  if(is.vector(pm)) { # convert to matrix and array so can use loop below
    Pi <- array(Pi,dim=c(2,2,1))
    pm <- matrix(pm,ncol=1)
    delta <- matrix(delta,ncol=1)
  }
  
  nw <- dim(pm)[2]
  PiPoiss <- matrix(c(zero,1-zero,zero,1-zero),byrow=TRUE,nrow=2)
  
  for(i in 1:nw) {
    Pi[,,i] <- PiPoiss
    delta[,i] <- compdelta(PiPoiss)
  }
  
  E <- Estate <- matrix(rep(0,2*nw),nrow=2)
  
  if(nw>1) {
    for(w in 1:nw) {
      Estate[,w] <- c(1/availhmm$Pi[1,2,w],1/availhmm$Pi[2,1,w])
      events <- apply(Estate*pm,2,sum)
      duration <- apply(Estate,2,sum)
      eventrate <- rep(0,length(duration))
      eventrate[duration>0] <- events[duration>0]/duration[duration>0]
      Pois.availhmm$Pi <- Pi
      
      if(pcu) 
        Pois.availhmm$pcu <- matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      else 
        Pois.availhmm$pm <- matrix(c(rep(0,nw),eventrate),nrow=2,byrow=TRUE)
      
      Pois.availhmm$delta <- delta
    }
  } else {
    Estate <- c(1/availhmm$Pi[1,2],1/availhmm$Pi[2,1])
    events <- sum(Estate*pm)
    duration <- sum(Estate)
    eventrate <- events/duration
    Pois.availhmm$Pi <- Pi[,,1]
    
    if(pcu) 
      Pois.availhmm$pcu <- c(0,eventrate)
    else 
      Pois.availhmm$pm <- c(0,eventrate)
    
    Pois.availhmm$delta <- as.vector(delta)
  }
  
  return(Pois.availhmm)
}

#' @title Logit function.
#'
#' @description
#'  Logit function
#'  
#' @param p probability (scalar or vector).
logit <- function(p) return(log(p/(1-p)))  # returns logit of p

#' @title Inverse logit function.
#'
#' @description
#'  Inverse logit function
#'  
#' @param x scalar or .
inv.logit <- function(x) return(1/(1+exp(-x))) # returns p from x=(logit of p)
