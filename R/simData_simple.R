
#' Simulation (one observer)
#' 
#' Simulate detections in a hidden state line transect model with one observer.
#' 
#' @param nbObs Number of observations to simulate
#' @param pcu Vector of state dependent availability parameters (Pr(available|state))
#' @param gamma Transition probability matrix of the hidden state process
#' @param b Vector of parameters of the detection hazard function
#' @param hfun Name of the detection hazard function
#' @param xmax Maximum perpendicular distance in view
#' @param ymax Maximum forward distance in view
#' @param ystep Discretization step for the forward distance
#' 
#' @return A data frame of coordinates x and y of detected animals.
#' 
#' @examples
#' pcu <- c(0.9,0)
#' gamma <- matrix(c(0.9,0.1,0.1,0.9),nrow=2)
#' b <- c(0.9,0.75,10)
#' hfun <- "h.EP2.0"
#' data <- simData_simple(nbObs=100,pcu=pcu,gamma=gamma,b=b,hfun=hfun,xmax=80,ymax=80,ystep=1)
#' 
#' @export

simData_simple <- function(nbObs,pcu,gamma,b,hfun,xmax,ymax,ystep)
{
  # data frame of observations
  data <- data.frame(x=rep(NA,nbObs),y=rep(NA,nbObs))
  
  # initial distribution of the Markov chain
  delta <- solve(t(diag(2)-gamma+1),rep(1,2))
  
  obs <- 1 # observation index
  
  # while less than nbObs detections have occurred...
  while(obs<=nbObs) {
    t <- 1 # time index
    S <- sample(1:2,size=1,prob=delta) # state process
    A <- sample(0:1,size=1,prob=c(1-pcu[S],pcu[S])) # availability process
    x <- runif(1,0,xmax)
    ystart <- ymax-runif(1,0,ystep)
    # ystart <- ymax
    y <- seq(ystart,0,by=-ystep)
    
    if(A==1 & runif(1)<h_rcpp(x,y[t],b,hfun))
      D <- TRUE # detection process
    else
      D <- FALSE
    
    # while the animal is not detected and hasn't left the observed area
    while(!D & t<length(y)) {
      t <- t+1
      S <- sample(1:2,size=1,prob=gamma[S,])
      A <- sample(0:1,size=1,prob=c(1-pcu[S],pcu[S]))
      
      r <- runif(1)
      p <- pcu[S]*h_rcpp(x,y[t],b,hfun)
      
      if(A==1 & r<p)
        D <- TRUE
    }
    
    # if detected
    if(D) {
      data[obs,] <- c(x,y[t])
      obs <- obs+1 
    }
    # else, the animal has gone through the observed area without being detected
  }
  
  return(data)
}
