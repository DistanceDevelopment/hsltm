
#' Simulation (one observer)
#' 
#' Simulate detections in a hidden state line transect model with one observer.
#' 
#' @param nbAnimals Number of animals
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
#' data <- simData_simple(nbAnimals=100,pcu=pcu,gamma=gamma,b=b,hfun=hfun,xmax=80,ymax=80,ystep=1)
#' 
#' @export

simData_simple <- function(nbAnimals,pcu,gamma,b,hfun,xmax,ymax,ystep)
{
  # data frame of observations
  data <- data.frame(x=numeric(),y=numeric())
  
  # initial distribution of the Markov chain
  delta <- solve(t(diag(2)-gamma+1),rep(1,2))
  
  # loop over the animals
  for(zoo in 1:nbAnimals) {
    t <- 1 # time index
    S <- sample(1:2,size=1,prob=delta) # state process
    A <- sample(0:1,size=1,prob=c(1-pcu[S],pcu[S])) # availability process
    x <- runif(1,0,xmax)
    # ystart <- ymax+runif(1,0,ystep)
    ystart <- ymax
    y <- seq(ystart,0,by=-ystep) # grid of y values
    
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
    if(D)
      data <- rbind(data,c(x,y[t]))
    
    # else, the animal has gone through the observed area without being detected
  }
  
  colnames(data) <- c("x","y")
  return(data)
}


#' Simulation (two observers)
#' 
#' Simulate detections in a hidden state line transect model with two observers.
#' 
#' @param nbAnimals Number of animals
#' @param pcu Vector of state dependent availability parameters (Pr(available|state))
#' @param gamma Transition probability matrix of the hidden state process
#' @param b List of vectors of parameters of the detection hazard function for each
#' observer
#' @param hfun Name of the detection hazard function (or vector of two names, if the
#' two observers have different detection functions)
#' @param xmax Maximum perpendicular distance in view
#' @param ymax Maximum forward distance in view
#' @param ystep Discretization step for the forward distance
#' 
#' @return A data frame with the following columns:
#' \item{id}{Index of the observer}
#' \item{d}{1 if observed, 0 otherwise}
#' \item{x}{Perpendicular distance of the observation}
#' \item{y}{Forward distance of the observation}
#' 
#' @details If xmax and ymax are not specified, they are computed to be the points
#' in which h(xmax,0) (resp. h(0,ymax)) is approximately 0.01*h(0,0). If ystep is
#' not specified, the y-dimension is discretized to 100 intervals.
#' 
#' @examples
#' pcu <- c(0.9,0)
#' gamma <- matrix(c(0.9,0.1,0.1,0.9),nrow=2)
#' b <- list(c(0.9,0.75,10),c(0.5,0.8,8))
#' hfun <- "h.EP2.0"
#' data <- simData_double(nbAnimals=100,pcu=pcu,gamma=gamma,b=b,hfun=hfun,xmax=100,ymax=100,ystep=1)
#' 
#' @export

simData_double <- function(nbAnimals,pcu,gamma,b,hfun,xmax,ymax,ystep)
{
  if(length(hfun)==1)
    hfun <- c(hfun,hfun)
  
  # data frame of observations
  data <- data.frame(id=numeric(),d=numeric(),x=numeric(),y=numeric())
  
  # initial distribution of the Markov chain
  delta <- solve(t(diag(2)-gamma+1),rep(1,2))
  
  # loop over animals
  for(zoo in 1:nbAnimals) {
    t <- 1 # time index
    k <- 1 # index used to indicate which row of the data should be filled next
    S <- sample(1:2,size=1,prob=delta) # state process
    A <- sample(0:1,size=1,prob=c(1-pcu[S],pcu[S])) # availability process
    D <- c(FALSE,FALSE) # detection process
    x <- runif(1,0,xmax)
    ystart <- ymax+runif(1,0,ystep)
    y <- seq(ystart,0,by=-ystep) # grid of y values
    
    if(A==1 & runif(1)<h_rcpp(x,y[t],b[[1]],hfun[1])) {
      D[1] <- TRUE
      data <- rbind(data,c(1,1,x,y[t]))
      k <- k+1
    }
    if(A==1 & runif(1)<h_rcpp(x,y[t],b[[2]],hfun[2])) {
      D[2] <- TRUE
      data <- rbind(data,c(2,1,x,y[t]))
      k <- k+1
    }
    
    # while the animal is not detected by both observers and hasn't left the observed area
    while(length(which(D))<2 & t<length(y)) {
      t <- t+1
      S <- sample(1:2,size=1,prob=gamma[S,])
      A <- sample(0:1,size=1,prob=c(1-pcu[S],pcu[S]))
      
      if(A==1 & !D[1] & runif(1)<h_rcpp(x,y[t],b[[1]],hfun[1])) {
        D[1] <- TRUE
        data <- rbind(data,c(1,1,x,y[t]))
        k <- k+1
      }
      
      if(A==1 & !D[2] & runif(1)<h_rcpp(x,y[t],b[[2]],hfun[2])) {
        D[2] <- TRUE
        data <- rbind(data,c(2,1,x,y[t]))
        k <- k+1
      }
    }
    
    # if not detected by one observer, add corresponding row
    if(D[1] & !D[2])
      data <- rbind(data,c(2,0,NA,NA))
    if(D[2] & !D[1])
      data <- rbind(data,c(1,0,NA,NA))
  }
  
  colnames(data) <- c("id","d","x","y")
  return(data)
}

