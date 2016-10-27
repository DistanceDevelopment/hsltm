#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Parameters transformation
//' 
//' Scale parameters from their "natural" interval to their "working" interval, 
//' the set of real numbers (for unconstrained optimization).
//' Written in C++.
//' 
//' @param par Parameters on their natural scale
//' @param hfun Name of hazard detection function
//' 
//' @return Unconstrained parameters (scaled to the set of real numbers)
//' 
//' @export
// [[Rcpp::export]]
arma::rowvec n2w_rcpp(arma::rowvec par, std::string hfun)
{
  int nbShapePar;
  if(hfun=="h.IP" || hfun=="h.EP1" || hfun=="h.EP1x")
    nbShapePar = 1;
  else
    nbShapePar = 2;
  
  arma::rowvec res(par.size());
  for(int i=0 ; i<nbShapePar ; i++)
    res[i] = log(par[i]); // shape > 0
  for(int i=nbShapePar ; i<par.size() ; i++)
    res[i] = par[i]; // betas real-valued
  
  return res;
}

//' Parameters inverse transformation
//' 
//' Scale parameters back from the set of real numbers ("working" interval)
//' to their "natural" interval.
//' Written in C++.
//' 
//' @param par Unconstrained parameters (scaled to the set of real numbers)
//' @param hfun Name of hazard detection function
//' 
//' @return Parameters on their natural scale
//' 
//' @export
// [[Rcpp::export]]
arma::rowvec w2n_rcpp(arma::rowvec par, std::string hfun)
{
  int nbShapePar;
  if(hfun=="h.IP" || hfun=="h.EP1" || hfun=="h.EP1x")
    nbShapePar = 1;
  else
    nbShapePar = 2;
  
  arma::rowvec res(par.size());
  for(int i=0 ; i<nbShapePar ; i++)
    res[i] = exp(par[i]); // shape > 0
  for(int i=nbShapePar ; i<par.size() ; i++)
    res[i] = par[i]; // betas real-valued
  
  return res;
}

//' Hazard detection function
//' 
//' Used within the C++ code, and in the simulation function simData.
//' 
//' @param x Perpendicular distance of the detection
//' @param y Forward distance of the detection
//' @param b Parameters of the hazard function
//' @param hfun Name of hazard detection function
//' 
//' @return The value of the hazard function at (x,y)
//' 
//' @export
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
double h_rcpp(double x, double y, arma::rowvec b, std::string hfun)
{
  double a1, a2; // intermediate variables
  double res; // result
  
  if(hfun=="h.EP1") {
    a1 = pow(std::abs(x),b[1])+pow(y,b[1]);
    a2 = pow(b[2],b[1]);
    res = b[0]*exp(-a1/a2);
  } 
  else if(hfun=="h.EP1x") {
    a1 = pow(std::abs(x)/b[2],b[0]);
    a2 = pow(std::abs(y)/b[1],b[0]);
    res = exp(-(a1+a2));
  }
  else if(hfun=="h.EP2x") {
    a1 = pow(std::abs(x)/b[3],b[0]);
    a2 = pow(std::abs(y)/b[2],b[1]);
    res = exp(-(a1+a2));
  }
  else if(hfun=="h.EP2") {
    a1 = pow(std::abs(x)/b[2],b[0]);
    a2 = pow(std::abs(y)/b[2],b[1]);
    res = exp(-(a1+a2));
  }
  else if(hfun=="h.IP") {
    a1 = b[0]*log(b[1]);
    a2 = b[1]/2*log(b[1]*b[1]+x*x+y*y);
    res = exp(a1-a2);
  }
  
  return res;
}

//' Forward distances in view
//' 
//' Computes all discrete forward distances in view. 
//' Written in C++.
//' 
//' @param x Perpendicular distance
//' @param null_yobs ???
//' @param yobs Forward distance
//' @param ymax Maximum forward distance in view
//' @param dy Forward distance step increment
//' 
//' @return A vector of the forward distances in view.
//' 
//' @export
// [[Rcpp::export]]
arma::vec gety_obs_rcpp(double x, bool null_yobs, double yobs, double ymax, double dy)
{
  double ymin = 0;
  
  if(!null_yobs)
    ymin = std::max(ymin,yobs);
  
  double yrange = ymax-ymin;
  double ndiv = floor(yrange/dy);
  
  if(ymax<=ymin) {
    arma::vec yi(1);
    yi[0] = ymin;
    return yi;
  }
  
  arma::vec yi(ndiv+1);
  for(int i=0 ; i<=ndiv ; i++)
    yi[i] = ymin+i*dy;
  
  return yi;
}

//' Probability of detection (single observer)
//'
//' Calculate probability of detection at (x,y), or if cdf=TRUE, by (x,y), in the case
//' of a model with one observer.
//' Written in C++.
//'  
//' @param x Perpendicular distance.
//' @param y Forward distance.
//' @param hfun Hazard function name.
//' @param b Hazard function parameter vector.
//' @param pcu Bernoulli state-dependent probability parameters.
//' @param Pi Markov model transition probability matrix.
//' @param delta Markov model stationary distribution.
//' @param ymax Maximum forward distance in view.
//' @param dy Forward distance step increment.
//' @param ally Flag for whether or not to return probabilities at all forward distances in view.
//' @param cdf Flag for whether or not to return cumulative distribution function in forward dimension.
//' This differs from specifying ally=TRUE in that ally=TRUE calculates the cdf from ymax 
//' to y=0, whereas cdf=TRUE calculates the cdf from ymax to y.
//' 
//' @export
// [[Rcpp::export]] 
arma::vec pxy_simple_rcpp(arma::vec x, arma::vec y, std::string hfun, arma::rowvec b, 
                          arma::vec pcu, arma::mat Pi, arma::rowvec delta, double ymax, 
                          double dy, bool ally, bool cdf)
{
  int nbStates = pcu.size();
  int nbObs = x.size();
  int nbPar = b.size()/nbObs;
  
  arma::vec p(nbObs); // result
  arma::rowvec bb(nbPar);
  arma::vec yi;
  arma::mat prodB(nbStates,nbStates);
  arma::mat Lambda(nbStates,nbStates);
  double hval;
  
  // loop over the observations
  for(int i=0 ; i<nbObs ; i++) {
    prodB.eye();
    
    for(int j=0 ; j<nbPar ; j++)
      bb[j] = b[i*nbPar+j]; // select appropriate row of b
    
    // transform parameters to natural scale
    bb = w2n_rcpp(bb,hfun);
    
    if(ally) // get all "y"s from min(y in view) to ymax
      yi = gety_obs_rcpp(x[i],true,0,ymax,dy);
    else // get "y"s from y[i] to ymax
      yi = gety_obs_rcpp(x[i],false,y[i],ymax,dy);
    
    int nT = yi.size();
    
    // probability of not detecting the animal until T_i-1
    for(int j=(nT-1) ; j>=1 ; j--) {
      hval = h_rcpp(x[i],yi[j],bb,hfun); // detection hazard function
      
      Lambda = arma::eye(nbStates,nbStates);
      Lambda.diag() = Lambda.diag() - hval*pcu; 
      
      prodB = prodB*Pi*Lambda;
    }
    
    hval = h_rcpp(x[i],yi[0],bb,hfun);
    
    if(cdf || ally) { // if not seen, calculation for last non-detection event
      Lambda = arma::eye(nbStates,nbStates);
      Lambda.diag() = Lambda.diag() - pcu*hval;
    }
    else { // if seen, calculation for last detection event
      Lambda.zeros();
      Lambda.diag() = pcu*hval;            
    }
    prodB = prodB*Pi*Lambda;
    
    p(i) = sum(delta*prodB);
  }
  
  if(cdf || ally) {
    // so far, p is actually Pr(not detected)
    for(int i=0 ; i<p.size() ; i++)
      p(i) = 1-p(i);
  }
  
  return p;
}

//' Probability of detection (double observer)
//'
//' Calculate probability of detection at (x,y), or if cdf=TRUE, by (x,y), in the case
//' of a model with two observers.
//' Written in C++.
//'  
//' @param x Matrix of perpendicular distances (two columns: one for each observer)
//' @param y Matrix of forward distances (two columns: one for each observer)
//' @param d Matrix of detections (two columns: one for each observer). d[i,j] is
//' 1 if observer j detected animal i, 0 otherwise.
//' @param hfun Hazard function name.
//' @param b Hazard function parameter matrix (two columns: one for each observer).
//' @param pcu Bernoulli state-dependent probability parameters.
//' @param Pi Markov model transition probability matrix.
//' @param delta Markov model stationary distribution.
//' @param ymax Maximum forward distance in view.
//' @param dy Forward distance step increment.
//' @param ally Flag for whether or not to return probabilities at all forward distances in view.
//' @param cdf Flag for whether or not to return cumulative distribution function in forward dimension.
//' This differs from specifying ally=TRUE in that ally=TRUE calculates the cdf from ymax 
//' to y=0, whereas cdf=TRUE calculates the cdf from ymax to y.
//' 
//' @export
// [[Rcpp::export]]
arma::vec pxy_double_rcpp(arma::mat x, arma::mat y, arma::mat d, std::string hfun, arma::mat b, 
                          arma::vec pcu, arma::mat Pi, arma::rowvec delta, double ymax, double dy,
                          bool ally, bool cdf)
{
  int nbStates = pcu.size();
  int nbObs = x.n_rows;
  int nbPar = b.n_rows/nbObs;
  
  arma::vec p(nbObs); // result
  arma::rowvec bb1(nbPar), bb2(nbPar);
  arma::vec yi;
  arma::mat prodB(nbStates,nbStates);
  arma::mat Lambda(nbStates,nbStates);
  double hval, hval1, hval2, h_either;
  double x1, y1, x2, y2;
  double pp;
  
  int cond1, cond2;
  int nT;
  
  // loop over the observations
  for(int i=0 ; i<nbObs ; i++) {
    
    // set conditions which are then used to distinguish the different possible cases
    if(d(i,0)==1 && d(i,1)==1)
      cond1 = 1; // if both observers detected the animal
    else
      cond1 = 2; // if only one observer detected the animal
    
    if(cond1==1 && y(i,0)==y(i,1))
      cond2 = 0; // if both observers detected the animal at the same time
    else if(cond1==1 && y(i,0)>=y(i,1))
      cond2 = 1; // if observer 1 detected the animal first
    else if(cond1==1 && y(i,0)<y(i,1))
      cond2 = 2; // if observer 2 detected the animal first
    else if(d(i,1)==0)
      cond2 = 1; // if observer 1 was the only one to detect the animal
    else if(d(i,0)==0)
      cond2 = 2; // if observer 2 was the only one to detect the animal
    
    // initialization
    prodB.eye();
    
    for(int j=0 ; j<nbPar ; j++) {
      bb1[j] = b(i*nbPar+j,0); // select parameters for observer 1
      bb2[j] = b(i*nbPar+j,1); // select parameters for observer 2
    }
    
    // transform parameters to natural scale
    bb1 = w2n_rcpp(bb1,hfun);
    bb2 = w2n_rcpp(bb2,hfun);
    
    //================================================//
    // 1. Compute probability up to first observation //
    //================================================//
    
    // coordinates of the first observation
    if(cond2==2) {
      x1 = x(i,1);
      y1 = y(i,1);
    }
    else {
      x1 = x(i,0);
      y1 = y(i,0);
    }
    
    yi = gety_obs_rcpp(x1,false,y1,ymax,dy);
    nT = yi.size();
    
    // probability of neither observer detecting the animal until first observation
    for(int j=(nT-1) ; j>=1 ; j--) {
      hval1 = h_rcpp(x1,yi[j],bb1,hfun); // detection function for observer 1
      hval2 = h_rcpp(x1,yi[j],bb2,hfun); // detection function for observer 2
      h_either = hval1+hval2-hval1*hval2; // detection function for either observer
      
      Lambda = arma::eye(nbStates,nbStates);
      Lambda.diag() = Lambda.diag() - h_either*pcu;
      
      prodB = prodB*Pi*Lambda;
    }
    
    hval1 = h_rcpp(x1,yi[0],bb1,hfun); // detection function for observer 1
    hval2 = h_rcpp(x1,yi[0],bb2,hfun); // detection function for observer 2
    h_either = hval1+hval2-hval1*hval2; // detection function for either observer
    
    // probability of either observer detecting the animal
    Lambda.zeros();
    Lambda.diag() = pcu*h_either;
    prodB = prodB*Pi*Lambda;
    
    //=================================================//
    // 2. Compute probability up to second observation //
    //=================================================//
    // if both observers detected the animal, but not at the same time
    if(cond1==1 && cond2!=0) {
      // coordinates of the second observation
      if(cond2==2) {
        x2 = x(i,0);
        y2 = y(i,0);
      }
      else {
        x2 = x(i,1);
        y2 = y(i,1);
      }
      
      yi = gety_obs_rcpp(x2,false,y2,y1,dy);
      nT = yi.size();
      
      // probability of other observer not detecting animal until second observation
      for(int j=(nT-1) ; j>=1 ; j--) {
        if(cond2==1)
          hval = h_rcpp(x2,yi[j],bb2,hfun); // detection function for second observer
        else
          hval = h_rcpp(x2,yi[j],bb1,hfun); // detection function for first observer
        
        Lambda = arma::eye(nbStates,nbStates);
        Lambda.diag() = Lambda.diag() - hval*pcu;
        
        prodB = prodB*Pi*Lambda;
      }
      
      if(cond2==1)
        hval = h_rcpp(x2,yi[0],bb2,hfun); // detection function for second observer
      else
        hval = h_rcpp(x2,yi[0],bb1,hfun); // detection function for first observer
      
      // probability of second observation
      Lambda.zeros();
      Lambda.diag() = pcu*hval;
      prodB = prodB*Pi*Lambda;
    }
    
    //=========================================//
    // 3. Probability of no second observation //
    //=========================================//
    // if only one observer detected the animal
    if(cond1==2) {
      yi = gety_obs_rcpp(x1,false,0,y1,dy);
      nT = yi.size();
      
      // probability of other observer not detecting animal at all
      for(int j=(nT-1) ; j>=1 ; j--) {
        if(cond2==1)
          hval = h_rcpp(x1,yi[j],bb2,hfun); // detection function for second observer
        else
          hval = h_rcpp(x1,yi[j],bb1,hfun); // detection function for first observer
        
        Lambda = arma::eye(nbStates,nbStates);
        Lambda.diag() = Lambda.diag() - hval*pcu;
        
        prodB = prodB*Pi*Lambda;
      }
    }
    
    // left-multiply by initial distribution, and sum
    p(i) = sum(delta*prodB);
    
    //====================================//
    // 4. Probability of observers' order //
    //====================================//
    // probability of both detecting the animal at the same time, given that someone detected it
    if(cond2==0) {
      x1 = x(i,0);
      y1 = y(i,0);
      hval1 = h_rcpp(x1,y1,bb1,hfun);
      hval2 = h_rcpp(x1,y1,bb2,hfun);
      h_either = hval1 + hval2 - hval1*hval2;
      p(i) = p(i)*hval1*hval2/h_either;
    }
    
    // probability that observer 1 detected the animal first, given that someone detected it
    if(cond2==1) {
      x1 = x(i,0);
      y1 = y(i,0);
      hval1 = h_rcpp(x1,y1,bb1,hfun);
      hval2 = h_rcpp(x1,y1,bb2,hfun);
      h_either = hval1 + hval2 - hval1*hval2;
      p(i) = p(i)*hval1*(1-hval2)/h_either;
    }
    
    // probability that observer 2 detected the animal first, given that someone detected it
    if(cond2==2) {
      x2 = x(i,1);
      y2 = y(i,1);
      hval1 = h_rcpp(x2,y2,bb1,hfun);
      hval2 = h_rcpp(x2,y2,bb2,hfun);
      h_either = hval1 + hval2 - hval1*hval2;
      p(i) = p(i)*(1-hval1)*hval2/h_either;
    }
  }
  
  if(cdf || ally) {
    // so far, p is actually Pr(not detected)
    for(int i=0 ; i<p.size() ; i++)
      p(i) = 1-p(i);
  }
  
  return p;
}
