# Changed all pcu to pm 6th May 2013. Not yet tested with bootstrap.

#                        --------------------------------
#----------------------- Start detection hazard functions --------------------------------------
#                      (C++ code makes this all redundant)

#h.EP1=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.EP1(b)
#  return(par[1]*exp(-(abs(x)^par[2] + y^par[2])/(par[3]^par[2])))
#}
#tfm.EP1=function(par) return(c(logit(par[1]),log(par[2:3]))) # Transforms EP1 parameters from natural to link scale
#invtfm.EP1=function(b) return(c(inv.logit(b[1]),exp(b[2:3]))) # Transforms EP1 parameters from link to natural scale

#h.EP1.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #----------------------------------------------------------
#{
#  par=invtfm.EP1.0(b)
#  return(exp(-(abs(x)^par[1] + abs(y)^par[1])/(par[2]^par[1])))
#}
#tfm.EP1.0=function(par) return(log(par))
#invtfm.EP1.0=function(b) return(exp(b))

#h.EP1x.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #
  # Modified to have separate x- and y- scale functions
  #----------------------------------------------------------
#{
#  par=invtfm.EP1x.0(b)
#  return(exp(-((abs(x)/par[3])^par[1] + (abs(y)/par[2])^par[1])))
#}
#tfm.EP1x.0=function(par) return(log(par))
#invtfm.EP1x.0=function(b) return(exp(b))

#h.EP2=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.EP2(b)
#  return(par[1]*exp(-((abs(x)/par[4])^par[2] + (abs(y)/par[4])^par[3])))
#}
#tfm.EP2=function(par) return(c(logit(par[1]),log(par[2:4])))
#invtfm.EP2=function(b) return(c(inv.logit(b[1]),exp(b[2:4])))

#h.EP2.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #----------------------------------------------------------
#{
#  par=invtfm.EP2.0(b)
#  return(exp(-((abs(x)/par[3])^par[1] + (abs(y)/par[3])^par[2])))
#}
#tfm.EP2.0=function(par) return(log(par))
#invtfm.EP2.0=function(b) return(exp(b))

#h.EP2x.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999), but with p(0)=1
  #
  # Modified to have separate x- and y- scale functions
  #----------------------------------------------------------
#{
#  par=invtfm.EP2x.0(b)
#  return(exp(-((abs(x)/par[4])^par[1] + (abs(y)/par[3])^par[2])))
#}
#tfm.EP2x.0=function(par) return(log(par))
#invtfm.EP2x.0=function(b) return(exp(b))

#h.IP=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.IP(b)
#  return(par[1]*exp(par[2]*log(par[3])-(par[2]/2)*log(par[3]^2+x^2+y^2)))
#}
#tfm.IP=function(par) return(c(logit(par[1]),log(par[2:3])))
#invtfm.IP=function(b) return(c(inv.logit(b[1]),exp(b[2:3])))
                                                            
#h.IP.0=function(x,y,b)
  #----------------------------------------------------------
  # Detection hazard function prob(detect | available at x,y).
  # From Equation (10) of Skkaug & Schweder (1999).
  #----------------------------------------------------------
#{
#  par=invtfm.IP.0(b)
#  return(exp(par[1]*log(par[2])-(par[1]/2)*log(par[2]^2+x^2+y^2)))
#}
#tfm.IP.0=function(par) return(log(par))
#invtfm.IP.0=function(b) return(exp(b))

# Derivative stuff below from earlier version: now redundant
#dinvt.dpar=function(b,fun)
#{
#  switch(fun,
#         h.EP1=d.dp.EP1(b),
#         h.EP1.0=d.dp.EP1.0(b),
#         h.EP1x.0=d.dp.EP1x.0(b),
#         h.EP2=d.dp.EP2(b),
#         h.EP2.0=d.dp.EP2.0(b),
#         h.EP2x.0=d.dp.EP2x.0(b),
#         h.IP=d.dp.IP(b),
#         h.IP.0=d.dp.IP.0(b),
#  )
#}
#d.dp.EP1=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:3]))) # Returns derivative w.r.t p of link functions at parameter estimates (parameters on link scale)
#d.dp.EP1.0=function(b) return(exp(b))
#d.dp.EP1x.0=function(b) return(exp(b))
#d.dp.EP2=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:4])))
#d.dp.EP2.0=function(b) return(exp(b))
#d.dp.EP2x.0=function(b) return(exp(b))
#d.dp.IP=function(b) return(c(d.dx.inv.logit(b[1]),exp(b[2:3])))
#d.dp.IP.0=function(b) return(exp(b))

#d.dx.inv.logit=function(x) return(exp(-x)/(1+exp(-x))^2) # dervivative of inverse logit w.r.t. x=logit(p)

#----------------------- End detection hazard functions --------------------------------------
#                        ------------------------------
