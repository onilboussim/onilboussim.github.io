
# R scripts for implementing calculations in
# Hahn, Hirano, and Karlan, "Adaptive Experimental Design Using the Propensity Score"
# Last updated: 9 July 2009
#
# Please send any bug reports or suggestions to Keisuke Hirano (hirano@u.arizona.edu)


datasum <- function(y,x,t) {
  # summarize data
  # assumes X is discrete

  # calculate unique values of X
  uniqx <- sort(unique(x))

  # calculate probability mass function of X
  calcfx <- function(c,x) mean(x==c)
  fx <- sapply(uniqx, calcfx, x)

  # calculate conditional variances
  calcvar <- function(c,y,x) var(y[x==c])
  v0 <- sapply(uniqx, calcvar, y[t==0], x[t==0])
  v1 <- sapply(uniqx, calcvar, y[t==1], x[t==1])

  # calculate conditional means
  calcm0 <- function(c,y,x,t){ mean(y[x==c & t==0]) }
  calcm1 <- function(c,y,x,t){ mean(y[x==c & t==1]) }  
  m0 <- sapply(uniqx, calcm0, y,x,t)
  m1 <- sapply(uniqx, calcm1, y,x,t)
  
  # calculate conditional treatment effects and ATE
  # calcbeta <- function(c,y,x,t) {mean(y[x==c & t==1]) - mean(y[x==c & t==0])}
  # betavec <- sapply(uniqx, calcbeta, y, x, t)
  betavec <- m1 - m0
  ATE <- sum(fx*betavec)

  # report results
  list(uniqx=uniqx, fx=fx, v0=v0, m0 = m0, v1=v1, m1=m1, beta=betavec, ATE=ATE)
}

avar <- function(p,fx,v0,v1,beta){
  # calculate asymptotic variance of estimators
  # p = vector of treatment probabilities
  # fx = probability mass vector for X
  # v0 = vector of conditional variance given x, t=0
  # v1 = vector of conditonal variance given x, t=1
  # beta = vector of conditional average treatment effects
  tau <- sum(fx*beta)
  sum( fx* ((v1/p)+(v0/(1-p))+(beta-tau)^2) )
}

ssprob <- function(popt,p1,kappa=.5){
  # calculate optimal second stage probabilities
  # popt = vector of optimal overall probs
  # p1 = 1st stage probs
  # kappa = fraction in 1st stage
  (1/(1-kappa))*(popt - kappa*p1)
}

hhk.calcpunc <- function(v0,v1){
  # calculate unconstrained optimal probabilities
  # v0, v1 are vectors of conditional variances
  # returns vector of treatment assignment probabilities
  s1 <- sqrt(v1)
  s0 <- sqrt(v0)
  s1/(s0+s1)  
}

hhk.unc <- function(x,y,t,treatmentprob=.5){
  # input data and calculate optimal unconstrained treatment probabilities
  # treatmentprob is the "original" treatment probabilities, used for variance comparison
  # calls hhk.calcpunc routine
  
  hhkpars <- datasum(y=y,x=x,t=t)
  nonoptvb <- avar(p=rep(treatmentprob,length(hhkpars$uniqx)), fx=hhkpars$fx, v0=hhkpars$v0, v1=hhkpars$v1, beta=hhkpars$beta)
  
  optimalp <- hhk.calcpunc(v0=hhkpars$v0, v1=hhkpars$v1)
  optimalvb <- avar(p=optimalp, fx=hhkpars$fx, v0=hhkpars$v0, v1=hhkpars$v1, beta=hhkpars$beta)
  list(popt=optimalp, uniqx=hhkpars$uniqx, fx=hhkpars$fx, v0=hhkpars$v0, m0=hhkpars$m0, v1=hhkpars$v1, m1=hhkpars$m1, beta=hhkpars$beta, ATE=hhkpars$ATE, vb=optimalvb, nonoptvb = nonoptvb)
}

Vcomp <- function(AVopt,AVnon,n){
  # compare optimal and nonoptimal variances
  # AVopt = optimal variance
  # AVnon = nonoptimal variance
  # n = sample size associated with AVnon
  # note: assumes AVnon > AVopt without checking
  # uses uniroot to numerically find equivalent sample size 
  
  vreduction <- 1-(AVopt/AVnon)
  print(paste(100*vreduction, "% reduction in variance.",sep=""),quote=FALSE)

  # difference in variance 
  Vdiff <- function(nopt,AVopt,n,AVnon){
    (AVopt/nopt)-(AVnon/n)
  }

  # use uniroot to find the root of Vdiff 
  nopt <- uniroot(Vdiff,AVopt=AVopt,n=n,AVnon=AVnon,lower=1,upper=n+1)

  # report results
  print(paste("Equivalent sample size: ", ceiling(nopt$root),sep=""),quote=FALSE)
  list(vreduction=vreduction,nopt=nopt$root)
}


hhk.calcpcon <- function(v0,v1,fx,pbar){
  # core routine to calculate optimal treatment probabilities
  # subject to constraint on overall prob
  #
  # v0,v1 are vectors of conditional variances
  # fx is vector of cell probabilities
  # pbar is overall treatment probability
  # calls the R routine constrOptim after substituting for one of the probabilities
  
  # objective function to minimize 
  constrObjective <- function(shortP,pbar,fx,v0,v1){
    K <- length(fx)   # should check if length(shortp) = K-1
    fK <- fx[K]
    shortFx <- fx[1:(K-1)]
    pK <- (pbar - sum(shortFx*shortP))/fK
    longP <- c(shortP,pK)
    sum( fx* ((v1/longP)+(v0/(1-longP))) )
  }

  # gradient of objective function (useful for constrOptim)
  constrGrad <- function(shortP,pbar,fx,v0,v1){
    K <- length(fx)
    fK <- fx[K]
    shortFx <- fx[1:(K-1)]
    pK <- (pbar - sum(shortFx*shortP))/fK
    shortV0 <- v0[1:(K-1)]
    shortV1 <- v1[1:(K-1)]

    Kterm <- (v1[K]/(pK^2)) - (v0[K]/( (1-pK)^2 ))    
    gradient <- shortFx*( -shortV1/(shortP^2) + (shortV0/((1-shortP)^2) ) + Kterm) 
  }

  K <- length(fx)   # should check if length(shortp) = K-1
  fK <- fx[K]
  shortFx <- fx[1:(K-1)]

  # constraints on probabilities: ui %*% shortP >= ci
  # we want all probabilities between 0 and 1 
  ui <- rbind(diag(K-1),-diag(K-1),-shortFx,shortFx)
  ci <- c( rep(0,times=K-1), rep(-1,times=K-1), -pbar, pbar-fK)

  shortPi <- rep(pbar,times=K-1)    #initial prob vector
  
  minvar <- constrOptim(theta=shortPi, f=constrObjective, grad=constrGrad, ui=ui,ci=ci, pbar=pbar,fx=fx,v0=v0,v1=v1)

  # after solving for optimum, add pK to the probability vector
  c(minvar$par, (pbar - sum(shortFx*minvar$par))/fK)
}

hhk.con <- function(x,y,t,pbar=.5,kappa=.5){
  # Given data, calculate optimal treatment probs subject to constraint on overall treatment prob
  # x = discretized covariate
  # y = outcome
  # t = treatment
  # pbar = overall treatment probability
  # kappa = sample fraction in first stage (this is not currently used)
  # calls the routine hhk.calcpcon

  # first summarize data and calculate variance under pure randomization
  hhkpars <- datasum(y=y,x=x,t=t)
  nonoptvb <- avar(p=rep(pbar,length(hhkpars$uniqx)), fx=hhkpars$fx, v0=hhkpars$v0, v1=hhkpars$v1, beta=hhkpars$beta)

  # calculate optimal treatment probabilities and calculate associated variance
  optimalp <- hhk.calcpcon(v0=hhkpars$v0,v1=hhkpars$v1,fx=hhkpars$fx,pbar=pbar)
  optimalvb <- avar(p=optimalp, fx=hhkpars$fx, v0=hhkpars$v0, v1=hhkpars$v1, beta=hhkpars$beta)

  # report results
  list(popt=optimalp, uniqx=hhkpars$uniqx, fx=hhkpars$fx, v0=hhkpars$v0, m0=hhkpars$m0, v1=hhkpars$v1, m1=hhkpars$m1, beta=hhkpars$beta, ATE=hhkpars$ATE, vb=optimalvb, nonoptvb=nonoptvb)
}

maketable <- function(results,x,t){
  # construct table of results in the form appearing in HHK paper
  # results = output from either hhk.unc or hhk.con
  
  count0 <- function(c,x,t){ sum( x==c & t==0 ) }
  count1 <- function(c,x,t){ sum( x==c & t==1 ) }
  n0 <- sapply(results$uniqx,count0,x,t)
  n1 <- sapply(results$uniqx,count1,x,t)
  data.frame(uniqx=results$uniqx,n0,m0=results$m0,v0=results$v0,n1,m1=results$m1,v1=results$v1,popt=results$popt)
}

