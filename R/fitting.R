#Fitting functions.

fisherscore.nhm <- function(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=TRUE,hessian,linesearch=FALSE,damped=FALSE,damppar=1,rtol,atol,score_test=FALSE,intens,gensolve,ncores=1) {
  #Fisher scoring algorithm

  f0 <- genlikelihood.nhm(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=TRUE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
  der <- attr(f0,"gradient")
  fish <- attr(f0,"fisher")
  if (score_test) return(list(g=-der,I=fish))
  if (damped) fish <- fish + diag(npar,damppar)
  if (min(eigen(fish)$values)<1e-5) stop("Fisher information matrix singular or near singular: Try a damped version of the scoring algorithm or reduce the number of parameters.")
  print(c(f0,der))
  g <- sum(abs(der))
  ldiff <- 100
  lold <- f0
  while (g>1e-3 & ldiff>1e-6) { #Stops if sum of absolute gradients smaller than 1e-3 or stops improving
    #May want different tolerance for models with many parameters.
    dir <- solve(fish)%*%der
    if (linesearch) {
      step <- optimize(linesearch.nhm,c(0,1),dir=dir,x0=x0,npar=npar,ncov=ncov,nstate=nstate,ncovvals=ncovvals,covs=covs,covtable=covtable,timesets=timesets,fromtimeindex=fromtimeindex,totimeindex=totimeindex,fromtimeset=fromtimesets,totimesets=totimesets,fromsets=fromsets,tosets=tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tcrit,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
      print(c(step$value,step$minimum))
      x1<- x0 - dir*step$minimum
    }else{
      x1 <-  x0 - dir
    }
    xold <- x0
    x0 <- x1
    f0 <- genlikelihood.nhm(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=TRUE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
    der <- attr(f0,"gradient")
    fish <- attr(f0,"fisher")
    if (damped) fish <- fish + diag(npar,damppar)
    print(c(f0,der))
    gold <- g
    g <- sum(abs(der))
    ldiff <- lold-f0
    if (ldiff>0) {
      lold <- f0
    }else{ #Revert to old values if function value has got worse
      x0 <- xold
      f0 <- lold
      der <- attr(lold,"gradient")
      fish <- attr(lold,"fisher")
      if (damped) fish <- fish + diag(npar,1)
    }
  }
  if (g>5e-1) warning("Final gradient large: Alter starting values or try a more robust optimization algorithm")
  if (hessian) { #Compute observed Fisher information
    print("Computing Hessian by finite-differences")
    hess <- finitedifferences.nhm(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,rtol,atol,intens,gensolve,ncores)
    if (min(Re(eigen(hess)$values))<0) warning("Hessian not positive definite")
  }else{
    hess <- NULL #Use the expected Fisher information
  }
  if (damped) fish<-fish - diag(npar,damppar)
  list(par=x0,value=f0[[1]],gradient=der,fisher_info=fish,hess=hess)
}


linesearch.nhm<-function(x,dir,x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,rtol,atol,intens,gensolve,ncores) {
  #Performs a line search for best step length for Fisher scoring algorithm
  a<-genlikelihood.nhm(x0 - x*dir,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
  print(a)
  a
}


#Parallel version of finitedifferences.nhm

finitedifferences.nhm <- function(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,rtol,atol,intens,gensolve,ncores,fixedpar=NULL,fixval=NULL) {
  #Computes the Hessian at the final parameter value found by the Fisher scoring algorithm
  #Uses a central differences approximation requiring 2*npar gradient evaluations
  
  #Need to add an additional control parameter to disable the parallelization here if needed.
  parallel_hess <- attr(intens,"parallel_hess")
  if (parallel_hess) {
    ncores_p <- ncores
    ncores_l <- 1
  }else{
    ncores_p <- 1
    ncores_l <- ncores
  }
  
  
  if (!is.null(fixedpar)) {
    y <- rep(0,npar)
    nonfixed <- (1:npar)[-fixedpar]
    y[nonfixed]<-x0
    y[fixedpar]<-fixval
  }else{
    nonfixed<-1:npar
  }
  
  #Should really skip any with fixed parameter values...
  stepsize <- 1e-5
  hess <- matrix(0,npar,npar)
  up <- down <- matrix(0,npar,npar)
  
  uplist <- mclapply(nonfixed, function(k) { step <- rep(0,npar);step[k] <- stepsize; attr(genlikelihood.nhm(x0+step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores_l),"gradient")},mc.cores=ncores_p)
  downlist <- mclapply(nonfixed, function(k) { step <- rep(0,npar);step[k] <- stepsize; attr(genlikelihood.nhm(x0-step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores_l),"gradient")},mc.cores=ncores_p)
  
  for (l in 1:length(nonfixed)) {
    k <- nonfixed[l]
    up[k,] <- uplist[[l]]
    down[k,] <- downlist[[l]]
    hess[k,] <-  c((up[k,] - down[k,])/(4*stepsize)) + c((up[,k] - down[,k])/(4*stepsize)) #Exploit that hessian should be symmetric
  }
  if (!is.null(fixedpar)) hess <- hess[-fixedpar,-fixedpar]
  return(hess)
}


#Version for misclassification models
finitedifferences_misc.nhm <- function(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar=fixedpar,fixval=fixval) {
  #Computes the Hessian at the final parameter value found by the Fisher scoring algorithm
  #Uses a central differences approximation requiring 2*npar gradient evaluations
  
  #Need to add an additional control parameter to disable the parallelization here if needed.
  parallel_hess <- attr(intens,"parallel_hess")
  if (parallel_hess) {
    ncores_p <- ncores
    ncores_l <- 1
  }else{
    ncores_p <- 1
    ncores_l <- ncores
  }
  
  
  if (!is.null(fixedpar)) {
    y <- rep(0,npar)
    nonfixed <- (1:npar)[-fixedpar]
    y[nonfixed]<-x0
    y[fixedpar]<-fixval
  }else{
    nonfixed<-1:npar
  }
  
  
  stepsize <- 1e-5
  hess <- matrix(0,npar,npar)
  up <- down <- matrix(0,npar,npar)
  
  uplist <- mclapply(nonfixed, function(k) { step <- rep(0,npar);step[k] <- stepsize; -attr(genlikelihood_misc.nhm(x0+step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,obslist,state,subject,time,cov2,init_state,indout=FALSE,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores=ncores_l,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ),"gradient")},mc.cores=ncores_p)
  downlist <- mclapply(nonfixed, function(k) { step <- rep(0,npar);step[k] <- stepsize; -attr(genlikelihood_misc.nhm(x0-step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,obslist,state,subject,time,cov2,init_state,indout=FALSE,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores=ncores_l,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ),"gradient")},mc.cores=ncores_p)
  
  for (l in 1:length(nonfixed)) {
    k <- nonfixed[l]
    up[k,] <- uplist[[l]]
    down[k,] <- downlist[[l]]
    hess[k,] <-  c((up[k,] - down[k,])/(4*stepsize)) + c((up[,k] - down[,k])/(4*stepsize)) #Exploit that hessian should be symmetric
  }
  if (!is.null(fixedpar)) hess <- hess[-fixedpar,-fixedpar]
  return(hess)
}



# 
# 
# finitedifferences.nhm <- function(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,rtol,atol,intens,gensolve,ncores,fixedpar=NULL,fixval=NULL) {
#   #Computes the Hessian at the final parameter value found by the Fisher scoring algorithm
#   #Uses a central differences approximation requiring 2*npar gradient evaluations
# 
# 
#   if (!is.null(fixedpar)) {
#     y <- rep(0,npar)
#     nonfixed <- (1:npar)[-fixedpar]
#     y[nonfixed]<-x0
#     y[fixedpar]<-fixval
#   }else{
#     nonfixed<-1:npar
#   }
# 
#   #Should really skip any with fixed parameter values...
#   stepsize <- 1e-5
#   hess <- matrix(0,npar,npar)
#   up <- down <- matrix(0,npar,npar)
#   for (k in nonfixed) {
#     step <- rep(0,npar)
#     step[k] <- stepsize
#     f1 <- genlikelihood.nhm(x0+step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
#     up[k,] <- attr(f1,"gradient")
#     print(up[k,])
#     f2 <- genlikelihood.nhm(x0-step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
#     down[k,] <- attr(f2,"gradient")
#     print(down[k,])
#   }
#   for (k in nonfixed) {
#     hess[k,] <-  c((up[k,] - down[k,])/(4*stepsize)) + c((up[,k] - down[,k])/(4*stepsize)) #Exploit that hessian should be symmetric
#   }
#   if (!is.null(fixedpar)) hess <- hess[-fixedpar,-fixedpar]
#   return(hess)
# }
# 
# #Version for misclassification models
# finitedifferences_misc.nhm <- function(x0,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar=fixedpar,fixval=fixval) {
#   #Computes the Hessian at the final parameter value found by the Fisher scoring algorithm
#   #Uses a central differences approximation requiring 2*npar gradient evaluations
# 
#   if (!is.null(fixedpar)) {
#     y <- rep(0,npar)
#     nonfixed <- (1:npar)[-fixedpar]
#     y[nonfixed]<-x0
#     y[fixedpar]<-fixval
#   }else{
#     nonfixed<-1:npar
#   }
# 
# 
#   stepsize <- 1e-5
#   hess <- matrix(0,npar,npar)
#   up <- down <- matrix(0,npar,npar)
#   for (k in 1:npar) {
#     step <- rep(0,npar)
#     step[k] <- stepsize
#     f1 <- genlikelihood_misc.nhm(x0+step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,obslist,state,subject,time,cov2,init_state,indout=FALSE,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ)
#     up[k,] <- -attr(f1,"gradient")
#     print(up[k,])
#     f2 <- genlikelihood_misc.nhm(x0-step,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore=FALSE,obslist,state,subject,time,cov2,init_state,indout=FALSE,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ)
#     down[k,] <- -attr(f2,"gradient")
#     print(down[k,])
#   }
#   for (k in nonfixed) {
#     hess[k,] <-  c((up[k,] - down[k,])/(4*stepsize)) + c((up[,k] - down[,k])/(4*stepsize)) #Exploit symmetry.
#   }
#   if (!is.null(fixedpar)) hess <- hess[-fixedpar,-fixedpar]
#   return(hess)
# }


#Not sure why the signs are different for misc versus non-misc.
# bhhh.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,initp_nhm=NULL,nparQ=NULL) {
#   #Wrapper function for maxLik: requires log-likelihood rather than negative log-likelihood and individual level contributions.
#
#   #out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,use.deriv,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve)
#   if (is.null(emat_nhm)) {
#     out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,ncores=ncores)
#     val <- -attr(out,"indlik")
#     attr(val,"gradient") <- -attr(out,"indder")
#   }else{
#     out <- genlikelihood_misc.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,initp_nhm,nparQ)
#     val <- attr(out,"indlik")
#     attr(val,"gradient") <- attr(out,"indder")
#   }
#   print(sum(val))
#   return(val)
# }

#Version to allow the possibility of fixed values.
bhhh.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL) {
  #Wrapper function for maxLik: requires log-likelihood rather than negative log-likelihood and individual level contributions.
  if (!is.null(fixedpar)) {
   y <- rep(0,npar)
   y[-fixedpar] <- x
   y[fixedpar] <- fixval
   x <- y
  }
  #out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,use.deriv,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve)
  if (is.null(emat_nhm)) {
    out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve,ncores=ncores)
    val <- -attr(out,"indlik")
    grad <- -attr(out,"indder")
    if (!is.null(fixedpar)) {
      grad <- grad[,-fixedpar]
    }
    attr(val,"gradient") <- grad
  }else{
    out <- genlikelihood_misc.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ)
    val <- attr(out,"indlik")
    grad <- attr(out,"indder")
    if (!is.null(fixedpar)) {
      grad <- grad[,-fixedpar]
    }
    attr(val,"gradient") <- grad
  }
  #print(sum(val))
  return(val)
}


nlminb_obj.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL,opt_env) {
  cur_config <- get("opt_config", envir = opt_env)
  if (!identical(cur_config$par, x)) {
    out <- nlminb.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar,fixval)
    assign("opt_config", out, envir=opt_env)
    return(out$obj)
  }else{
    return(cur_config$obj)
  }
}

nlminb_grad.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL,opt_env) {
  cur_config <- get("opt_config", envir = opt_env)
  if (!identical(cur_config$par, x)) {
    out <- nlminb.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar,fixval)
    assign("opt_config", out, envir=opt_env)
    return(out$grad)
  }else{
    return(cur_config$grad)
  }
}


nlminb_hess.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL,opt_env) {
  cur_config <- get("opt_config", envir = opt_env)
  if (!identical(cur_config$par, x)) {
    out <- nlminb.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar,fixval)
    assign("opt_config", out, envir = opt_env)
    return(out$hess)
  }else{
    return(cur_config$hess)
  }
}


#
# nlminb_grad.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL) {
#   nlminb.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar,fixval,otype=2)
# }
#
# nlminb_hess.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=TRUE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL) {
#   nlminb.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ,fixedpar,fixval,otype=3)
# }


nlminb.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist=NULL,state=NULL,subject=NULL,time=NULL,cov2=NULL,init_state=NULL,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm=NULL,ncores,firstobs,initp=NULL,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ=NULL,fixedpar=NULL,fixval=NULL) {
  if (!is.null(fixedpar)) {
    y <- rep(0,npar)
    y[-fixedpar] <- x
    y[fixedpar] <- fixval
    x <- y
  }
  indout<-TRUE
  #out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,use.deriv,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve)
  if (is.null(emat_nhm)) {
    out <- genlikelihood.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,indout,sublist,nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve,ncores=ncores)
    lik <- out[[1]]
    grad <- attr(out,"gradient")
    igrad <- attr(out,"indder")
    hess <- apply(array(apply(igrad,1,function(x) outer(x,x)),c(npar,npar,dim(igrad)[1])),c(1,2),sum)
    if (!is.null(fixedpar)) {
        grad <- grad[-fixedpar]
        hess <- hess[-fixedpar,-fixedpar]
    }
  }else{
    out <- genlikelihood_misc.nhm(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,obslist,state,subject,time,cov2,init_state,indout,sublist,nsub,rtol,atol,intens,gensolve,emat_nhm,ncores,firstobs,initp,initcovs,finalcovs,initp_nhm,nparQ)
    lik <- -out[[1]]
    grad <- -attr(out,"gradient")
    igrad <- -attr(out,"indder")
    hess <- apply(array(apply(igrad,1,function(x) outer(x,x)),c(npar,npar,dim(igrad)[1])),c(1,2),sum)
    if (!is.null(fixedpar)) {
      grad <- grad[-fixedpar]
      hess <- hess[-fixedpar,-fixedpar]
    }
    }
  return(list(par=x, obj=lik, grad=grad, hess=hess))
}


genlikelihood.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,finalcovs=NULL,gensolve,ncores) {
  #Main likelihood function
  if (is.null(censor)) censor <- 0
  cdlik <- array(0,c(ncovvals,npar))
  if (fishscore) cd2lik <- array(0,c(ncovvals,npar,npar))
  clik <- 0
  if (indout) {
    indlik <- rep(0,nsub)
    indder <- array(0,c(nsub,npar))
  }
  par <- x
  if (ncov==0) covvalue <- NULL
  if (ncores==1) {
    genout <- lapply(1:ncovvals,genlikelihoodwrapper.nhm, par=par,npar=npar,ncov=ncov,nstate=nstate,ncovvals=ncovvals,covs=covs,covtable=covtable,timesets=timesets,fromtimeindex=fromtimeindex,totimeindex=totimeindex,fromtimesets=fromtimesets,totimesets=totimesets,fromsets=fromsets,tosets=tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tcrit,fishscore=fishscore,indout=indout,sublist=sublist,nsub=nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve)
  }else{
    genout <- mclapply(1:ncovvals,genlikelihoodwrapper.nhm, mc.cores = ncores,par=par,npar=npar,ncov=ncov,nstate=nstate,ncovvals=ncovvals,covs=covs,covtable=covtable,timesets=timesets,fromtimeindex=fromtimeindex,totimeindex=totimeindex,fromtimesets=fromtimesets,totimesets=totimesets,fromsets=fromsets,tosets=tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tcrit,fishscore=fishscore,indout=indout,sublist=sublist,nsub=nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve)
  }

  #Retrieve the relevant parts
  clik <- sapply(genout,function(x) x$lval)
  cdlik <- t(sapply(genout,function(x) x$dlval))
  if (fishscore) cd2lik <- aperm(array(sapply(genout,function(x) x$d2lval),c(npar,npar,ncovvals)),c(3,1,2))
  if (indout) {
    indlikV <- sapply(genout,function(x) x$indlik)
    indlik <- apply(indlikV,1,sum)
    indderV <- array(sapply(genout,function(x) x$indder),c(nsub,npar,ncovvals))
    indder <- apply(indderV,c(1,2),sum)
  }
  out <- sum(clik)
  attr(out,"gradient") <- apply(cdlik,2,sum)
  if (fishscore) attr(out,"fisher") <- apply(cd2lik,c(2,3),sum)
  if (indout) {
    attr(out,"indlik") <- indlik
    attr(out,"indder") <- indder
  }
  out
}


genlikelihoodwrapper.nhm <- function(j, par,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit,fishscore,indout,sublist,nsub=NULL,rtol,atol,intens,finalcovs,gensolve,ncores) {
  #Cycles through the set of unique covariate patterns
  if (ncov>0) covvalue <- unlist(covs[j,])
  nevent <- covtable[j]
  times <- timesets[[j]]
  low <- fromtimeindex[[j]]
  up <- totimeindex[[j]]
  lowt <- fromtimesets[[j]]
  upt <- totimesets[[j]]
  from <- fromsets[[j]]
  to <- tosets[[j]]
  subs <- sublist[[j]]
  likcal <- sublik_deriv.nhm(nevent,times,low,up,lowt,upt,from,to,subs,npar,ncov,nstate,covvalue,par,tcrit,censor,censor.states,death,death.states,fishscore,rtol,atol,transitions=FALSE,intens,finalcovs,gensolve)
  lval <- sum(likcal$lik)
  dlval <- apply(likcal$dlik,2,sum)
  if (fishscore) {
    d2lval <- array(apply(likcal$d2lik,c(2,3),sum),c(npar,npar))
  }else{
    d2lval <- NULL
  }
  if (indout) {
    indlik <- rep(0,nsub)
    indder <- array(0,c(nsub,npar))
    for (k in 1:length(subs)) {
      indlik[subs[k]] <- indlik[subs[k]]+likcal$lik[k]
      indder[subs[k],] <- indder[subs[k],]+likcal$dlik[k,]
    }
  }else{
    indlik <- indder <- NULL
  }
  return(list(lval=lval,dlval=dlval,d2lval=d2lval,indlik=indlik,indder=indder))
}


sublik_deriv.nhm <- function(nevent,times,low,up,lowt,upt,from,to,subs,npar,ncov,nstate,covvalue,par,tcrit,censor,censor.states,death,death.states,fishscore,rtol,atol,transitions,intens,finalcovs,gensolve,misc_version=FALSE) {
  #Function to obtain transition probability matrices (and their derivatives) for a particular covariate value
  pmat <- array(0,c(nstate,nstate,nevent))
  if (is.null(censor)) censor<-0
  if (death | !identical(censor,0)) fishscore <- FALSE
  if (ncov>0) {
    parms <- c(npar,ncov,nstate,covvalue,par)
  }else{
    parms <- c(npar,ncov,nstate,par)
  }
  #Solve extended system of differential equations from time t_0
  if (length(times)>1) {
    res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*npar)),times=times,func=genintens_deriv.nhm,rtol=rtol,atol=atol,tcrit=tcrit,parms=parms,intens=intens)
  }else{
    #If the unique set contains only one time point it is not necessary to solve any equations since the initial conditions will hold.
    res <- array(0,c(1, 1+ nstate^2*(npar +1)))
    res[1,] <- c(times[1], c(c(diag(nstate)),rep(0,nstate^2*npar)))
  }
  p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
  lik <- 0
  dlik <- array(0,c(nevent,npar))
  if (fishscore) d2lik <- array(0,c(nevent,npar,npar))
  dp0 <- array(0,c(length(times),nstate,nstate,npar))
  for (l in 1:npar) {
    q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times),nstate,nstate))
    dp0[,,,l] <- q0
  }
  #Obtain transition probabilities and their derivatives.
  dpmat <- array(0,c(nstate,nstate,nevent,npar))
  for (i in 1:nevent) {
    if (abs(det(p0[low[i],,]))>2e-16) {
      if (low[i] > 1) {
        invp <- gensolve(p0[low[i],,])
      }else{
        invp <- diag(nstate)
      }
      pmat[,,i] <- invp%*%p0[up[i],,]
      for (l in 1:npar) {
        dpmat[,,i,l] <- invp%*%(dp0[up[i],,,l] - dp0[low[i],,,l]%*%pmat[,,i])
      }
    }else{
      error_mess <- paste("Matrix singular when inverting transition probability matrix from ",min(times)," to ",times[low[i]],". Add a split between these points using the splits option, or choose better starting values.",sep="")
      stop(error_mess)
      #pmat[,,i] <- 1e-100
      #dpmat[,,i,] <- 1e-100
    }
    if (!misc_version) {
      cen <- FALSE
      for (k in 1:length(censor)) {
        if (to[i]==censor[k]) {
          lik[i] <- -log(sum(pmat[from[i],censor.states[[k]],i]))
          for (l in 1:npar) {
            dlik[i,l] <- -sum(dpmat[from[i],censor.states[[k]],i,l])/exp(-lik[i])
          }
          cen <- TRUE
        }
      }
      if (to[i]%in%death.states & death) {
        ################
        #covvalue needs to refer to final covvalue not covvalue during the interval
        covvalueD <- c(unlist(finalcovs[subs[i],]))
        Q <- intens(upt[i],covvalueD,par)
        ################
        q <- Q$q
        qp <- Q$qp
        p <- pmat[from[i],1:(nstate-1),i]%*%q[1:(nstate-1),to[i]]
        lik[i] <- -log(p)
        for (l in 1:npar) {
          dp <- dpmat[from[i],-to[i],i,l]%*%q[-to[i],to[i]] + pmat[from[i],-to[i],i]%*%qp[-to[i],to[i],l]
          dlik[i,l] <-  -dp/p
        }
      }
      if ((!to[i]%in%death.states | !death) & !cen) {
        lik[i] <- -log(pmat[from[i],to[i],i])
        dlik[i,] <-  -dpmat[from[i],to[i],i,]/pmat[from[i],to[i],i]
        if (fishscore) {
          ks <- array(0,c(npar,npar))
          for (l in 1:nstate) {
            ks <- ks + dpmat[from[i],l,i,]%*%t(dpmat[from[i],l,i,])/(pmat[from[i],l,i] + 1*(pmat[from[i],l,i]==0))
          }
          d2lik[i,,] <- ks
        }
      }
    }
  }
  if (misc_version) {
    out <- list(pmat=pmat, dpmat = dpmat)
    return(out)
  }else{
    lik <- replace(lik,is.na(lik),1e10)
    dlik <- replace(dlik,is.na(dlik),1e10)
    if (fishscore) {
      d2lik <- replace(d2lik,is.na(d2lik),1e10)
      out <- list(lik=lik,dlik=dlik,d2lik=d2lik)
    } else{
      out <- list(lik=lik,dlik=dlik)
    }
    if (transitions) out <- list(pmat=pmat, dpmat = dpmat)
  }
  out
}


genintens.nhm <- function(t,y,parms,intens) {
  #Wrapper function for the user supplied intens when derivatives not required.
  npar <- parms[1]
  ncov <- parms[2]
  nstate <- parms[3]
  p <- array(y,dim=c(nstate,nstate))
  if (ncov>0) {
    z <- parms[4:(3+ncov)]
  }else{
    z <- NULL
  }
  x <- parms[(4+ncov):(3+ncov+npar)]
  Q <- intens(t,z=z,x=x)
  q <- Q$q
  vec <- c((p%*%q))
  list(vec,NULL)
}

genintens_deriv.nhm <- function(t,y,parms,intens) {
  #Wrapper function for the user supplied intens
  npar <- parms[1]
  ncov <- parms[2]
  nstate <- parms[3]
  p <- array(y[1:(nstate^2)],dim=c(nstate,nstate))
  pp <- array(0,c(nstate,nstate,npar))
  for (j in 1:npar) {
    pp[,,j] <- array(y[(1 + j*nstate^2):((j+1)*nstate^2)],dim=c(nstate,nstate))
  }
  if (ncov>0) {
    z <- parms[4:(3+ncov)]
  }else{
    z <- NULL
  }
  x <- parms[(4+ncov):(3+ncov+npar)]
  Q <- intens(t,z=z,x=x)
  q <- Q$q
  qp <- Q$qp
  vec <- c((p%*%q))
  for (j in 1:npar) {
    vec <- c(vec,c(p%*%qp[,,j] + pp[,,j]%*%q))
  }
  list(vec,NULL)
}

#Routine to find the product of pmats and the corresponding first derivatives
#Only used for left-truncation (especially if left-truncation is from a later time point)
prodmat <- function(A,dA) {
  #Assume have an R x R x J array
  #And that
  if (dim(A)[3] ==1) return(list(P=A[,,1],dP=dA[,,1,]))
  P <- A[,,1]
  dP <- dA[,,1,]
  for (j in 2:dim(A)[3]) {
    dPnew <- array(0,dim(dP))
    for (k in 1:dim(dA)[4]) {
      dPnew[,,k] <- P%*%dA[,,j,k] + dP[,,k]%*%A[,,j]
    }
    dP <- dPnew
    P <- P%*%A[,,j]
  }
  return(list(P=P,dP=dP))
}



forward_indiv <- function(k, subject, time, state, pmat, dpmat, cov2, x, intens, nostate, nstate, npar, init, dinit, ctrunc, dctrunc, emat_nhm) {
  #subset <- which(subject==k) #Better to have stored this before?
  subset <- (attr(subject,"substart")[k]):(attr(subject,"substop")[k])
  #Relies on all unique subjects being in the subjectl vector when it is created.
  obs <- state[subset]
  pmat0 <- pmat[,,subset]
  dpmat0 <- dpmat[,,subset,]
  es <- emat_nhm(state[subset],time[subset],cov2[subset,,drop=FALSE],x,intens)
  #e<-es$e
  #de<-es$de
  e0 <-es$e #e[,subset]
  de0<-es$de #de[,subset,]
  pass <- forward_alg_nhm(pmat0,dpmat0,init[,k],dinit[,,k],e0,de0,length(obs),nostate,nstate,npar)
  CLIK <-log(pass$prob) - ctrunc[k]
  CDLIK<-pass$dprob/pass$prob - dctrunc[k,]
  return(list(CLIK=CLIK, CDLIK=CDLIK))
}


genlikelihood_misc.nhm <- function(x,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit=NULL,fishscore,obslist,state,subject,time,cov2,init_state,indout=FALSE,sublist=NULL,nsub=NULL,rtol,atol,intens,gensolve,emat_nhm,ncores=1,firstobs,initp,initcovs=NULL,finalcovs=NULL,initp_nhm=NULL,nparQ) {

  if (is.null(censor)) censor <- 0
  par <- x
  if (ncov==0) covvalue <- NULL
  nobs <- length(state)
  pmat <-array(0,c(nstate,nstate,nobs))
  dpmatQ<-array(0,c(nstate,nstate,nobs,nparQ))

  if (ncores==1) {
    genout <- lapply(1:ncovvals,genlikelihoodwrapper_misc.nhm, par=par,npar=nparQ,ncov=ncov,nstate=nstate,ncovvals=ncovvals,covs=covs,covtable=covtable,timesets=timesets,fromtimeindex=fromtimeindex,totimeindex=totimeindex,fromtimesets=fromtimesets,totimesets=totimesets,fromsets=fromsets,tosets=tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tcrit,fishscore=fishscore,indout=indout,sublist=sublist,nsub=nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve)
  }else{
    genout <- mclapply(1:ncovvals,genlikelihoodwrapper_misc.nhm, mc.cores = ncores,par=par,npar=nparQ,ncov=ncov,nstate=nstate,ncovvals=ncovvals,covs=covs,covtable=covtable,timesets=timesets,fromtimeindex=fromtimeindex,totimeindex=totimeindex,fromtimesets=fromtimesets,totimesets=totimesets,fromsets=fromsets,tosets=tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tcrit,fishscore=fishscore,indout=indout,sublist=sublist,nsub=nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve)
  }
  for (j in 1:ncovvals) {
    pmat[,,obslist[[j]]]<-genout[[j]]$pmat
    dpmatQ[,,obslist[[j]],]<-genout[[j]]$dpmat
  }

  #Now expand out the dpmat to include the unused parameters
  dpmat <- array(0,c(nstate,nstate,nobs,npar))
  dpmat[,,,1:nparQ] <- dpmatQ

  if (firstobs=="exact") {
  init <- array(0,c(nstate,nsub))
  init[cbind(init_state,1:nsub)]<-1
  #init<-c(1,rep(0,nstate-1))
  dinit<-array(0,c(nstate,npar,nsub))
  #es<-emat_nhm(x)
  }else{
    #No difference here on whether the first obs is absent or used.
    #Depends on whether the initp is to be estimated and if so whether covariates apply to it.
      if (is.null(initp_nhm)) {
        init <- array(rep(initp,nsub),c(nstate,nsub))
        dinit<-array(0,c(nstate,npar,nsub))
      }else{
        initP <- initp_nhm(nsub,initcovs,par)
        init <- initP$init
        dinit <- initP$dinit
      }
  }


  #Note that emat_nhm also depends on the main covariates...
  #In general we need emat_nhm to return:
  #e: as an R x R x nsub matrix
  #de: as an R x R x npar x nsub matrix

  nostate <- attr(emat_nhm,"nostate")

  #STORING ALL EMATRICES IN ONE GO IS POTENTIALLY VERY MEMORY INTENSIVE.
  #es <- emat_nhm(state,time,cov2,x,intens)
  #e<-es$e
  #de<-es$de
  clik<-0
  cdlik<-array(0,c(nsub,npar))

  ctrunc <- rep(0,nsub)
  dctrunc <- array(0,c(nsub, npar))

  if (!is.null(attr(emat_nhm,"ltrunc"))) {
    ltrunc <- attr(emat_nhm,"ltrunc")
   #Routine to get the relevant left-truncation terms...
   for (k in 1:nsub) {
     #if (ltrunc$ninitial[k]) {      #NEED TO BE CAREFUL SINCE SOME SUBJECTS WILL HAVE FIRST OBSERVATION AT THE INITIAL TIME #DOes that matter?
     relob <- which(subject==k)[1:ltrunc$ltruncation_entry]
     Ps <- prodmat(pmat[,,relob,drop=FALSE], dpmat[,,relob,,drop=FALSE])
     pmatv <- Ps$P
     dpmatv <- Ps$dP
     #Then need to compute
     dptrunc <- rep(0,npar)
     ptrunc <- sum((init[,k]%*%pmatv)[ltrunc$ltruncation_states])
     #and also
     for (j in 1:npar) {
     #This needs to also depend on dinit...
     dptrunc[j] <- sum((init[,k]%*%dpmatv[,,j])[ltrunc$ltruncation_states]) + sum((dinit[,j,k]%*%pmatv)[ltrunc$ltruncation_states])
     ###############
     ###############


     }
     ctrunc[k] <- log(ptrunc)
     dctrunc[k,] <- dptrunc/ptrunc
     #}
   }
  }


  CLIK_lists <- mclapply(1:nsub, forward_indiv, mc.cores=ncores, subject=subject, time=time, state=state, pmat=pmat, dpmat=dpmat, cov2=cov2, x=x, intens=intens, nostate=nostate,nstate=nstate,npar=npar,init=init,dinit=dinit,ctrunc=ctrunc,dctrunc=dctrunc,emat_nhm=emat_nhm)
  for (k in 1:nsub) {
    clik[k] <- CLIK_lists[[k]]$CLIK
    cdlik[k,] <- CLIK_lists[[k]]$CDLIK
  }

  out <- sum(clik)
  attr(out,"gradient") <- apply(cdlik,2,sum)
  if (indout) {
    attr(out,"indlik") <- clik
    attr(out,"indder") <- cdlik
  }
  out
}

genlikelihoodwrapper_misc.nhm <- function(j, par,npar,ncov,nstate,ncovvals,covs,covtable,timesets,fromtimeindex,totimeindex,fromtimesets,totimesets,fromsets,tosets,death,death.states,censor,censor.states,tcrit=NULL,fishscore,indout,sublist,nsub=NULL,rtol,atol,intens,finalcovs,gensolve) {
  #Cycles through the set of unique covariate patterns
  if (ncov>0) covvalue <- unlist(covs[j,])
  nevent <- covtable[j]
  times <- timesets[[j]]
  low <- fromtimeindex[[j]]
  up <- totimeindex[[j]]
  lowt <- fromtimesets[[j]]
  upt <- totimesets[[j]]
  from <- fromsets[[j]]
  to <- tosets[[j]]
  subs <- sublist[[j]]
  out <- sublik_deriv.nhm(nevent,times,low,up,lowt,upt,from,to,subs,npar,ncov,nstate,covvalue,par,tcrit,censor,censor.states,death,death.states,fishscore=FALSE,rtol,atol,transitions=FALSE,intens,finalcovs,gensolve,misc_version=TRUE)
  #out <- sublik_deriv_misc.nhm(nevent,times,low,up,lowt,upt,from,to,npar,ncov,nstate,covvalue,par,tcrit,censor,censor.states,death,death.states,rtol,atol,intens,gensolve)
  return(out)
}

#This is the proper package version...

forward_alg_nhm <- function(pmat0,dpmat0,init,dinit,e,de,nob,nost,nlst,npar) {
  dprob<-rep(0,npar)
  prob<-0
  alp<-array(0,c(nob+1,nlst))
  dalp<-array(0,c(nob+1,nlst,npar))
  v1<-as.double(init)
  v2<-as.double(dinit)
  v4<-as.double(e)
  v5<-as.double(de)
  v6<-as.double(alp)
  v7<-as.double(dalp)
  v8<-as.double(pmat0)
  v9<-as.double(dpmat0)
  v10<-as.integer(nob)
  v11<-as.integer(nlst-1)
  v12<-as.integer(npar-1)
  v13<-as.double(prob)
  v14<-as.double(dprob)
  forward <- .C("nhm_forwardalg",v1,v2,v4,v5,v6,v7,v8,v9,v10,v11,v12,prob=v13,
             dprob=v14,PACKAGE="nhm")
  list(prob=forward$prob,dprob=forward$dprob)
}
#
#
# forward_alg_nhm <- function(pmat0,dpmat0,init,dinit,e,de,nob,nost,nlst,npar) {
#   dprob<-rep(0,npar)
#   prob<-0
#   alp<-array(0,c(nob+1,nlst))
#   dalp<-array(0,c(nob+1,nlst,npar))
#   v1<-as.double(init)
#   v2<-as.double(dinit)
#   v4<-as.double(e)
#   v5<-as.double(de)
#   v6<-as.double(alp)
#   v7<-as.double(dalp)
#   v8<-as.double(pmat0)
#   v9<-as.double(dpmat0)
#   v10<-as.integer(nob)
#   v11<-as.integer(nost-1)
#   v12<-as.integer(nlst-1)
#   v13<-as.integer(npar-1)
#   v14<-as.double(prob)
#   v15<-as.double(dprob)
#   forward <- .C("forwardalg_extend2",v1,v2,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,prob=v14,
#                 dprob=v15)
#   list(prob=forward$prob,dprob=forward$dprob)
# }
#
