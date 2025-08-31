find_initial <- function(model, covvalue=NULL, tstart, initp=NULL, ltrunc=NULL,rtol=1e-6,atol=1e-6,mode="main") {
  
  #model: fitted nhm object
  #covvalue: Specific parameter vector for the prediction (should be a vector of length equal to ncov). If omitted, covvalue the covariate means from the data will be used
  #tstart: Time at which to calculate the initial state occupation probabilities (conditional on left truncation)
  #initp: Optional vector of state occupation probabilities at the origin time. This only needs to be specified if it is different from what is fitted in the model.
  #ltrunc: Optional list containing ltruncation_time and ltruncation_states. If supplied will replace the values in the original model fit object. 
  #rtol, atol: relative and absolute tolerance for the solution of the differential equations.
  #mode: Argument for internal use to faciliate parametric bootstrapping: "main" ensures standard errors and calculated, if mode="boot" then standard errors are not calculated.
  
  truncation_states <- NULL
  if (!is.null(model$fixedpar)) stop("Function not currently available for models with fixedpar")
  if (tstart==0 & model$model_object$type=="weibull") tstart<-1e-5 #Avoid singularity at 0.
  npar <- model$model_object$npar
  nparQ <- model$model_object$nparQ
  parQ <- model$par[1:npar]
  par <- model$par
  ncov <- model$ncov
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  nstate <- model$nstate
  tcrit <- max(model$tcrit, model$maxtime+2)
  if (is.null(tcrit)) tcrit<-model$maxtime+2
  if (ncov>0) {
    parms <- c(nparQ,ncov,nstate,covvalue,parQ)
  }else{
    parms <- c(nparQ,ncov,nstate,parQ)
  }
  
  #i) Known state occupation probabilities at time tstart
  #ii) Estimated state occupation probabilities at time tstart
  #iii) If tstart!=NULL then need to calculate the initial state probability vector as well.
  
  
  #Cases:
  #initp supplied by the function call
  #with ltrunc supplied and t0 \neq tstart => define dinitp=0 and then use the ltrunc routine
  #without ltrunc supplied, or with t0=tstart. => define dinitp=0 and use the non-ltrunc routine
  #initp not supplied
  #emat model
  #with ltrunc => use ltrunc supplied in model.
  # t0 = tstart: extract initp from the model and use the non-ltrunc routine
  # t0 < tstart: extract initp from the model and use the ltrunc routine
  #with ltrunc, but ltrunc also supplied by function
  # Give warning then 
  # t0 = tstart: extract initp from the model and use the non-ltrunc routine
  # t0 < tstart: extract initp from the model and use the ltrunc routine
  #not an emat model
  ##Set initp=c(1,0,..,0) and dinitp=0 + give warning
  #with ltrunc supplied and t0 \neq tstart
  #Use ltrunc routine
  #without ltrunc, or with t0=tstart
  #Use non-ltrunc routine
  
  
  
  #Determine the method type:
  if (!is.null(initp)) {
    if (length(initp)!=nstate) stop("initp should be a vector of length equal to the number of states in the model.")
    if (!is.null(ltrunc)) {
      if (ltrunc$ltruncation_time < tstart) {
        routine <- "ltrunc"
        t0 <- ltrunc$ltruncation_time
        dinitp <- array(0,c(nstate,npar))
      }
      if (ltrunc$ltruncation_time >= tstart) {
        warning("Left truncation details ignored")
        routine <- "normal"
        dinitp <- array(0,c(nstate,npar))
      }
    }else{
      routine <- "normal"
      dinitp <- array(0,c(nstate,npar))
    }
  }else{
    if (!is.null(model$model_object$emat_nhm)) {
      if(!is.null(ltrunc)) {
        if (!is.null(attr(model$model_object$emat_nhm,"ltrunc"))) warning("Supplied ltrunc object used rather than that stored in the model")
        if (ltrunc$ltruncation_time < tstart) {
          routine <- "ltrunc"
          t0 <- ltrunc$ltruncation_time
        }else{
          routine <- "normal"
        }
      }else{
        ltrunc <- attr(model$model_object$emat_nhm,"ltrunc")
        if (!is.null(ltrunc)) {
          if (ltrunc$ltruncation_time < tstart) {
            routine <- "ltrunc"
            t0 <- ltrunc$ltruncation_time
          }else{
            routine <- "normal"
          }
        }else{
          routine <- "normal"
        }
      }
      #In all cases use initp_nhm to get the initial values...
      initp_nhm <- model$model_object$initp_nhm
      I <- initp_nhm(nsub = 1, z = covvalue, x = par)
      initp <- c(I$initp)
      dinitp <- I$dinitp[,,1]
      if (is.null(dinitp)) dinitp <- array(0,c(nstate,npar))
    }else{
      initp <- c(1,rep(0,nstate-1))
      dinitp <- array(0,c(nstate,npar))
      if (!is.null(ltrunc)) {
        if (ltrunc$ltruncation_time < tstart) {
          routine <- "ltrunc"
          t0 <- ltrunc$ltruncation_time
        }else{
          routine <- "normal"
          t0 <- tstart
        }
      }else{
        routine <- "normal"
        t0 <- 0
      }
      Wmessage <- paste("Patients assumed to be in state 1 at time", t0, "with probability 1.",sep=" ")
      warning(Wmessage)
    }
  }
  ####Left-truncation correction of the initial probs
  if (routine=="ltrunc") {
    truncation_states <- ltrunc$ltruncation_states
    if (is.null(truncation_states)) stop("ltrunc must include the ltruncation_states.")
    
    #Compute the left-truncation probabilities:
    times0 <- c(t0, tstart)
    if (mode=="main") {
      res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*nparQ)),times=times0,func=genintens_deriv.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
      p0 <- array(res[,2:(nstate^2+1)],c(length(times0),nstate,nstate))
      dp0 <- array(0,c(length(times0),nstate,nstate,npar))
      for (l in 1:nparQ) {
        q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times0),nstate,nstate))
        dp0[,,,l] <- q0
      }
    }else{
      res <- deSolve::lsoda(y=c(c(diag(nstate))),times=times0,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
      p0 <- array(res[,2:(nstate^2+1)],c(length(times0),nstate,nstate))
      dp0 <- array(0,c(length(times0),nstate,nstate,npar)) #Ignore the rest.
    }
    #Estimate the specific initial probabilities.
    itrunc <- rep(0,nstate)
    itrunc[truncation_states]<-1
    probA <- (t(initp)%*%p0[2,,])
    init_star <- rep(0,nstate)
    init_star[truncation_states] <- probA[truncation_states]/c((probA)%*%(itrunc))
    
    #Also need the corresponding derivatives...
    dinit_star <- array(0,c(nstate,npar))
    if (mode=="main") {
      for (i in 1:npar) {
        dprobA <- (dinitp[,i]%*%p0[2,,] + t(initp)%*%dp0[2,,,i])
        dinit_star[truncation_states,i] <- dprobA[truncation_states]/c(probA%*%itrunc) - probA[truncation_states] * c(dprobA%*%itrunc)/c(probA%*%itrunc)^2
      }
    }    
    initp <- init_star
    dinitp <- dinit_star
  }
  return(list(initp=initp, dinitp=dinitp,truncation_states=truncation_states))
}


state_life_expectancy <- function(model, covvalue=NULL, tstart=0, tmax=NULL, initp=NULL, npt=500, discount=NULL,utilities=NULL,ltrunc=NULL,rtol=1e-6,atol=1e-6,ci=TRUE, sim=FALSE, mode="main", B=1000,coverage=0.95) {
  
  #model: fitted nhm object
  #covvalue: Specific parameter vector for the prediction (should be a vector of length equal to ncov). If omitted, covvalue the covariate means from the data will be used
  #tstart: Time at which to calculate the initial state occupation probabilities (conditional on left truncation)
  #tmax: Maximum time to integrate when calculating the expectation (upper point for restricted life years)
  #initp: Optional vector of state occupation probabilities at the origin time. This only needs to be specified if it is different from what is fitted in the model.
  #npt: Total number of points to evaluate the occupation probabilities when estimating the integral
  #discount: Discounting function (if discounted LYs are required for a health economic evaluation). This function can also incorporate age-specific health utilities if required.
  #utilities: either a nstate length vector of state specific utilities in order to produce a single QALY, or a nstate x m matrix of utilities to give different quantities.
  #ltrunc: Optional list containing ltruncation_time and ltruncation_states. If supplied will replace the values in the original model fit object. 
  #rtol, atol: relative and absolute tolerance for the solution of the differential equations.
  #ci: Logical for whether confidence intervals should be calculated for the quantities.
  #sim: Whether CIs should be estimated using parametric bootstrap (if TRUE) or via the delta method (if FALSE).
  #mode: Argument for internal use to faciliate parametric bootstrapping: "main" ensures standard errors and calculated, if mode="boot" then standard errors are not calculated.
  #B: Number of simulations for the parametric bootstrap
  #coverage: Coverage probability for the confidence intervals.
  
  
  if (is.null(tmax)) tmax <- max(model$model_object$time)
  
  if (ci & !sim) {
    ci_type <- "asymptotic"
  }
  if (ci & sim) {
    ci_type <- "boot"
  }
  if (!ci) ci_type <- "none"
  
  
  if (mode=="main" & ci_type=="boot") {
    #Run B samples to get distribution 
    ests <- array(0,c(B, model$nstate))
    covmat <- solve(model$hess)
    if (!is.null(utilities)) {
      if (is.vector(utilities)) { 
        qaly_ests <- rep(0,B) #Assumes a single vector of utiltiies
      }else{
        qaly_ests <- array(0,c(B,dim(utilities)[2]))
      }
    }
    for (b in 1:B) {
      samp <- mvtnorm::rmvnorm(1,mean = model$par,sigma = covmat)
      modelstar <- model
      modelstar$par <- c(samp)
      A <- state_life_expectancy(modelstar, covvalue, tstart, tmax, initp, npt, discount,utilities,ltrunc,rtol,atol, ci=TRUE,sim=TRUE,mode="boot")
      ests[b,] <- A$est
      if (!is.null(utilities)) {
        if (is.vector(utilities))  {
          qaly_ests[b] <- A$qaly_est
        }else{
          qaly_ests[b,] <- A$qaly_est
        }
      }
      # if (!is.null(progress)) {
      #    if (b%%progress ==0) print(b)
      #  }
    }
  }
  if (mode=="boot" | ci_type!="asymptotic") {
    compute_se <- FALSE
  }else{
    compute_se <- TRUE
  }
  
  #NB: currently the utility weights are assumed known and no uncertainty accounted for in the variance.
  #However, should be possible to accommodate assuming utilities are estimated entirely independently of the model.
  init_object <- find_initial(model, covvalue, tstart, initp, ltrunc, rtol, atol, mode)
  initp <- init_object$initp
  dinitp <- init_object$dinitp
  
  npar <- model$model_object$npar
  nparQ <- model$model_object$nparQ
  parQ <- model$par[1:npar]
  par <- model$par
  ncov <- model$ncov
  fisher <- model$hess
  if (is.null(fisher)) fisher <- model$fisher
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  nstate <- model$nstate
  tcrit <- max(model$tcrit, max(times)+2)
  if (is.null(tcrit)) tcrit<-max(times)+2
  if (ncov>0) {
    parms <- c(nparQ,ncov,nstate,covvalue,parQ)
  }else{
    parms <- c(nparQ,ncov,nstate,parQ)
  }
  
  ###################
  times <-seq(tstart,tmax, length.out=npt+1)
  covmat<-solve(fisher)[1:npar,1:npar]
  if (compute_se) {
    res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*nparQ)),times=times,func=genintens_deriv.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
    for (l in 1:nparQ) {
      q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times),nstate,nstate))
      dp0[,,,l] <- q0
    }
  }else{
    res <- deSolve::lsoda(y=c(c(diag(nstate))),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
  }
  
  delta <- times[2]-times[1]
  d <- c(delta/2,rep(delta,length(times)-2),delta/2)
  Dmat0 <- array(rep(d,nstate),c(length(times),nstate))
  if (is.null(discount)) { 
    Dmat <- Dmat0
  }else{
    discounts <- discount(times)
    if (is.vector(discounts)) {
      Dmat <- Dmat0 * array(rep(discounts,nstate),c(length(times),nstate))
    }else{
      if (!is.matrix(discounts) | !identical(dim(discounts),dim(Dmat0))) stop("discounts function should either return a vector, or a matrix with nstate columns")
      Dmat <- Dmat0 * discounts
    }
  }
  
  #Dmat is a ntimes x nstate matrix of the weight to be given to each time point in the summation that approximates the integral
  #Dmat in general would be a combination of the trapezium rule weights plus the discounting weight plus the utility weights (i.e. only the latter will vary between states)
  #Compute the estimates...
  est <- rep(0,nstate)
  for (j in 1:length(times)) {
    est <- est + Dmat[j,] * c(initp%*%p0[j,,])
  }
  
  if (compute_se) {
    #####
    #Is it better to do the algebra for the derivatives first before applying the covmat transformation?
    #i.e. get d(\sum_j D_j * \pi %*% Pmat(t_0,t_j))/dthet_k 
    ders <- array(0,c(nstate,npar))
    for (m in 1:npar) { 
      for (j in 1:length(times)) {
        ders[,m] <- ders[,m] + Dmat[j,] * ((initp %*% dp0[j,,,m] + dinitp[,m]%*%p0[j,,]))
      }
    }
    
    #Then just apply each of those directly to covmat.
    
    est_cov <- ders%*%covmat%*%t(ders)
    
    est_low <- est - qnorm( 1 - (1-coverage)/2)*diag(est_cov)^0.5
    est_high <- est + qnorm( 1 - (1-coverage)/2)*diag(est_cov)^0.5
    
  }else{
    est_low <- est_high <- est_cov <- ders <- NULL
  }
  
  #If utilities would then do
  if (!is.null(utilities)) {
    qaly_est <- t(utilities)%*%est
    if (compute_se) {
      qaly_var <- t(utilities)%*%est_cov%*%utilities
      qaly_low <- qaly_est - qnorm( 1 - (1-coverage)/2)*diag(qaly_var)^0.5
      qaly_high<- qaly_est + qnorm( 1 - (1-coverage)/2)*diag(qaly_var)^0.5
    }else{
      qaly_var <- qaly_low <- qaly_high <-  NULL
    }
  }else{
    qaly_est <- qaly_var <- qaly_low <- qaly_high <- NULL
  }
  
  if (mode=="main" & ci_type=="boot") {
    est_cov <- var(ests)
    if (!is.null(utilities)) {
      qaly_var <- var(qaly_ests)
      qaly_low <- qaly_high <- rep(0,length(qaly_est))
      if (length(qaly_est) >1) {
        for (i in 1:length(qaly_est)) {
          q <- quantile(qaly_ests[,i],c((1-coverage)/2,1 - (1-coverage)/2))
          qaly_low[i] <- q[1]
          qaly_high[i] <- q[2]
        }
      }else{
        q <- quantile(qaly_ests,c((1-coverage)/2,1 - (1-coverage)/2))
        qaly_low <- q[1]
        qaly_high <- q[2]
      }
    }
    est_low <- est_high <- rep(0,length(est))
    for (i in 1:length(est)) {
      q <- quantile(ests[,i],c((1-coverage)/2,1 - (1-coverage)/2))
      est_low[i] <- q[1]
      est_high[i] <- q[2]
    }
  }
  if (mode=="main" & ci_type=="none") {
    qaly_low <- qaly_high <- est_low <- est_high <- NULL 
  }
  
  return(list(est=est, est_cov = est_cov, est_low=est_low, est_high=est_high, qaly_est=qaly_est, qaly_var=qaly_var, qaly_low=qaly_low, qaly_high=qaly_high, initp=initp, ders=ders))
}


#Function to calculate the expected hitting time to a given state assuming an initial starting distribution.
expected_hitting_time <- function(model, state, covvalue=NULL, tstart=0, tmax=NULL, initp=NULL, npt=500, ltrunc=NULL,rtol=1e-6,atol=1e-6, ci=TRUE, sim=FALSE, mode="main", B=1000,coverage=0.95) {
  
  #model: fitted model
  #state: the hitting state of interest. NB: This function only gives valid results for progressive models.
  #covvalue: covariate vector
  #tstart, tmax: Here tmax defines the conditional distribution 
  #initp: optional set of initial probability values at time t0.
  #npt: number of evaluation points used in the trapezium rule
  #ltrunc: optional left-truncation object. NB: standard left-truncation won't work since if tstart>t0 some patients will already have reached hitting state.
  #rtol, atol: For the solver
  #ci: Logical for whether confidence intervals should be calculated for the quantities.
  #sim: Whether CIs should be estimated using parametric bootstrap (if TRUE) or via the delta method (if FALSE).
  #mode: Argument for internal use to faciliate parametric bootstrapping: "main" ensures standard errors and calculated, if mode="boot" then standard errors are not calculated.
  #B: Number of simulations for the parametric bootstrap
  #coverage: Coverage probability for the confidence intervals.
  
  if (is.null(tmax)) tmax <- max(model$model_object$time)
  if (ci & !sim) {
    ci_type <- "asymptotic"
  }
  if (ci & sim) {
    ci_type <- "boot"
  }
  if (!ci) ci_type <- "none"
  if (mode=="main" & ci_type=="boot") {
    #Run B samples to get distribution 
    ests0 <- ests <- rep(0,B)
    covmat <- solve(model$hess)
    for (b in 1:B) {
      samp <- mvtnorm::rmvnorm(1,mean = model$par,sigma = covmat)
      modelstar <- model
      modelstar$par <- c(samp)
      #A <- expected_hitting_time(modelstar, state, covvalue, tstart, tmax, initp, npt, ltrunc,rtol,atol, ci_type="boot",mode="boot")
      A <- expected_hitting_time(modelstar, state, covvalue, tstart, tmax, initp, npt, ltrunc,rtol,atol, ci=TRUE, sim=TRUE, mode="boot")
      ests[b] <- A$est
      ests0[b] <- A$est0
    }
  }
  if (mode=="boot" | ci_type!="asymptotic") {
    compute_se <- FALSE
  }else{
    compute_se <- TRUE
  }
  
  
  if (length(state)>1) stop("state must be an integer scalar from 1 to nstate")
  
  #In general, need to know for which states, j, q_{jr}>0.  In theory can infer this from "trans" in the model object. 
  
  #Potentially, the whole first part of the function is the same as for the life expectancy function.
  #the only difference is in warning whether initp ends up with no zero entries for "state" itself
  
  init_object <- find_initial(model, covvalue, tstart, initp, ltrunc, rtol, atol, mode)
  initp <- init_object$initp
  dinitp <- init_object$dinitp
  
  if (initp[state] >0) warning("Initial probabilities give non-zero chance of onset already occurring. Estimates may not be meaningful.")
  
  
  npar <- model$model_object$npar
  nparQ <- model$model_object$nparQ
  parQ <- model$par[1:npar]
  par <- model$par
  ncov <- model$ncov
  fisher <- model$hess
  if (is.null(fisher)) fisher <- model$fisher
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  nstate <- model$nstate
  tcrit <- max(model$tcrit, tmax+2)
  if (ncov>0) {
    parms <- c(nparQ,ncov,nstate,covvalue,parQ)
  }else{
    parms <- c(nparQ,ncov,nstate,parQ)
  }
  
  ###################
  times <-seq(tstart,tmax, length.out=npt+1)
  covmat<-solve(fisher)[1:npar,1:npar]
  
  if (compute_se) {
    res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*nparQ)),times=times,func=genintens_deriv.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
    for (l in 1:nparQ) {
      q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times),nstate,nstate))
      dp0[,,,l] <- q0
    }
  }else{
    res <- deSolve::lsoda(y=c(c(diag(nstate))),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
  }
  
  ###################
  #Need to extract the intensity function and corresponding derivative values at all the timepoints
  ###################
  
  trans <- model$model_object$trans
  fromstates <- which(trans[,state]!=0)
  
  #Extract the relevant intens and dintens values at the times.
  intens_val <- lapply(times, model$model_object$intens, z=covvalue,x=parQ)
  intenslist <- dintenslist <- list()
  for (k in 1:length(fromstates)) {
    intenslist[[k]] <- rep(0,length(times))
    dintenslist[[k]] <- array(0,c(length(times),npar))
    for (i in 1:length(times)) {
      intenslist[[k]][i] <- intens_val[[i]]$q[fromstates[k],state]
      dintenslist[[k]][i,1:nparQ] <- intens_val[[i]]$qp[fromstates[k],state,]
    }
  }
  
  delta <- times[2]-times[1]
  d <- c(delta/2,rep(delta,length(times)-2),delta/2)
  D <- d*times
  D0 <- d
  
  
  #Dmat is a ntimes x nstate matrix of the weight to be given to each time point in the summation that approximates the integral
  #Dmat in general would be a combination of the trapezium rule weights plus the discounting weight plus the utility weights (i.e. only the latter will vary between states)
  
  #Compute the estimates...
  est0 <- est1 <- 0
  for (j in 1:length(times)) {
    h <- c(initp%*%p0[j,,])
    for (k in 1:length(fromstates)) {
      g <- h[fromstates[k]] * intenslist[[k]][j] 
      est1 <- est1 + D[j] * g
      est0 <- est0 + D0[j] * g
    }
  }
  
  if (compute_se) {
    ders0 <- ders1 <- rep(0, npar)
    for (j in 1:length(times)) {
      h <- c(initp%*%p0[j,,])
      for (m in 1:npar) { 
        hd <- c(dinitp[,m]%*%p0[j,,])
        dh <- c(initp%*%dp0[j,,,m])
        for (k in 1:length(fromstates)) {
          ders1[m] <- ders1[m] + D[j] *(h[fromstates[k]]*dintenslist[[k]][j,m] + hd[fromstates[k]]*intenslist[[k]][j] + dh[fromstates[k]]*intenslist[[k]][j])
          ders0[m] <- ders0[m] + D0[j] *(h[fromstates[k]]*dintenslist[[k]][j,m] + hd[fromstates[k]]*intenslist[[k]][j] + dh[fromstates[k]]*intenslist[[k]][j])
        }
      }
    }
    
    #Then just apply each of those directly to covmat.
    ders <- rbind(ders1,ders0)
    
    est_cov <- ders%*%covmat%*%t(ders)
    
    #log-estimate
    var_lest <- (c(1/est1,-1/est0))%*%est_cov%*%(c(1/est1,-1/est0))
    
    var_est0 <- est_cov[2,2]
    
    
    #direct estimate
    var_est <- c(1/est0, -est1/est0^2)%*%est_cov%*%c(1/est0, -est1/est0^2)
    
    est_low <- est1/est0 - qnorm( 1 - (1-coverage)/2)*(var_est)^0.5
    est_high <- est1/est0 + qnorm( 1 - (1-coverage)/2)*(var_est)^0.5
    
    
    est_low2 <- est1/est0*exp( - qnorm( 1 - (1-coverage)/2)*(var_lest)^0.5)
    est_high2 <- est1/est0*exp( qnorm( 1 - (1-coverage)/2)*(var_lest)^0.5)
    
  }else{
    if (mode!="main" | ci_type=="none") {
      var_est <- var_lest <- var_est0 <- est_low <- est_high <- est_low2 <- est_high2 <- NULL
    }else{
      var_est <- var(ests)
      var_lest <- var(log(ests))
      var_est0 <- var(ests0)
      est_low <- quantile(ests,(1-coverage)/2)
      est_high <- quantile(ests,1-(1-coverage)/2)
      est_low2 <- est_high2 <- NULL
    }
  }
  return(list(est=est1/est0, var_est = var_est, var_lest=var_lest, est_low=est_low, est_high=est_high, est_low2=est_low2, est_high2=est_high2, est0=est0, var_est0=var_est0, initp=initp))
}


state_occupation_probability.nhm <- function(model, covvalue=NULL, time0 = 0, times=NULL, initp=NULL, ltrunc=NULL,rtol=1e-6,atol=1e-6,ci=TRUE, sim=FALSE,mode="main", B=1000,coverage=0.95, statemerge=FALSE) {
  
  #Structure should broadly follow state_life_expectancy except don't need to do the integration.
  
  if (statemerge) {
    phasemap <- model$model_object$phasemap
    nostate <- model$model_object$nostate
    if (is.null(phasemap) | is.null(nostate)) {
      warning("No phasemap specified in the original fitted model, change statemerge to FALSE.")
      phasemap <- 1:(model$nstate)
      nostate <- model$nstate
    }
    phasetrans <- phasemaptrans(phasemap)
  }else{
    nostate <- model$nstate
    phasetrans <- diag(nostate)
  }
  
  if (ci & !sim) {
    ci_type <- "asymptotic"
  }
  if (ci & sim) {
    ci_type <- "boot"
  }
  if (!ci) {
    ci_type <- "none"
  }
  
  #Determine times using the same method as predict.nhm
  if (is.null(times)) {
    times <- seq(time0, model$maxtime, length.out=201)
  }else{
    if (times[1] < time0) stop("Prediction times must be greater than or equal to time0")
    if (length(times)>1) {
      if (min(diff(times)) < 0) times <- sort(times)
    }
  }
  
  
  if (mode=="main" & ci_type=="boot") {
    #Run B samples to get distribution 
    
    Bests <- array(0,c(B, length(times), nostate))
    covmat <- solve(model$hess)
    for (b in 1:B) {
      samp <- mvtnorm::rmvnorm(1,mean = model$par,sigma = covmat)
      modelstar <- model
      modelstar$par <- c(samp)
      A <- state_occupation_probability.nhm(modelstar, covvalue, time0, times, initp, ltrunc,rtol,atol, ci=TRUE,sim=TRUE,mode="boot",statemerge=FALSE)
      if (statemerge) {
        As <- A$ests%*%phasetrans
      }else{
        As <- A$ests
      }
      Bests[b,,] <- As
      # if (!is.null(progress)) {
      #    if (b%%progress ==0) print(b)
      #  }
    }
  }
  if (mode=="boot" | ci_type!="asymptotic") {
    compute_se <- FALSE
  }else{
    compute_se <- TRUE
  }
  
  init_object <- find_initial(model, covvalue, time0, initp, ltrunc, rtol, atol, mode)
  initp <- init_object$initp
  dinitp <- init_object$dinitp
  
  npar <- model$model_object$npar
  nparQ <- model$model_object$nparQ
  parQ <- model$par[1:npar]
  par <- model$par
  ncov <- model$ncov
  fisher <- model$hess
  if (is.null(fisher)) fisher <- model$fisher
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  nstate <- model$nstate
  tcrit <- max(model$tcrit, max(times)+2)
  if (is.null(tcrit)) tcrit<-max(times)+2
  if (ncov>0) {
    parms <- c(nparQ,ncov,nstate,covvalue,parQ)
  }else{
    parms <- c(nparQ,ncov,nstate,parQ)
  }
  
  ###################
  #times <-seq(tstart,tmax, length.out=npt+1)
  if (compute_se) {
    covmat<-solve(fisher)[1:npar,1:npar]
    res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*nparQ)),times=times,func=genintens_deriv.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
    for (l in 1:nparQ) {
      q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times),nstate,nstate))
      dp0[,,,l] <- q0
    }
  }else{
    res <- deSolve::lsoda(y=c(c(diag(nstate))),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
    p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
    dp0 <- array(0,c(length(times),nstate,nstate,npar))
  }
  
  ests <- array(0,c(length(times),nostate))
  for (j in 1:length(times)) ests[j,] <- (initp%*%p0[j,,])%*%(phasetrans)
  
  if (compute_se) {
    #####
    ests_low <- ests_high <- array(0,c(length(times),nostate)) ###
    ests_cov <- array(0,c(length(times),nostate,nostate)) ###
    for (j in 1:length(times)) {
      ders <- array(0,c(nstate,npar))
      for (m in 1:npar) { 
        ders[,m] <- initp %*% dp0[j,,,m] + dinitp[,m]%*%p0[j,,]
      }
      #Then just apply each of those directly to covmat.
      DD <- ders%*%covmat%*%t(ders)
      if (statemerge) DD <- t(phasetrans)%*%DD%*%(phasetrans)
      ests_cov[j,,] <- DD
      ests_low[j,] <- ests[j,] - qnorm( 1 - (1-coverage)/2)*diag(ests_cov[j,,])^0.5
      ests_high[j,] <- ests[j,] + qnorm( 1 - (1-coverage)/2)*diag(ests_cov[j,,])^0.5
    }
  }else{
    ests_low <- ests_high <- ests_cov <-NULL
  }
  
  
  if (mode=="main" & ci_type=="boot") {
    ests_low <- ests_high <- array(0,c(length(times),nostate))
    ests_cov <- array(t(apply(Bests,2,cov)),c(length(times),nostate,nostate))
    for (i in 1:nstate) {
      for (j in 1:length(times)) {
        q <- quantile(Bests[,j,i],c((1-coverage)/2,1 - (1-coverage)/2))
        ests_low[j,i] <- q[1]
        ests_high[j,i] <- q[2]
      }
    }
  }
  if (mode=="main" & ci_type=="none") {
    ests_low <- ests_high <- NULL 
  }
  out <- list(times= times, ests=ests, ests_cov = ests_cov, ests_low=ests_low, ests_high=ests_high, initp=initp)
  return(out)
}


#Auxiliary function to convert the phasemap into a transformation matrix
phasemaptrans <- function(phasemap) {
  nstate <- length(phasemap)
  nostate <- max(phasemap)
  trans <- array(0,c(nstate,nostate))
  trans[cbind(1:nstate, phasemap)] <- 1
  return(trans)
}

