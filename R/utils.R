#Function to plot prevalence estimates and asymptotic confidence intervals.
#To do: allow covvalue to be a data frame so that multiple values can be computed on the same plot.
plot.nhm <- function(x, #nhm fitted object
                     what="probabilities", #what to plot - either probabilities or intensities or stateoccup
                     time0=0, #Time from which to compute transition probabilities
                     state0=1, #Starting state (ignored if what=stateoccup)
                     times=NULL, #Times to compute the transition probability.
                     covvalue=NULL, #Covariate value. Should be a vector of length ncov
                     ci=TRUE, #Whether confidence intervals are required
                     sim=FALSE, #Whether method of simulating from multivariate Normal distribution should be used for CIs
                     coverage=0.95, #Nominal coverage for confidence limits.
                     B=1000, #Number of simulations if simulation method to be used.
                     rtol=1e-6, #Relative error tolerance to be passed to lsoda
                     atol=1e-6, #Absolute error tolerance to be passed to lsoda
                     initp=NULL, #optional vector of initial state occupation probabilities. If NULL then will use the estimates from the model. If original model was left-truncated will assume probabilities at tstart correspond to those implied by the left-truncation model. If ltrunc supplied will similarly calculate based on left-truncation from the value of t0 supplied. If model does not include misclassification, will assume entry in state 1.
                     ltrunc=NULL, #optional list containing "ltruncation_states": a vector specifying the indices of the non-absorbing states. "ltruncation_time": the time/age at which patients are left-truncation. Note: this list can be omitted if original model involved left-truncation.
                     main_arg=NULL, #Start to give to the plot titles.
                     xlab="Time",
                     statemerge=FALSE, #whether the states should be merged if phasemap is used in the model fitting.
                     ...
) {
  if (statemerge) {
    phasemap <- x$model_object$phasemap
    nostate <- x$model_object$nostate
    if (is.null(phasemap) | is.null(nostate)) {
      warning("No phasemap specified in the original fitted model, change statemerge to FALSE.")
      phasemap <- 1:(x$nstate)
      nostate <- x$nstate
    }
    phasetrans <- phasemaptrans(phasemap)
  }else{
    nostate <- x$nstate
    phasetrans <- diag(nostate)
  }
  
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  model <- x
  if (!is.null(model$fixedpar)) stop("Plot not currently available for models with fixedpar")
  if (ci & model$singular) {
    warning("Confidence intervals cannot be calculated due to a singular Hessian")
    ci <- FALSE
  }
  if (what=="probabilities") {
    if (is.null(main_arg)) main_arg <- "Probability in state "
    preds <- predict(model,time0,state0,times,covvalue,ci,sim,coverage,B,rtol,atol,statemerge)
    if (!statemerge) {
      nstate <- model$nstate
    }else{
      nstate <- nostate
    }
    prevalence <- preds$probabilities
    times <- preds$times
    prevL <- preds$lower
    prevU <- preds$upper
    par(mfrow=c(2,ceiling(nstate/2)))
    for (i in 1:nstate) {
      plot(times,prevalence[,i],type="l",xlab=xlab,ylab="Probability",main=paste(main_arg,i,sep=""),ylim=c(0,1))
      if (ci) {
        lines(times,prevL[,i],lty=2)
        lines(times,prevU[,i],lty=2)
      }else{
        prevL<-prevU<-NULL
      }
    }
  }else{
    if (what=="intensities") {
      if (is.null(main_arg)) main_arg <- "Intensity "
      preds <- qmatrix.nhm(model,time0,times,covvalue,ci,sim,coverage,B)
      ntrans <- dim(preds$intensities)[2]
      intensnames <- names(preds$intensities)
      intensities <- preds$intensities
      times <- preds$times
      prevL <- preds$lower
      prevU <- preds$upper
      if (ci) {
        maxintens <- apply(prevU,2,stats::quantile,0.98)
      }else{
        maxintens <- apply(intensities,2,stats::quantile,0.98)
      }
      par(mfrow=c(2,ceiling(ntrans/2)))
      for (i in 1:ntrans) {
        plot(times,intensities[,i],type="l",xlab=xlab,ylab="Intensity",main=paste(main_arg,intensnames[i],sep=" "),ylim=c(0,maxintens[i]))
        if (ci) {
          lines(times,prevL[,i],lty=2)
          lines(times,prevU[,i],lty=2)
        }else{
          prevL<-prevU<-NULL
        }
      }
    }else{
      if (what=="stateoccup") {
        if (is.null(main_arg)) main_arg <- "Probability in state "
        preds <- state_occupation_probability.nhm(model,covvalue,time0, times, initp, ltrunc, rtol,atol,ci,sim,mode="main",B,coverage, statemerge)
        if (statemerge) {
          nstate <- model$model_object$nostate
        }else{
         nstate <- model$nstate
        }
        prevalence <- preds$ests
        times <- preds$times
        prevL <- preds$ests_low
        prevU <- preds$ests_high
        par(mfrow=c(2,ceiling(nstate/2)))
        for (i in 1:nstate) {
          plot(times,prevalence[,i],type="l",xlab=xlab,ylab="Probability",main=paste(main_arg,i,sep=""),ylim=c(0,1))
          if (ci) {
            lines(times,prevL[,i],lty=2)
            lines(times,prevU[,i],lty=2)
          }else{
            prevL<-prevU<-NULL
          }
        }
      }
      if (what!="stateoccup") stop("Only transition probabilities, state occupation probabilities or intensities can be plotted.")
    }
  }
  
}

predict.nhm <- function(object,
                        time0=0, #Time from which to compute transition probabilities
                        state0=1, #Starting state
                        times=NULL, #Times to compute the transition probability for
                        covvalue=NULL, #Covariate value. Should be a vector of length ncov
                        ci=TRUE, #Whether confidence intervals are required?n
                        sim=FALSE, #Whether method of simulating from multivariate Normal distribution should be used for CIs
                        coverage=0.95, #Nominal coverage for confidence limits.
                        B=1000, #Number of simulations if simulation method to be used.
                        rtol=1e-6, #Relative error tolerance to be passed to lsoda
                        atol=1e-6, #Absolute error tolerance to be passed to lsoda
                        statemerge=FALSE, #whether the states should be merged if phasemap is used in the model fitting.
                        ...){
  model <- object
  
  if (state0 > model$nstate) stop("Invalid starting state")
  
  if (statemerge) {
    phasemap <- model$model_object$phasemap
    nostate <- model$model_object$nostate
    if (is.null(phasemap) | is.null(nostate)) stop("No phasemap specified in the original fitted model, change statemerge to FALSE.")
    phasetrans <- phasemaptrans(phasemap)
  }else{
    nostate <- model$nstate
    phasetrans <- diag(nostate)
  }  
  if (!is.null(model$fixedpar)) stop("predict not currently available for models with fixedpar")
  if (time0==0 & model$model_object$type=="weibull") time0<-1e-5 #Avoid singularity at 0.
  if (is.null(times)) {
    times <- seq(time0,model$maxtime,length.out=100)
  }
  npar <- model$model_object$nparQ
  par <- model$par[1:npar]
  ncov <- model$ncov
  fisher <- model$hess
  intens <- model$intens
  nstate <- model$nstate
  if (is.null(fisher)) fisher <- model$fisher
  if (is.null(fisher) & ci) stop("Fisher information required for confidence intervals")
  if (model$singular & ci) {
    warning("Cannot calculate confidence intervals due to singular Hessian")
    ci <- FALSE
  }
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  tcrit <- max(model$tcrit, max(times)+2)
  if (is.null(tcrit)) tcrit<-max(times)+2
  if (ncov>0) {
    parms <- c(npar,ncov,nstate,covvalue,par)
  }else{
    parms <- c(npar,ncov,nstate,par)
  }
  if (min(times)!=time0) times<-c(time0,times)
  res <- deSolve::lsoda(y=c(diag(nstate)),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
  res2 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
  prevalence <- res2[,state0,]%*%phasetrans
  
  if (ci) {
    if (sim) {
      sig <- solve(fisher)[1:npar,1:npar]
      #Ensure symmetric
      sig <- 0.5*(sig + t(sig))
      newpar<-rmvnorm(n=B,mean=par,sigma=sig)
      prevB<-array(0,c(length(times),nostate,B))
      for (i in 1:B) {
        if (ncov>0) {
          parms <- c(npar,ncov,nstate,covvalue,newpar[i,])
        }else{
          parms <- c(npar,ncov,nstate,newpar[i,])
        }
        res <- deSolve::lsoda(y=c(diag(nstate)),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
        res2 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
        prevB[,,i] <- res2[,state0,]%*%phasetrans
      }
      
      const <- 0.5*(1-coverage)
      prevL<-apply(prevB,c(1,2),function(x) sort(x)[B*const])
      prevU<-apply(prevB,c(1,2),function(x) sort(x)[B*(1-const)])
    }else{
      #Delta method version.
      #Need first derivatives also...
      covmat<-solve(fisher)[1:npar,1:npar]
      res <- deSolve::lsoda(y=c(c(diag(nstate)),rep(0,nstate^2*npar)),times=times,func=genintens_deriv.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
      p0 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
      lik <- 0
      dp0 <- array(0,c(length(times),nstate,nstate,npar))
      for (l in 1:npar) {
        q0 <- array(res[,(2 + l*nstate^2):((l+1)*nstate^2+1)],c(length(times),nstate,nstate))
        dp0[,,,l] <- q0
      }
      dprev<-dp0[,state0,,]
      if (!statemerge) {
      vars<-array(0,c(length(times),nstate))
      for (k in 1:length(times)) {
        for (l in 1:nstate) {
          vars[k,l] <- (t(dprev[k,l,])%*%covmat%*%(dprev[k,l,]))
        }
      }
      }else{
        vars<-array(0,c(length(times),nostate))
        for (k in 1:length(times)) {
          pp <- t(phasetrans)%*%dprev[k,,]
          for (l in 1:nostate) {
            vars[k,l] <- (t(pp[l,])%*%covmat%*%(pp[l,]))
          }
        }
      }
      
      #Use a logit transformation for now.
      vars <- vars/(prevalence * (1 - prevalence) + 1*(prevalence==0) + 1*(prevalence==1))^2
      lprevalence <- log(prevalence) - log(1 - prevalence)
      
      const <- qnorm(1- 0.5*(1-coverage))
      prevL <- 1/(1 + exp( -(lprevalence - const*vars^0.5)))
      prevU <- 1/(1 + exp( -(lprevalence + const*vars^0.5)))
    }
  }else{
    prevL<-prevU<-NULL
  }
  list(times=times,probabilities=prevalence,lower=prevL,upper=prevU)
}


print.nhm <- function(x,ci=TRUE,...) {
  model <- x
  par <- model$par
  fisher <- model$hess
  #Add something that gives parameter names
  par_names <- model$parnames

  if (is.null(fisher)) fisher <- model$fisher_info
  if (is.null(fisher) & ci) stop("Fisher information required for confidence intervals")
  if (is.null(par_names)) par_names <- sapply(1:length(par),function(x) paste("Parameter",x))
  if (model$singular & ci) {
    warning("Cannot calculate confidence intervals due to a singular Hessian")
    ci <- FALSE
  }
  if (ci) {
  std.err <- diag(solve(fisher))^0.5
    mat <- round(cbind(par, par - qnorm(0.975)*std.err,par + qnorm(0.975)*std.err),4)
    dimnames(mat)[[2]] <- c("Est","Low 95%","Up 95%")
  }else{
    mat <- round(cbind(par),4)
    dimnames(mat)[[2]] <- c("Est")
  }
  if (!is.null(model$fixedpar)) {
    mat2 <- array(NA,c(model$npar,dim(mat)[2]))
    mat2[-model$fixedpar,] <- mat
    mat2[model$fixedpar,1] <- model$fixval
    dimnames(mat2)[[2]] <- dimnames(mat)[[2]]
    mat <- mat2
  }
  dimnames(mat)[[1]] <- par_names
  print(mat)
  cat(paste("\nDeviance:",round(2*model$value,2)))
  if (!is.null(model$fixedpar)) {
  cat(paste("\nParameter(s):",paste(par_names[model$fixedpar],collapse=","),"treated as fixed.",sep=" "))
  }
}

AIC.nhm <- function(object,...,k=2) {
  if (!is.null(object$fixedpar)) stop("AIC not currently available for models with fixedpar")
  dotargs <- list(object,...)
  named <- if (is.null(names(dotargs)))
    rep_len(FALSE, length(dotargs)) else (names(dotargs) != "")
  if(any(named)) {
    warning("the following arguments to 'AIC.nhm' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse=", "))
  }
  dotargs <- dotargs[!named]
  is.nhm <- vapply(dotargs,function(x) inherits(x,"nhm"), NA)
  dotargs <- dotargs[is.nhm]
  AICs <- sapply(dotargs,function(x) 2*x$value + k*x$npar)
  AICs
}

anova.nhm <- function(object, ...) {
  if (!is.null(object$fixedpar)) stop("anova not currently available for models with fixedpar")
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep_len(FALSE, length(dotargs)) else (names(dotargs) != "")
  if(any(named)) {
    warning("the following arguments to 'anova.nhm' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse=", "))
  }
  dotargs <- dotargs[!named]
  is.nhm <- vapply(dotargs,function(x) inherits(x,"nhm"), NA)
  dotargs <- dotargs[is.nhm]
  if (length(dotargs)==0) stop("More than one nhm fitted object must be supplied.")
  if (length(dotargs)>=1) {
    deviances <- c(2*object$value, sapply(dotargs,function(x) 2*(x$value)))
    dev_difference <- sapply(dotargs,function(x) 2*(x$value - object$value))
    par_difference <- sapply(dotargs,function(x) (x$npar - object$npar))
    dev_difference <- -dev_difference * sign(par_difference)
    par_difference <- abs(par_difference)
    params <- c(object$npar,sapply(dotargs,function(x) x$npar ))
    output <- data.frame("Dev"=deviances,"npar"=params,"Lik Ratio"=c(NA,dev_difference),"Df"=c(NA,par_difference),"p val"=c(NA,1-pchisq(dev_difference,par_difference)))
    row.names(output)<-paste("Model",1:length(deviances),sep=" ")
    class(output) <- c("anova","data.frame")
    return(output)
  }
}

########################################
###Work in progress#####################
########################################
logLik.nhm <- function(object, by.subject=FALSE, ...) {
  #Function to obtain the log-likelihood of the model
  if (!inherits(object,"nhm")) {
    stop("Object should be a fitted nhm object")
  }
  if (by.subject) {
    #Need to call bhhh.nhm with the right model outputs.
    #First need to re-process the data.
    control <- object$control
    tmax  <- control$tmax
    if (is.null(tmax)) tmax <- max(object$model_object$time)+1
    if (is.null(control)) control <- nhm.control()
    gensolve <- object$gensolve
    if (is.null(gensolve)) gensolve <- solve
    hidden <- object$model_object$hidden
    if (hidden) {
      emat_nhm <- object$model_object$emat_nhm
      ltrunc <- attr(emat_nhm,"ltrunc")
    }else{
      ltrunc <- NULL
    }
    processed <- dataprocess.nhm(object$model_object$state,object$model_object$time,object$model_object$subject,object$model_object$covariates,object$model_object$ncov,object$splits,object$model_object$firstobs,object$safe,ltrunc) #Include firstobs in the data processing step.
    finalcovs <- processed$finalcovs
    if (hidden) {
      init_state <- processed$init_state
      initcovs <- processed$initcovs
      #Additional processing for models with hidden states
      hiddenout <- hidden_process(object$model_object$state, object$model_object$time, object$model_object$subject, object$model_object$covariates, processed, object$model_object$firstobs, object$ncov)
      cov2 <- hiddenout$cov2
      statel <- hiddenout$statel
      timel <- hiddenout$timel
      subjectl <- hiddenout$subjectl
      out <- bhhh.nhm(object$par,npar=object$npar,ncov=object$ncov,nstate=object$nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=object$model_object$death,death.states=object$model_object$death.states,censor=object$model_object$censor,censor.states=object$model_object$censor.states,tcrit=tmax,fishscore=FALSE,obslist=processed$obslist,state=statel,subject=subjectl,time=timel,init_state=init_state,indout=TRUE,cov2=cov2,sublist=processed$sublist,nsub=processed$nsub,rtol=control$rtol,atol=control$atol,intens=object$model_object$intens,gensolve=gensolve,emat_nhm=emat_nhm,ncores=control$ncores,firstobs=object$model_object$firstobs,initp=object$model_object$initp,initcovs=initcovs,finalcovs=finalcovs,initp_nhm=object$model_object$initp_nhm,nparQ=object$model_object$nparQ,fixedpar=object$fixedpar,fixval=object$fixval)
      val <- -as.numeric(out)
      grad <- -attr(out,"gradient")
    }else{
      out <- bhhh.nhm(object$par,npar=object$npar,ncov=object$ncov,nstate=object$nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=object$model_object$death,death.states=object$model_object$death.states,censor=object$model_object$censor,censor.states=object$model_object$censor.states,tcrit=tmax,fishscore=FALSE,indout=TRUE,sublist=processed$sublist,nsub=processed$nsub,rtol=control$rtol,atol=control$atol,intens=object$model_object$intens,gensolve=gensolve,finalcovs=finalcovs,ncores=control$ncores,fixedpar=object$fixedpar,fixval=object$fixval)
      val <- as.numeric(out)
      grad <- attr(out,"gradient")
    }
    return(list(val=val, grad=grad))
  }else{
    val <- -object$value
    attr(val,"df") <- object$npar
    class(val) <- "logLik"
    return(val)
  }
}

scoreresid.nhm <- function(x, plot=FALSE) {
  if (!inherits(x,"nhm")) {
    stop("Object should be a fitted nhm object")
  }
  hess <- x$hess
  if (is.null(hess)) {
    stop("Unable to compte score residuals since no Hessian given.")
  }else{
    covmat <- solve(hess)
  }
  req <- logLik.nhm(x, by.subject=TRUE)
  derivs <- req$grad
  sres <- colSums(t(derivs) * covmat %*%t(derivs))
  if (plot) {
    plot(sres,type="n")
  }
  sres
}


print.nhm_score <- function(x,which_comp=NULL, test_name=NULL, ...) {
  #x: Fitted nhm_score object
  #which_comp: Optional vector identifying which parameters to test as null. By default will extract the parameters marked as "Nonhom" from the model.
  #test_name: Optional character string to give title to test output. Primarily used by informative_test.nhm
  object <- x
  if (!is.null(object$fixedpar)) {
    parclass <- object$parclass[-object$fixedpar]
    parnames <- object$parnames[-object$fixedpar]
  }else{
    parclass <- object$parclass
    parnames <- object$parnames
  }
  if (is.null(which_comp)) {
    if (!is.null(parclass)) {
      which_comp <- which(parclass=="Nonhom")
      if (length(which_comp)==0) stop("No parameters relating to non-homogeneity. Supply which_comp")
    }else{
      stop("Cannot identify parameters relating to non-homogeneity. Supply which_comp")
    }
    cat("Score test for non-homogeneity components \n")
  }else{
    if (!is.null(test_name)) {
      cat(paste(test_name,"\n",sep=""))
    }else{
      cat("Score test of the supplied components \n")
    }
  }
  if (length(which_comp)==0) {
    cat("No parameters to test.")
  }else{
    Istar <- (solve(object$I))[which_comp,which_comp,drop=FALSE]
    Iother <- (solve(object$I))[-which_comp,-which_comp,drop=FALSE]
    s_other <- c(t(object$g[-which_comp])%*%Iother%*%object$g[-which_comp])
    if (s_other/length(object$g[-which_comp]) > 0.1) warning("Gradient indicates the null model may not have been maximized. Score test assumes all parameters except the ones governing non-homogeneity have already been maximized")

    individual_z <- rep(0,length(which_comp))
    for (k in 1:length(which_comp)) {
      include <- c((1:dim(object$I)[1])[-which_comp],which_comp[k])
      #Rearrange I so k is last component
      I_k <- object$I[include,include,drop=FALSE]
      #Extract the diagonal component of I^-1 corresponding to k
      Istar_k <- diag((solve(I_k)))[length(include)]
      individual_z[k] <- object$g[which_comp[k]]*Istar_k^0.5
    }

    #individual_z <- object$g[which_comp]* diag(Istar)^0.5 #Implemented in v0.1.0
    combined_z <- c(t(object$g[which_comp])%*%Istar%*%(object$g[which_comp]))

    mat <- round(cbind(c(individual_z,combined_z),c(2*(1-pnorm(abs(individual_z))),1-pchisq(combined_z,df=length(which_comp)))),4)
    dimnames(mat)[[1]] <- c(parnames[which_comp],"Overall")
    dimnames(mat)[[2]] <- c("Stat","p-val")
    print(mat)
  }
}


qmatrix.nhm <- function(object, time0=0, times=NULL, covvalue=NULL, ci=TRUE, sim=FALSE, coverage=0.95, B=1000) {
  #Function to compute the transition intensities
  model <- object
  if (!is.null(object$fixedpar)) stop("qmatrix.nhm not currently available for models with fixedpar")
  if (is.null(times)) {
    times <- seq(time0,model$maxtime,length.out=100)
  }
  npar <- model$model_object$nparQ
  par <- model$par[1:npar]
  ncov <- model$ncov
  fisher <- model$hess
  intens <- model$intens
  nstate <- model$nstate
  trans <- model$model_object$trans
  R <- nstate
  #Should do every non-zero transition because covariates effects may be different.
  n1 <- (array(paste(rep(1:R,R),rep(1:R,each=R),sep="->"),c(R,R)))
  n1 <- n1[c(trans)!=0]
  if (is.null(fisher)) fisher <- model$fisher
  if (is.null(fisher) & ci) stop("Fisher information required for confidence intervals")
  if (ci & model$singular) {
    warning("Cannot compute confidence intervals with a singular Hessian.")
    ci <- FALSE
  }
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (length(covvalue)!=model$ncov) stop("covvalue must be a vector of length equal to the number of covariates")
  nstate < model$nstate
  Q <- sapply(times,function(t) intens(t,covvalue,par))
  q1 <- t(sapply(Q[1,],function(x) x))
  intensities <- data.frame(q1[,c(trans)!=0])
  names(intensities) <- n1
  row.names(intensities) <- round(times,3)
    if (ci) {
    if (sim) {
      sig <- solve(fisher)[1:npar,1:npar]
      #Ensure symmetric
      sig <- 0.5*(sig + t(sig))
      newpar<-rmvnorm(n=B,mean=par,sigma=sig)
      intB<-array(0,c(length(times),length(n1),B))
      for (i in 1:B) {
        Q <- sapply(times,function(t) intens(t,covvalue,newpar[i,]))
        q1 <- t(sapply(Q[1,],function(x) x))
        intB[,,i] <- q1[,c(trans)!=0]
      }
      const <- 0.5*(1-coverage)
      intensL<-data.frame(apply(intB,c(1,2),function(x) sort(x)[B*const]))
      intensU<-data.frame(apply(intB,c(1,2),function(x) sort(x)[B*(1-const)]))
    }else{
      #Delta method version.
      #Need first derivatives also...
      covmat<-solve(fisher)[1:npar,1:npar]
      dq1 <- array(sapply(Q[2,],function(x) x),c(nstate^2,npar,length(times)))
      dq1 <- dq1[c(trans)!=0,,]
      vars<-array(0,c(length(times),length(n1)))
      for (k in 1:length(times)) {
        for (l in 1:length(n1)) {
          vars[k,l] <- t(dq1[l,,k])%*%covmat%*%(dq1[l,,k])
        }
      }
      #Transform to the log-scale.
      vars <- vars/(intensities^2 + 1*(intensities==0))
      const <- qnorm(1- 0.5*(1-coverage))
      intensL <- intensities*exp(-const*vars^0.5)
      intensU <- intensities*exp(const*vars^0.5)
    }
    names(intensL) <- names(intensU) <- n1
    row.names(intensL) <- row.names(intensL) <- round(times,3)
  }else{
    intensL<-intensU<-NULL
  }
  list(times=times,intensities=intensities,lower=intensL,upper=intensU)
}


ematrix.nhm <- function(object, covvalue=NULL) {
  #Function to compute the misclassification probabilities and their SEs for a given covariate value
  model <- object
  if (!is.null(object$fixedpar)) stop("ematrix.nhm not currently available for models with fixedpar")
  if (object$model_object$hidden==FALSE) stop("Model does not include misclassification")
  par <- model$par
  npar <- model$npar
  ncov <- model$ncov
  fisher <- model$hess
  intens <- model$intens
  nstate <- model$nstate
  trans <- model$model_object$trans
  nostate <- max(model$model_object$phasemap)
  if (is.null(fisher)) fisher <- model$fisher
  if (is.null(fisher)) stop("Fisher information required for confidence intervals")
  if (model$singular) {
    stop("Cannot compute standard errors with a singular Hessian.")
  }
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (!is.null(covvalue) & ncov==0) warning("Model has no covariates: supplied values ignored.")
  R <- nostate
  #Should do every non-zero transition because covariates effects may be different.
  emat_nhm <- model$model_object$emat_nhm
  if (ncov>0) {
  E <- emat_nhm(state=1:R, t=0, z=array(rep(covvalue,each=R),c(R,length(covvalue))),x=par,intens=NULL,override=TRUE)
  }else{
  E <- emat_nhm(state=1:R, t=0, z=NULL,x=par,intens=NULL,override=TRUE)
  }
  #Just do a basic Delta method SE for now...
  emat <- E$e
  SEemat <- array(0,dim(emat))
  for (i in 1:R) {
    for (j in 1:R) {
       SEemat[i,j] <- sqrt(t(E$de[i,j,])%*%solve(model$hess)%*%(E$de[i,j,]))
    }
  }
  list(emat=emat,SEemat=SEemat)
}

initialprob.nhm <- function(object, covvalue=NULL) {
  #Function to compute the initial probabilities and their SEs for a given covariate value
  if (!is.null(object$fixedpar)) stop("initialprob.nhm not currently available for models with fixedpar")
  if (object$model_object$hidden==FALSE) stop("Model does not include misclassification")
  if (is.null(object$model_object$initp_nhm)) {
    warning("No parameters estimated for initial probabilities")
    return(list(initp=object$model_object$initp))
  }
  model <- object
  par <- model$par
  npar <- model$npar
  ncov <- model$ncov
  fisher <- model$hess
  nstate <- model$nstate
  trans <- model$model_object$trans

  if (is.null(fisher)) fisher <- model$fisher
  if (is.null(fisher)) stop("Fisher information required for confidence intervals")
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  R <- nstate
  #Should do every non-zero transition because covariates effects may be different.
  initp_nhm <- model$model_object$initp_nhm
  #Is somewhat annoying that have to deal with the intens function when shouldn't have to
  I <- initp_nhm(nsub=1, z=covvalue,x=par)

  #Just do a basic Delta method SE for now...
  initp <- I$initp[,1]
  if (model$singular) {
    SEinitp <- NULL
    warning("Cannot compute standard errors due to singular Hessian.")
  }else{
    SEinitp <- rep(0,R)
    for (i in 1:R) {
      SEinitp[i] <- sqrt(t(I$dinitp[i,,1])%*%solve(model$hess)%*%(I$dinitp[i,,1]))
    }
  }
  list(initp=initp,SEinitp=SEinitp)
}


print.nhm_model <- function(x,...) {
  model <- x
  if (model$hidden) {
  cat(paste("\nnhm model object for a misclassification hidden Markov model.\n"))
  cat(paste("\nModel has",model$npar,"unknown parameters."))
  }else{
  cat(paste("\nnhm model object for a Markov model.\n"))
  }
  cat(paste("\nModel for intensities has",model$nparQ,"unknown parameters and is of",model$type,"type.\n \n"))
  #Also add a parameter table.
  print(data.frame(Name=model$parnames,Type=model$parclass))
}





check_intens <- function(intens, ctimes, parval, cov_list) {
  #intens: supplied intensities function
  #ctimes: Times at which to evaluate the derivatives
  #parval: x value at which to evaluate the derivatives
  #cov_list: covariate values at which to evaluate the derivatives
  discr <- array(0,c(length(cov_list),length(ctimes),length(parval)))
  for (g in 1:length(cov_list)) {
    for (k in 1:length(ctimes)) {
      check_env <- new.env()
      check_env$t <- ctimes[k]
      check_env$z <- cov_list[[g]]
      check_env$x <- parval
      dQ <- attr(stats::numericDeriv(quote(intens(t,z,x)$q),"x",check_env),"gradient")
      dQ2 <- intens(ctimes[k],cov_list[[g]],parval)$qp
      discr[g,k,] <- apply(dQ-dQ2, 3, function(x) max(abs(x)))
    }
  }
  if (max(discr) > 0.01) {
    I <- which.max(apply(discr,1,max))
    J <- which.max(apply(discr,2,max))
    dpars <- which(discr[I,J,] > 0.01)
    if (length(dpars)>1) {
      str <- "parameters"
    }else{
      str <- "parameter"
    }
    if (length(cov_list)>1) {
      warntext <- paste("Bespoke intens function likely to be incorrect: Substantial discrepancy between supplied and numerical derivatives at starting values when using z=(", paste(cov_list[[I]],collapse=","),") at time t=",round(ctimes[J],3)," for ",str,":",paste(dpars,collapse=","),".",sep="")
    }else{
      warntext <- paste("Bespoke intens function likely to be incorrect: Substantial discrepancy between supplied and numerical derivatives at time t=",round(ctimes[J],3)," for ",str,":",paste(dpars,collapse=","),".",sep="")
    }
    warning(warntext)
  }
  return(discr)
}

