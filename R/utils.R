#Function to plot prevalence estimates and asymptotic confidence intervals.
#To do: allow covvalue to be a data frame so that multiple values can be computed on the same plot.
plot.nhm <- function(x, #nhm fitted object
                     what="probabilities", #what to plot - either probabilities or intensities
                     time0=0, #Time from which to compute transition probabilities
                     state0=1, #Starting state
                     times=NULL, #Times to compute the transition probability.
                     covvalue=NULL, #Covariate value. Should be a vector of length ncov
                     ci=TRUE, #Whether confidence intervals are required
                     sim=FALSE, #Whether method of simulating from multivariate Normal distribution should be used for CIs
                     coverage=0.95, #Nominal coverage for confidence limits.
                     B=1000, #Number of simulations if simulation method to be used.
                     rtol=1e-6, #Relative error tolerance to be passed to lsoda
                     atol=1e-6, #Absolute error tolerance to be passed to lsoda
                     main_arg=NULL, #
                     xlab="Time",
                     ...
) {
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
  preds <- predict.nhm(model,time0,state0,times,covvalue,ci,sim,coverage,B,rtol,atol)
  nstate <- model$nstate
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
    }
  if (what!="intensities") stop("Only probabilities or intensities can be plotted.")
  }

}

predict.nhm <- function(object,
                        time0=0, #Time from which to compute transition probabilities
                        state0=1, #Starting state
                        times=NULL, #Times to compute the transition probability for
                        covvalue=NULL, #Covariate value. Should be a vector of length ncov
                        ci=TRUE, #Whether confidence intervals are required
                        sim=FALSE, #Whether method of simulating from multivariate Normal distribution should be used for CIs
                        coverage=0.95, #Nominal coverage for confidence limits.
                        B=1000, #Number of simulations if simulation method to be used.
                        rtol=1e-6, #Relative error tolerance to be passed to lsoda
                        atol=1e-6, #Absolute error tolerance to be passed to lsoda
                        ...){
  model <- object
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
  nstate < model$nstate
  tcrit <- model$tcrit
  if (is.null(tcrit)) tcrit<-max(times)+2
  if (ncov>0) {
    parms <- c(npar,ncov,nstate,covvalue,par)
  }else{
    parms <- c(npar,ncov,nstate,par)
  }
  if (min(times)!=time0) times<-c(time0,times)
  res <- deSolve::lsoda(y=c(diag(nstate)),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
  res2 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
  prevalence <- res2[,state0,]
  if (ci) {
    if (sim) {
      sig <- solve(fisher)[1:npar,1:npar]
      #Ensure symmetric
      sig <- 0.5*(sig + t(sig))
      newpar<-rmvnorm(n=B,mean=par,sigma=sig)
      prevB<-array(0,c(length(times),nstate,B))
      for (i in 1:B) {
        if (ncov>0) {
          parms <- c(npar,ncov,nstate,covvalue,newpar[i,])
        }else{
          parms <- c(npar,ncov,nstate,newpar[i,])
        }
        res <- deSolve::lsoda(y=c(diag(nstate)),times=times,func=genintens.nhm,tcrit=tcrit,rtol=rtol,atol=atol,parms=parms,intens=intens)
        res2 <- array(res[,2:(nstate^2+1)],c(length(times),nstate,nstate))
        prevB[,,i] <- res2[,state0,]
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
      vars<-array(0,c(length(times),nstate))
      for (k in 1:length(times)) {
        for (l in 1:nstate) {
          vars[k,l] <- t(dprev[k,l,])%*%covmat%*%(dprev[k,l,])
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

  if (is.null(fisher)) fisher <- model$fisher
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

print.nhm_score <- function(x,which_comp=NULL,...) {
  #x: Fitted nhm_score object
  #which_comp: Optional vector identifying which parameters to test as null. By default will extract the parameters marked as "Nonhom" from the model.
  object <- x
  if (!is.null(object$fixedpar)) {
    parclass <- object$parclass[-object$fixedpar]
    parnames <- object$parnames[-object$fixedpar]
  }else{
    parclass <- object$parclass
    parnames <- object$parnames
  }
  cat("Score test for non-homogeneity components \n")
  if (is.null(which_comp)) {
    if (!is.null(parclass)) {
      which_comp <- which(parclass=="Nonhom")
    }else{
      stop("Cannot identify parameters relating to non-homogeneity. Supply which_comp")
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

  if (is.null(fisher)) fisher <- model$fisher
  if (is.null(fisher)) stop("Fisher information required for confidence intervals")
  if (model$singular) {
    stop("Cannot compute standard errors with a singular Hessian.")
  }
  if (is.null(covvalue) & ncov>0) {
    covvalue <- model$covmeans
  }
  if (!is.null(covvalue) & ncov==0) warning("Model has no covariates: supplied values ignored.")
  R <- nstate
  #Should do every non-zero transition because covariates effects may be different.
  emat_nhm <- model$model_object$emat_nhm
  #Is somewhat annoying that have to deal with the intens function when shouldn't have to
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
  if (model$singular) {
    stop("Cannot compute standard errors with a singular Hessian.")
  }
  R <- nstate
  #Should do every non-zero transition because covariates effects may be different.
  initp_nhm <- model$model_object$initp_nhm
  #Is somewhat annoying that have to deal with the intens function when shouldn't have to
  I <- initp_nhm(nsub=1, z=covvalue,x=par)

  #Just do a basic Delta method SE for now...
  initp <- I$initp[,1]
  SEinitp <- rep(0,R)
  for (i in 1:R) {
      SEinitp[i] <- sqrt(t(I$dinitp[i,,1])%*%solve(model$hess)%*%(I$dinitp[i,,1]))
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
