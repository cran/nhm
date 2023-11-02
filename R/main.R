
model.nhm <- function(formula,  data, subject, covariates=NULL,  type, trans,nonh=NULL, covm=NULL, centre_time=NULL, emat=NULL, ecovm=NULL, firstobs=NULL, initp=NULL, initp_value=NULL, initcovm=NULL, splinelist=NULL,degrees=NULL,censor=NULL,censor.states=NULL,death=FALSE,death.states=NULL,intens=NULL) {
  #formula: A formula identifying the state and time variable names within the data e.g. state ~ time
  #data: data frame containing including the observed states, observation times, subject identifers and any covariates. NB: covariates cannot be factors.
  #subject: The name of the subject (id) variable within the data
  #covariates: A vector of names of covariates in the data
  #type: type of intensity model. "bespoke": user supplied, "weibull"/"gompertz": weibull/gompertz models. "bspline": b-spline basis model.
  #trans: square matrix of viable transitions. Impossible transitions should be 0. Others should be labelled consecutively. Labelling transitions with the same value assumes the parameter is shared.
  #nonh: square matrix to indicate non-homogeneous transitions. Impossible transitions or homogeneous transitions should be 0. Otherwise label consecutively. Labelling the same value implies the same non-homogeneity.
  #covm: Either a nstate x nstate x ncov array to indicate covariate effects, where ncov is the same as the length of the covariates vector OR a named list of nstate x nstate matrices to indicate covariate effects for individual covariates.
  #centre_time: For gomertz models only: value by which to centre times in the gompertz model (to improve convergence)
  #emat: square matrix of viable misclassification errors. Impossible errors should be 0. Others should be labelled consecutively. Labelling transitions with the same value assumes the parameter is shared.
  #ecovm: Either a nstate x nstate x ncov array to indicate covariate effects on misclassification, where ncov is the same as the length of the covariates vector OR a named list of nstate x nstate matrices to indicate misclassification covariate effects for individual covariates
  #firstobs: For misclassification models; whether the first observation is exactly observed ("exact") or is absent and as such the state should be ignored. Defaults to "exact".
  #initp: For misclassification models where firstobs="absent" or "misc", the initial state probabilites. Only required if initial probabilites to be estimated: A vector of length nstate indicating which parameters to be estimated. Same convention as with qmat. First entry should be zero.
  #initp_value: For misclassification models where firstobs="absent" or "misc": Value of the initial state probabilities, treated as fixed if initp is missing.
  #initcovm: For misclassification models where firstobs="absent" or "misc" and est_initp=TRUE: Either a nstate x ncov matrix to indicate covariate effects, or a named list of vectors of length nstate to indicate covariate effects for the initial probability vector. Will use the covariates supplied at the initial time for this.
  #splinelist: For bspline models only: list (of length equal to the number of nonhomogeneous transitions) of knot point locations including the boundary knots.
  #degrees: For bspline models only: optional vector (of length equal to number of nonhomogeneous transitions) of degrees of splines. Defaults to 3 if not specified.
  #censor: Vector of censor state indicators in the data
  #censor.states: List of vectors of states in which subject occupy if censored by corresponding censor state indicator. Can be a vector if only one censor state marker is present.
  #death: Whether exact death times are present in the data set
  #death.states: Which states have exact death times. Otherwise ignored.
  #intens: Optional supplied intensity function. See elsewhere for requirements on the function.

  #Run checks:
  if (!is.data.frame(data)) stop("data must be a data.frame object")

  call <- match.call()
  indx <- match(c("subject"), names(call), nomatch = 0)
  if (indx!=0) {
    temp <- call[indx]
    subjectvar <- as.character(temp$subject)
  }else{
    stop("subject variable not specified")
  }
  indexes <- match(subjectvar,names(data), nomatch=0)
  if (indexes==0) stop("Subject variable not found")
  data$subject <- data[,indexes]
  #Check trans matrix
  if (!is.matrix(trans)) stop("trans must be a matrix")
  if (dim(trans)[1]!=dim(trans)[2]) stop("trans must be a square matrix")
  if (!identical(unique(diag(trans)),0)) warning("trans matrix contains diagonal entries which will be ignored.")
  diag(trans) <- 0
  if (type!="bespoke") {
  if (length(unique(c(trans)))!=max(trans)+1) {
    stop("Parameters in the trans matrix must be numbered consecutively from 1")
  }
  #Check the nonh matrix
  if (!is.null(nonh)) {
    if (!is.matrix(nonh)) stop("nonh must be a matrix")
    if (dim(nonh)[1]!=dim(nonh)[2]) stop("nonh must be a square matrix")
    if (!identical(unique(diag(nonh)),0)) warning("trans matrix contains diagonal entries which will be ignored.")
    diag(nonh) <- 0
    if (length(unique(c(nonh)))!=max(nonh)+1) {
      stop("Parameters in the nonh matrix must be numbered consecutively from 1")
    }
  }else{
    nonh <- array(0,dim(trans))
    type <- "gompertz"
    warning("No nonh matrix supplied. Assuming a time homogeneous model")
  }
  if (max(nonh)==0 & type!="gompertz") {
    type <- "gompertz"
    warning("nonh matrix implies a homogeneous model")
  }
  }else{
    if (!is.null(covm) | !is.null(nonh)) {
      warning("covm  and nonh not required for bespoke models. Argument(s) ignored.")
    }
  }
  if (is.null(covariates) & !is.null(covm)) warning("covm supplied but covariates argument missing: covm will be ignored")
  if (!is.logical(death)) stop("death must be TRUE/FALSE")
  if (death & is.null(death.states)) warning("No death.states specified for a model with exact death times.")
  if (!death & !is.null(death.states)) warning("death.states specified for a model without exact death times.")
  if (is.null(censor) & !is.null(censor.states)) stop("censor.states specified without specifying the censoring label(s)")
  if (!is.null(firstobs)) {
    if (is.null(emat)) {
      warning("firstobs only relevant for misclassification models")
    }else{
      if (!firstobs%in%c("absent","misc","exact")) stop("firstobs must be either \'absent\', \'misc\', \'exact\'.")
    }
  }
  initial_process <- process_inputs(formula,  data, covariates,covm,ecovm,initcovm)
  data <- initial_process$data
  covariates <- initial_process$covariates
  covm <- initial_process$covm
  ecovm <- initial_process$ecovm
  initcovm <- initial_process$initcovm

  #if (!"state"%in%names(data) | !"time"%in%names(data) | !"subject"%in%names(data)) {
  #  stop("data must contain 'state', 'time' and 'subject' columns")
  #}
  possiblestates <- c(1:dim(trans)[1],unlist(censor))
  if (any(!data$state%in%possiblestates)) stop("State vector includes values inconsistent with the model specification")
  if (any(!possiblestates%in%unique(data$state))) {
    absentstates <- paste(possiblestates[!possiblestates%in%unique(data$state)],collapse=", ")
    warning(paste("Some expected states not present in the data:",absentstates,sep=" "))
  }

  #Ensure censoring only occurs as last observation
  if (!is.null(censor) & sum(data$state%in%censor)>0) {
    sub <- tapply(1:length(data$subject),data$subject,max)
    cenmin <- tapply(1:length(data$subject),list(data$subject,1*(data$state %in% censor)),min)[,2]
    if (sum((sub - cenmin)[!is.na(cenmin)] !=0)>0) stop("Only final state for a subject can be censored")
  }


  if (!is.null(covariates)) {
    if (is.vector(covariates)) covariates <- cbind(covariates)
    if (dim(data)[1] != dim(covariates)[1]) stop("covariate not of same length as state/time vectors")
  }


  if (!is.null(emat)) {
    if (!identical(dim(trans),dim(emat))) stop("emat must be the same dimension as trans")
    if (!identical(unique(diag(emat)),0)) warning("emat matrix contains diagonal entries which will be ignored.")
    diag(emat) <- 0
    ematpars <- max(emat)
    if (length(unique(c(emat)))!=max(emat)+1) {
      stop("Parameters in the emat matrix must be numbered consecutively from 1")
    }
    if (!is.null(ecovm)) ematpars <- ematpars + max(ecovm)
    if (!is.null(initp)) ematpars <- ematpars + max(initp)
    if (!is.null(initcovm)) ematpars <- ematpars + max(initcovm)
  }else{
    if (!is.null(ecovm)) warning("ecovm only relevant to misclassification models. Term ignored. Add an emat term, if desired.")
    ematpars <- 0
  }

  if (!type%in%c("bespoke","weibull","gompertz","bspline")) {
    stop("Only types: bespoke, weibull, gompertz, bspline are currently supported.")
  }
  if (type!="bespoke") {
    if (is.null(trans) | is.null(nonh)) stop("trans and nonh must be supplied if not using a user supplied intens function")
    if (!identical(dim(trans),dim(nonh))) stop("trans and nonh should be the same dimension")
    if (!is.null(intens)) warning("Supplied intens ignored. Set type to bespoke if want a user supplied function")
    intens <- intens_generate.nhm(type=type,trans = trans,nonh = nonh,covm=covm,centre_time = centre_time,splinelist = splinelist,degrees=degrees, covnames=names(covariates))
    nparQ <- attr(intens,"npar")
  }else{
    nparQ <- attr(intens,"npar")
    if (is.null(nparQ)) {
      ###Attempt to infer the number of parameters
       nparQ <- dim(intens(0,rep(0,dim(covariates)[2]),rep(0,1000))$qp)[3]
      warning("Number of parameters in intensities model inferred from the dimension of the derivative matrix.")
    }
  }

  if (is.null(firstobs)) {
    firstobs <- "exact"
  }


  if (!is.null(emat)) {
    nparI <-0
    if (!is.null(initp)) {
      nparI <- max(initp)
      if (!is.null(initcovm)) nparI <- nparI + max(initcovm)
    }else{
      if (!is.null(initcovm)) warning("initcovm ignored: only relevant if initial probabilities estimated. Add an initp term, if required.")
    }
    emat_nhm <- emat_generate.nhm(emat, ecovm, censor, censor.states, death, death.states, intens, nparQ,nparI,covnames=names(covariates))
    nparE <- attr(emat_nhm,"npar")
    initp_nhm <- NULL
    if (is.null(initp)) {
      est_initp <- FALSE
    }else{
      est_initp <- (max(initp)>0)
    }
    if (is.null(initp_value)) {
      initp_value <- c(1,rep(0,dim(trans)[1]-1))
    }else{
      if (firstobs=="exact") warning("Supplied initp vector will be ignored. Specify firstobs=\"absent\" if unknown initial probabilities desired.")
      if (length(initp_value)!=dim(trans)[1]) stop("Supplied initial probability vector must match the number of states")
      if (min(initp_value) < 0) stop("Invalid initial probability vector")
      if (sum(initp_value)!=1) {
        warning("initp vector normalized to sum to 1")
        initp_value <- initp_value/sum(initp_value)
      }
    }
      if (est_initp) {
        initp_nhm <- initp_generate.nhm(initp, initcovm, nparQ,nparE,covnames=names(covariates))
        #Should a function still be made otherwise?
        nparI <- attr(initp_nhm,"npar")
      }else{
        nparI <- 0
      }
  }else{
    if (!is.null(initp) | !is.null(initcovm) | !is.null(initp_value)) warning("Initial probabilities only relevant to misclassification models.")
    nparI <- 0
    nparE <- 0
  }
  if (type=="bespoke") nstate <- dim(trans)[1]
  if (type=="bespoke" & is.null(attr(intens,"parnames"))) {
    attr(intens,"parnames")<-paste("Bespoke Q parameter",1:nparQ)
  }
  if (type=="bespoke" & is.null(attr(intens,"parclass"))) {
    attr(intens,"parclass")<-rep("Bespoke",nparQ)
  }
  #Convert the censor states to be nstate+1,nstate+2,...
  if (!is.null(censor)) {
    if (!is.null(censor.states)) {
      if (any(!unlist(censor.states)%in%(1:dim(trans)[1]))) stop("Invalid censor.states term")
    }
    for (i in 1:length(censor)) {
      data$state[data$state==censor[i]]<-i+dim(trans)[1]
    }
    censor <- dim(trans)[1] + 1:length(censor)

    if (!is.list(censor.states)) {
      censor.states <- list(censor.states)
    }
    if (length(censor.states)!=length(censor)) stop("censor.states must be a list of same length as the censor vector")
  }

  data$subject <- match(data$subject,unique(data$subject))

  output <- list()
  output$npar <- nparQ + nparE + nparI
  output$nparQ <- nparQ #Needed to check the intens function correctly.
  if (!is.null(covariates)) {
    output$ncov <- dim(covariates)[2]
  }else{
    output$ncov <- 0
  }
  output$nstate <- dim(trans)[1]
  output$state <- data$state
  output$time <- data$time
  output$subject <- data$subject
  output$covariates <- covariates
  output$type <- type
  output$censor <- censor
  output$censor.states <- censor.states
  output$death <- death
  output$death.states <- death.states
  output$hidden <- (!is.null(emat))
  output$firstobs <- firstobs
  output$initp <- initp_value
  output$intens <- intens
  output$trans <- trans
  output$parnames <- attr(intens,"parnames")
  output$parclass <- attr(intens,"parclass")
  if (output$hidden) {
    output$emat_nhm <- emat_nhm
    output$initp_nhm <- initp_nhm
    output$parnames <- c(output$parnames,attr(emat_nhm,"parnames"),attr(initp_nhm,"parnames"))
    output$parclass <- c(output$parclass,attr(emat_nhm,"parclass"),attr(initp_nhm,"parclass"))
  }
  class(output) <- "nhm_model"
  return(output)
}


nhm  <-  function(
  model_object, #nhm_model object created using model.nhm
  initial=NULL, #Vector of initial parameter values
  gen_inits=FALSE, #Whether to automatically generate initial values (only available for models without misclassification)
  control, #Object of class nhm.control specifying various settings for the solution of the KFEs and the optimization. Default is nhm.control(...)
  score_test=FALSE, #Whether a score test at the initial values should be performed.
  fixedpar=NULL
)
{
  if (!inherits(model_object,"nhm_model")) stop("model_object must be created using model.nhm")
  if (missing(control)) {
    control <- nhm.control(checks = (model_object$type=="bespoke"))
  }
  if (!inherits(control,"nhm_control")) stop("control object must be made using nhm.control")
  tmax <- control$tmax
  coarsen <- control$coarsen
  coarsen.vars <- control$coarsen.vars
  coarsen.lv <- control$coarsen.lv
  checks <- control$checks
  rtol <- control$rtol
  atol <- control$atol
  linesearch <- control$linesearch
  damped <- control$damped
  damppar <- control$damppar
  hessian <- control$obsinfo
  splits <- control$splits
  ncores <- control$ncores
  print.level <- control$print.level
  maxLikcontrol <- control$maxLikcontrol
  fishscore <- control$fishscore
  hessmeth <- ifelse(hessian,TRUE,"bhhh")
  npar <- model_object$npar
  nparQ <- model_object$nparQ
  ncov <- model_object$ncov
  nstate <- model_object$nstate
  state <- model_object$state
  time <- model_object$time
  subject <- model_object$subject
  covariates <- model_object$covariates
  censor <- model_object$censor
  censor.states <- model_object$censor.states
  death <- model_object$death
  death.states <- model_object$death.states
  hidden <- model_object$hidden
  intens <- model_object$intens
  emat_nhm <- model_object$emat_nhm
  initp <- model_object$initp
  firstobs <- model_object$firstobs
  initp_nhm <- model_object$initp_nhm
  if (is.null(initial) & hidden) stop("Initial values must be supplied for models with misclassification")
  if (is.null(initial) & !gen_inits) stop("Initial values must be supplied if gen_inits=FALSE")
  if (is.null(initial) & !is.null(fixedpar)) stop("Initial values must be supplied if fixedpar are specified")
  if (!is.null(fixedpar)) {
    if (gen_inits) stop("Initial values must be supplied if fixedpar are specified")
    if (max(fixedpar)>npar) stop("fixedpar includes values outside 1,2,...,npar")
    fixval <- initial[fixedpar]
  }else{
    fixval <- NULL
  }
  if (gen_inits) {
    initial <- generate_inits_nhm(data=data.frame(state=state,time=time,subject=subject),trans=model_object$trans,censor=censor, censor.states=censor.states,npar=npar)
    message("Initial values generated using crudeinits.msm, assuming time homogeneity and no covariate effects \n")
  }
  if (length(initial)!=npar) stop(paste("Initial parameter vector must be same length as number of parameters in the model. There should be",npar,"parameters corresponding to:",paste(model_object$parnames,collapse=", "),sep=" "))
  #Coarsen covariates if required:
  if (coarsen & ncov>0) {
    if (is.null(coarsen.vars)) coarsen.vars <- rep(1,ncov) #Which variables to coarsen
    if (is.null(coarsen.lv)) {
      coarsen.lv <- 10
      warning("Covariate coarsened to 10 unique values. Use coarsen.lv to choose a different number")
    }
    #Now perform the coarsening on the covariates vector...
    newcovariates <- covariates
    clust <- stats::kmeans(covariates[,coarsen.vars],coarsen.lv)
    covC <- clust$centers[clust$cluster,]
    newcovariates[,coarsen.vars] <- covC
    covariates <- newcovariates
  }

  #Check whether the trans matrix implies a progressive model
  if (is.null(model_object$trans)) {
    gensolve <- solve
    solvemethod <- "solve"
  }else{
    if (!identical(upper.tri(model_object$trans)*model_object$trans,model_object$trans)) {
      #Assume have a progressive model.
      gensolve <- solve
      solvemethod <- "solve"
    }else{
      gensolve <- function(x) backsolve(x, diag(1,dim(model_object$trans)[1]))
      solvemethod <- "backsolve"
    }
  }

  if (checks & !hidden) {
    #Sort data to ensure it is in correct form
    if (ncov>0) {
      fulldat <- data.frame(state=state,time=time,subject=subject,covariates=covariates)
    }else{
      fulldat <- data.frame(state=state,time=time,subject=subject)
    }
    fulldat <- fulldat[order(subject,time),]
    state <- fulldat$state
    time <- fulldat$time
    subject <- fulldat$subject
    if (ncov>0) covariates <- fulldat[,4:(3+ncov),drop=FALSE]

    #Check intens is consistent with values
    if (sum(is.na(initial))>0) stop("Missing or invalid entries in initial parameter vector")
    #print(initial)
    if (ncov>0) {
      check <- intens(1,rep(0,ncov),initial)
    }else{
      check <- intens(1,NULL,initial)
    }
    if (is.null(check$q) | is.null(check$qp)) stop("Generator matrix function intens not valid")
    if (!prod(dim(check$q)==nstate)) stop("Generator matrix function intens not valid. Number of states not consistent with the trans matrix")
    if (!prod(dim(check$qp)==c(nstate,nstate,nparQ))) stop("Derivatives in generator matrix function not valid. npar should match the number used in the generator matrix only. Parameters for ematrix and initp are excluded")
    if (is.null(death.states)) death.states <- nstate
    vs <- apply(check$q[death.states,,drop=FALSE],1,function(x) sum(abs(x)))
    if (death) {
      if (max(vs)>1e-8) stop("Exact death times can only be attributed to absorbing states!")
    }
  }


  processed <- dataprocess.nhm(state,time,subject,covariates,ncov,splits,firstobs) #Include firstobs in the data processing step.

  #Now need to pass across the final covariate values (at least if exact=TRUE)
  finalcovs <- processed$finalcovs
  if (hidden) {
    init_state <- processed$init_state
    initcovs <- processed$initcovs
    if (firstobs%in%c("exact","absent")) {
    fobs <- rep(tapply(1:length(subject),subject,min),table(subject))
    statel <- state[-fobs]
    timel <- time[-fobs]
    subjectl <- subject[-fobs]
    if (ncov>0) {
      cov2 <- covariates[-fobs,,drop=FALSE]
    }else{
      cov2 <- NULL
    }
    }else{
      statel <- state
      timel <- time
      subjectl <- subject
      if (ncov>0) {
        cov2 <- covariates
      }else{
        cov2 <- NULL
      }
    }
  }

  if (is.null(fishscore)) {
    fishscore <- FALSE
  }
  if (fishscore & (!is.null(censor) | death)) {
    fishscore <- FALSE
    warning("Fisher scoring not available for models with censoring or exact death times")
  }
  if (fishscore & !is.null(fixedpar)) {
    fishscore <- FALSE
    warning("Fisher scoring not currently available for models with fixedpar")
  }
  if (is.null(tmax)) tmax <- max(time)+1

  if (!is.null(fixedpar)) {
    initval <- initial[-fixedpar]
  }else{
    initval <- initial
  }

  if (!hidden) {
    if (fishscore) {
      fitmodel <- fisherscore.nhm(initial,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tmax,hessian=hessian,linesearch=linesearch,damped=damped,damppar=damppar,rtol=rtol,atol=atol,score_test=score_test,intens=intens,gensolve=gensolve,ncores=ncores)
      if (score_test) {
        fitmodel$parnames <- model_object$parnames
        fitmodel$parclass <- model_object$parclass
        fitmodel$par <- initial
        class(fitmodel) <- "nhm_score"
        return(fitmodel)
      }
    }else{
      #Note that this requires version 0.8 or above of maxLik to ensure the gradient is actually used in optimization.
      if (!score_test) fit <- maxLik(bhhh.nhm,start=initval,method="BHHH",print.level=print.level,finalHessian=hessmeth,control=maxLikcontrol,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tmax,fishscore=FALSE,indout=TRUE,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,finalcovs=finalcovs,gensolve=gensolve,ncores=ncores,fixedpar=fixedpar,fixval=fixval)
      if (score_test) {
        out <- bhhh.nhm(initval,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tmax,fishscore=FALSE,indout=TRUE,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,finalcovs=finalcovs,ncores=ncores,fixedpar=fixedpar,fixval=fixval)
        #Determine
        grad <- attr(out,"grad")
        g <- apply(grad,2,sum)
        I <- apply(array(apply(grad,1,function(x) outer(x,x)),c(length(g),length(g),processed$nsub)),c(1,2),sum)
        output <- list(g=g,I=I,parnames=model_object$parnames,parclass=model_object$parclass,par=initial,fixedpar=fixedpar)
        class(output) <-"nhm_score"
        return(output)
      }
      fitmodel <- list(par=c(fit$estimate),value=-fit$maximum,grad=fit$gradient,hess=-fit$hessian)
    }
  }else{
    #At the moment this is just the single value!
    if (!score_test) fit <- maxLik(bhhh.nhm,start=initval,method="BHHH",print.level=print.level,finalHessian=hessmeth,control=maxLikcontrol,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tmax,fishscore=FALSE,obslist=processed$obslist,state=statel,subject=subjectl,time=timel,indout=TRUE,init_state=init_state,cov2=cov2,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,emat_nhm=emat_nhm,ncores=ncores,firstobs=firstobs,initp=initp,initcovs=initcovs,finalcovs=finalcovs,initp_nhm=initp_nhm,nparQ=nparQ,fixedpar=fixedpar,fixval=fixval)
    if (score_test) {
      out <- bhhh.nhm(initval,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,tcrit=tmax,fishscore=FALSE,obslist=processed$obslist,state=statel,subject=subjectl,time=timel,init_state=init_state,indout=TRUE,cov2=cov2,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,emat_nhm=emat_nhm,ncores=ncores,firstobs=firstobs,initp=initp,initcovs=initcovs,finalcovs=finalcovs,initp_nhm=initp_nhm,nparQ=nparQ,fixedpar=fixedpar,fixval=fixval)
      grad <- attr(out,"grad")
      g <- apply(grad,2,sum)
      I <- apply(array(apply(grad,1,function(x) outer(x,x)),c(length(g),length(g),processed$nsub)),c(1,2),sum)
      output <- list(g=g,I=I,parnames=model_object$parnames,parclass=model_object$parclass,par=initial,fixedpar=fixedpar)
      class(output) <-"nhm_score"
      return(output)
    }
    #if (!score_test) fit <- maxLik(bhhh_misc.nhm,start=initial,method="BHHH",print.level=2,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,use.deriv=TRUE,tcrit=tmax,fishscore=FALSE,obslist=processed$obslist,state=statel,subject=subjectl,time=timel,indout=TRUE,init_state=init_state,cov2=cov2,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,emat_nhm=emat_nhm)
    #if (score_test) return(bhhh_misc.nhm(initial,npar=npar,ncov=ncov,nstate=nstate,ncovvals=processed$ncovvals,covs=processed$covs,covtable=processed$covtable,timesets=processed$timeset,fromtimeindex=processed$fromtimeindex,totimeindex=processed$totimeindex,fromtimesets=processed$fromtimesets,totimesets=processed$totimesets,fromsets=processed$fromsets,tosets=processed$tosets,death=death,death.states=death.states,censor=censor,censor.states=censor.states,use.deriv=TRUE,tcrit=tmax,fishscore=FALSE,obslist=processed$obslist,state=statel,subject=subjectl,time=timel,init_state=init_state,indout=TRUE,cov2=cov2,sublist=processed$sublist,nsub=processed$nsub,rtol=rtol,atol=atol,intens=intens,gensolve=gensolve,emat_nhm=emat_nhm))
    fitmodel <- list(par=c(fit$estimate),value=-fit$maximum,grad=fit$gradient,hess=-fit$hessian)

  }
  #Convert the fitted model into a proper standardized object?
  #Need enough information to directly apply plot.prevalence.msm
  if (!is.null(fitmodel$hess)) {
    pdcheck <- min(Re(eigen(fitmodel$hess,only.values=TRUE)$values))
    if (pdcheck < 1e-8) {
      fitmodel$singular<-TRUE
      warning("Hessian is singular or close to singular. Either optimization has not converged or the MLE is on the boundary of the parameter space.")
    }else{
      fitmodel$singular<-FALSE
    }
  }
  fitmodel$model_object <- model_object  #May as well pass everything across
  fitmodel$parnames <- model_object$parnames
  fitmodel$npar <- npar
  fitmodel$ncov <- ncov
  fitmodel$nstate <- nstate
  fitmodel$tcrit <- tmax
  fitmodel$intens <- intens
  fitmodel$splits <- splits
  fitmodel$maxtime <- max(time)
  fitmodel$solvemethod <- solvemethod
  fitmodel$fixedpar <- fixedpar
  fitmodel$fixval <- fixval
  #Also want the mean covariate values...
  if (ncov>0) {
    fitmodel$covmeans <- apply(covariates,2,mean)
  }
  class(fitmodel) <- "nhm"
  return(fitmodel)
}

nhm.control <- function(tmax=NULL, coarsen=FALSE, coarsen.vars=NULL, coarsen.lv=NULL, checks=FALSE,rtol=1e-6, atol=1e-6, fishscore=NULL, linesearch=FALSE, damped=FALSE, damppar=0,obsinfo=TRUE,splits=NULL,ncores=1,print.level=2,maxLikcontrol=NULL) {
  if (!is.logical(coarsen)) stop("coarsen must be a logical")
  if (coarsen) {
    if (is.null(coarsen.vars)) stop("coarsen.vars must be specified")
    if (is.null(coarsen.lv)) stop("coarsen.lv must be specified")
    if (!is.numeric(coarsen.lv)) stop("coarsen.lv must be a numeric value")
    if (length(coarsen.lv)>1) stop("coarsen.lv must be a single numeric value")
    if (coarsen.lv <= 0) stop("coarsen.lv must be a positive integer")
  }
  if (!is.logical(checks)) stop("checks must be a logical")
  if (!is.numeric(rtol)) stop("rtol must be a numeric value")
  if (!is.numeric(atol)) stop("atol must be a numeric value")
  if (!is.logical(linesearch)) stop("linesearch must be a logical")
  if (!is.logical(damped)) stop("linesearch must be a logical")
  if (!is.numeric(damppar)) stop("damppar must be a numeric value")
  if (damppar < 0) stop("damppar must be non-negative")
  if (!is.logical(obsinfo)) stop("obsinfo must be a logical")
  if (!is.null(splits)) {
    if (!is.numeric(splits)) stop("splits must be a numeric vector, if supplied")
  }
  if (!is.numeric(ncores)) stop("ncores must be numeric")
  if (length(ncores)>1) stop("ncores must be a single value")
  if (ncores <1) stop("ncores must be greater or equal to 1.")
  if (!is.null(maxLikcontrol)) {
    if (!is.list(maxLikcontrol)) stop("maxLikcontrol must be a list, if supplied")
  }
  if (!is.null(fishscore)) {
    if (!is.logical(fishscore)) stop("fishscore must be a logical")
  }else{
    fishscore <- FALSE
  }
  if (!is.numeric(print.level)) stop("print.level must be numeric")
  controll <- list(tmax=tmax, coarsen=coarsen, coarsen.vars=coarsen.vars, coarsen.lv=coarsen.lv, checks=checks,rtol=rtol, atol=atol, linesearch=linesearch, damped=damped, damppar=damppar,obsinfo=obsinfo,splits=splits,ncores=ncores,print.level=print.level,maxLikcontrol=maxLikcontrol,fishscore=fishscore)
  class(controll) <- "nhm_control"
  return(controll)
}
