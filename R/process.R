#Preprocessing functions.

#Initial processing to get data, covariates, covm and ecovm into a standard form
process_inputs <- function(formula,  data, covariates,covm,ecovm,initcovm) {
  #Function to convert the inputs into a standardized form

  #NB: subject no longer does anything useful...

  #Specifically:
  #a) Identify the state, time and subject variables (core data)
  #a) Split up the covariates from the core data
  #b) Standardize the names of the core data
  #d) Turn covm and ecovm into standardized nstate x nstate x ncov arrays
  subjectvar <- "subject"
  statevar <- all.vars(formula[[2]])
  timevar <- all.vars(formula[[3]])

  indexes <- match(c(statevar,timevar,subjectvar),names(data), nomatch=0)
  if (any(indexes==0)) {
    k <- c(statevar,timevar,subjectvar)[which(indexes==0)]
    stop(paste("Variable(s) ",k,"not found in the data."))
  }
  standard_data <- data[,indexes]
  names(standard_data) <- c("state","time","subject")
  if (is.factor(standard_data$state)) stop("state variable should be numeric, not a factor")
  if (any(is.na(standard_data))) stop("state, time and subject variables must not contain missing values")
  
  if (!is.null(covariates)) {

    covariates <- unique(covariates)

    #Next look at covm and ecovm

    #Make covm into an array rather than list
    if (!is.null(covm)) {
      if (is.list(covm)) {
        if (is.null(names(covm))) stop("covm must be a named list")
        indexes <- match(names(covm),covariates,nomatch=0)
        if (any(indexes==0)) {
          stop("names in covm must correspond to named covariates")
        }
        standard_covm <- array(0,c(dim(covm[[1]]),length(covariates)))
        for (i in 1:length(covm)) {
          standard_covm[,,indexes[i]] <- covm[[i]]
        }
      }else{
        if (dim(covm)[3]!=length(covariates)) stop("If covm is an array it must be of dimension nstate x nstate x ncov")
        standard_covm <- covm
      }
    }else{
      standard_covm <- NULL
    }

    #Make ecovm into an array rather than list
    if (!is.null(ecovm)) {
      if (is.list(ecovm)) {
        if (is.null(names(ecovm))) stop("ecovm must be a named list")
        indexes <- match(names(ecovm),covariates,nomatch=0)
        if (any(indexes==0)) {
          stop("names in ecovm must correspond to named covariates")
        }
        standard_ecovm <- array(0,c(dim(ecovm[[1]]),length(covariates)))
        for (i in 1:length(ecovm)) {
          standard_ecovm[,,indexes[i]] <- ecovm[[i]]
        }
      }else{
        if (dim(ecovm)[3]!=length(covariates)) stop("If ecovm is an array it must be of dimension nstate x nstate x ncov")
        standard_ecovm <- ecovm
      }
    }else{
      standard_ecovm <- NULL
    }
    
    
    #Make ecovm into an array rather than list
    if (!is.null(initcovm)) {
      if (is.list(initcovm)) {
        if (is.null(names(initcovm))) stop("ecovm must be a named list")
        indexes <- match(names(initcovm),covariates,nomatch=0)
        if (any(indexes==0)) {
          stop("names in initcovm must correspond to named covariates")
        }
        standard_initcovm <- array(0,c(length(initcovm[[1]]),length(covariates)))
        for (i in 1:length(initcovm)) {
          if (!is.vector(initcovm[[i]])) stop("initcovm should be a list of vectors of length nstate")
          standard_initcovm[,indexes[i]] <- initcovm[[i]]
        }
      }else{
        if (dim(initcovm)[2]!=length(covariates)) stop("If initcovm is an array it must be of dimension nstate x ncov")
        standard_initcovm <- initcovm
      }
    }else{
      standard_initcovm <- NULL
    }

    #Check the supplied covariates are present and are not the time, state or subject variables
    covmatch <- match(covariates,names(data),nomatch = 0)
    if (any(covmatch==0)) stop("covariates must match names in data frame")
    if (statevar%in%covariates) warning("state variable assigned as covariate within the model. Unlikely to be a sensible model.")
    if (timevar%in%covariates) warning("time variable assigned as covariate within the model. Unlikely to be a sensible model.")
    if (subjectvar%in%covariates) warning("subject variable assigned as covariate within the model. Unlikely to be a sensible model.")

    standard_covariates <- data[,covmatch,drop=FALSE]
    whichfact <- sapply(1:length(covmatch), function(x) is.factor(standard_covariates[,x]))
    if (any(whichfact)) stop("Factor covariates not currently supported. Create appropriate binary dummy variables.")
    
    
    ###Check for missing covariate values
    if (any(is.na(standard_covariates))) {
      #Establish which rows have missing values
      whichmiss <- which(apply(standard_covariates,1,function(x) any(is.na(x))))
      #Determine which subjects affected
      affect_sub <- unique(standard_data$subject[whichmiss])
      #Delete these patients and renumber
      todelete <- which(standard_data$subject %in% affect_sub)
      ndeleted <- length(affect_sub)
      nrowdeleted <- length(todelete)
      string <- paste("Covariates include missing values.",ndeleted,"subjects have been deleted, corresponding to",nrowdeleted,"observations",sep=" ")
      warning(string)
      standard_data <- standard_data[-todelete,]
      standard_covariates <- standard_covariates[-todelete,]
      #Renumber the subjects so they are in numerical order
      standard_data$subject <- match(standard_data$subject,unique(standard_data$subject))
    }
    ###########
    
  }else{
    standard_covariates <- standard_covm <- standard_ecovm <- standard_initcovm <- NULL
  }
  processed <- list(data=standard_data,covariates=standard_covariates, covm=standard_covm, ecovm=standard_ecovm,initcovm=standard_initcovm)
  return(processed)
}

#Data processing for finding the transition probabilities
dataprocess.nhm <- function(state,time,subject,covariates,ncov,splits,firstobs) {
  #Function to process data for speedier maximization of the likelihood.
  fobs <- rep(tapply(1:length(subject),subject,min),table(subject))
  lobs <- rep(tapply(1:length(subject),subject,max),table(subject))
 if (firstobs%in%c("exact","absent")) {
  #Essentially discard the first observation as a placeholder
  from <- c(NA,state[1:(length(state)-1)])
  init_state <- state[unique(fobs)]
  to <- state
  from <- from[-fobs]
  to <- to[-fobs]
  fromt <- c(NA,time[1:(length(state)-1)])
  tot <- time
  fromt <- fromt[-fobs]
  tot <- tot[-fobs]
  subject <- subject[-fobs]
 }else{
   ftim <- rep(tapply(time,subject,min),table(subject))
   from <- c(1,state[1:(length(state)-1)])
   init_state <- state[unique(fobs)] #NB: This won't be used if "exact" is not specified.
   to <- state
   from[unique(fobs)] <- 1 #NB: This will effectively be treated like a placeholder
   fromt <- c(0,time[1:(length(state)-1)])
   fromt[unique(fobs)] <- ftim[unique(fobs)] #Set to the minimum time so that the first pmat is the identity matrix
   tot <- time
 }
  subject <- match(subject,unique(subject)) #Force to be unique 1,2,...,N
  nsub <- length(unique(subject))

  if (!is.null(splits)) {
    splits <- sort(unique(c(min(fromt),splits,max(fromt))))
    splitcat <- as.numeric(cut(fromt, breaks=splits,include.lowest = TRUE))
  }else{
    splitcat <- rep(1,length(fromt))
  }

  if (ncov>0) {
    initcovs <- covariates[unique(fobs),,drop=FALSE] #NB: Not used in the "exact" case and only if there are initp covariates.
    finalcovs<- covariates[unique(lobs),,drop=FALSE]
    if (firstobs%in%c("exact","absent")) {
     rem<-(1:length(covariates[,1]))[-fobs]-1
     covariates <- covariates[rem,,drop=FALSE]
    }else{
      #Need to retain but repeat the first ones twice?
      rem<-(1:length(covariates[,1]))[-fobs]-1
      covariates0 <- covariates[rem,,drop=FALSE]
      covariates <- 0*covariates
      covariates[-unique(fobs),,drop=FALSE] <- covariates0
      covariates[unique(fobs),,drop=FALSE] <- initcovs 
    }
    uniq  <-  data.frame(unique(cbind(covariates,splitcat)))
    uniqcov <- uniq[,-dim(uniq)[2],drop=FALSE]
    nouniq  <-  dim(uniq)[1]
    pastedu  <-  do.call("paste",uniq)
    pastedc  <-  do.call("paste",data.frame(cbind(covariates,splitcat)))
    covindex  <-  match(pastedc,pastedu)
  }else{
    finalcovs <- initcovs <- uniqcov <- uniq <- NULL
    covindex <- splitcat
  }

  timeslist <- tapply(c(fromt,tot),rep(covindex,2),function(x) sort(unique(x)))
  fromlist <- tapply(from,covindex,list)
  obslist <- tapply(1:length(from),covindex,list)
  tolist <- tapply(to,covindex,list)
  fromtlist <- tapply(fromt,covindex,list)
  totlist <- tapply(tot,covindex,list)
  sublist <- tapply(subject,covindex,list)
  if (length(unique(covindex))>1) {
    covtable <- table(covindex)
    ncovvals <- length(covtable)
  }else{
    covvs <- NULL
    covtable <- length(from)
    ncovvals <- 1
  }
  fromtindex <- fromtlist
  totindex <- totlist
  subjectindex <- sublist
  for (i in 1:ncovvals) {
    fromtindex[[i]] <- apply(array(fromtlist[[i]],c(length(fromtlist[[i]]),1)),1,function(x) which(x==timeslist[[i]]))
    totindex[[i]] <- apply(array(totlist[[i]],c(length(totlist[[i]]),1)),1,function(x) which(x==timeslist[[i]]))
    subjectindex[[i]] <- apply(array(sublist[[i]],c(length(sublist[[i]]),1)),1,function(x) which(x==timeslist[[i]]))
  }
  list(ncovvals=ncovvals,covs=uniqcov,covtable=covtable,timesets=timeslist,fromtimeindex=fromtindex,totimeindex=totindex,fromtimesets=fromtlist,totimesets=totlist,fromsets=fromlist,tosets=tolist,sublist=sublist,nsub=nsub,obslist=obslist,init_state=init_state,initcovs=initcovs,finalcovs=finalcovs)
}




get_names <- function(trans,nonh,covm,covnames,model,nparper=NULL) {
  covlabel <- "Covariate:"
  if (model=="gompertz") {
    label <- c("Base:","Trend:")
  }
  if (model=="weibull") {
    label <- c("Rate:","Shape:")
  }
  if (model=="bspline") {
    label <- c("Base:","Spline pars:")
    if (is.null(nparper)) stop("nparper missing")
  }
  if (model=="emat") {
    label <- c("Emat:")
    covlabel <- "Emat Covariate:"
  }
  R <- dim(trans)[1]
  if (model=="initp") {
    label <- c("Initp:")
    covlabel <- "Initp Covariate:"
    R <- length(trans)
  }
  
  if (model!="initp") {
  n1 <- (array(paste(rep(1:R,R),rep(1:R,each=R),sep="->"),c(R,R)))
  parnamesA <- paste(label[1],sapply(1:max(trans),function(x) paste(c(n1)[which(c(trans)==x)],collapse="/")),sep=" ")
  }else{
  n1 <- 1:R
  parnamesA <- paste(label[1],sapply(1:max(trans),function(x) paste(c(n1)[which(c(trans)==x)],collapse="/")),sep=" ")
  }
  if (!model%in%c("emat","initp")) {
    if (max(nonh)>0) {
    n2 <- (array(paste(rep(1:R,R),rep(1:R,each=R),sep="->"),c(R,R)))
    parnamesB <- paste(label[2],sapply(1:max(nonh),function(x) paste(c(n2)[which(c(nonh)==x)],collapse="/")),sep=" ")
    }else{
      parnamesB <- NULL
    }
  }else{
    parnamesB<-NULL
  }
  if (model=="bspline") {
    parnamesBb <- unlist(sapply(nparper,function(x) 1:x))
    parnamesB <- paste(rep(parnamesB,nparper),parnamesBb,sep="...")
  }
  
  if (model!="initp") {
  if (!is.null(covm)) {
    Nc <- dim(covm)[3]

    if (is.null(covnames)) covnames <- paste(covlabel,1:Nc)

    n3 <- (array(paste(covlabel,rep(covnames,each=R^2),rep(rep(1:R,R),Nc),"->",rep(rep(1:R,each=R),Nc)),c(R,R,Nc)))
    parnamesC <- paste(sapply(1:max(covm),function(x) paste(c(n3)[which(c(covm)==x)],collapse="/")),sep=" ")
  }else{
    parnamesC <- NULL
  }
  }else{
    if (!is.null(covm)) {
    Nc <- dim(covm)[2]
    R <- dim(covm)[1]
    if (is.null(covnames)) covnames <- paste(covlabel,1:Nc)
    nc <- array(paste(covlabel,rep(covnames,each=R),"on",rep(1:R,Nc)),c(R,Nc))
    parnamesC <- paste(sapply(1:max(covm),function(x) paste(c(nc)[which(c(covm)==x)],collapse="/")),sep=" ")
    }else{
      parnamesC <-NULL
    }
    
  }
  parnames <- c(parnamesA,parnamesB,parnamesC)
  if (!model%in%c("emat","initp")) {
    parclass <- rep(c("Trans","Nonhom","Cov"),c(length(parnamesA),length(parnamesB),length(parnamesC)))
  }else{
    if (model=="emat") {
    parclass <- rep("Emat",length(parnames))
    }else{
    parclass <- rep("Initp",length(parnames))  
    }
  }
  return(list(parnames=parnames,parclass=parclass))
}

#Create an overall wrapper function
intens_generate.nhm <- function(type="gompertz",trans, nonh, covm=NULL, centre_time=0, splinelist=NULL, degrees=NULL,covnames=NULL) {
  #TO DO:
  #Way of incorporating covariate values relating to ematrix for misclassification
  #Minimally would just specify how many extra parameters come from ematrix.

  #Validity checks:
  if (is.null(centre_time)) centre_time<-0
  if (!type%in%c("gompertz","weibull","bspline")) stop("Invalid type")
  if (dim(trans)[1]!=dim(trans)[2] | length(dim(trans))!=2) stop("trans should be a square matrix")
  if (dim(nonh)[1]!=dim(nonh)[2] | length(dim(nonh))!=2) stop("trans should be a square matrix")
  if (!identical(dim(trans),dim(nonh))) stop("trans and nonh should be of the same dimension")
  if (!is.null(covm) & (!identical(dim(covm)[1:2],dim(trans))  | length(dim(covm))!=3)) stop("covm should be an R x R x Nc array")
  if (type!="bspline" & (!is.null(splinelist) | !is.null(degrees))) warning("splinelist and degrees are only for bsplines")
  if (type!="gompertz" & centre_time!=0) warning("centre_time is only used for gompertz")

  if (type=="gompertz") {
    intens <- intens_generate_gompertz(trans,nonh,covm,centre_time,covnames)
  }
  if (type=="weibull") {
    intens <- intens_generate_weibull(trans,nonh,covm,covnames)
  }
  if (type=="bspline") {
    intens <- intens_generate_bspline(trans,nonh,covm,splinelist,degrees,covnames)
  }
  return(intens)
}



intens_generate_gompertz <- function(trans,nonh, covm=NULL,centre_time=0,covnames=NULL) {
  #trans: R x R matrix of viable transitions (0= no transition, p = which of the parameter vector)
  #nonh: R x R matrix of viable transitions with a linear time effect
  #covm: R x R x Nc array ofviable transitions with a log-linear covariate effect
  #centre_time: Value to centre_time the times by: potentially improves convergence.
  gtname <- get_names(trans,nonh,covm,covnames,model="gompertz")
  parnames <- gtname$parnames
  parclass <- gtname$parclass

  #approach would be to first establish what the mapping from x to Q is
  npar1 <- max(trans)
  npar2 <- max(nonh)
  
  if (!is.null(covm)) {
    npar3 <- max(covm)
  }else{
    npar3 <- 0
  }
  npar <- npar1 + npar2 + npar3
  B <- array(0,c(dim(trans),npar1))
  for (b in 1:npar1) {
    B[,,b] <- array(1,dim(trans))*(trans==b)
  }
  if (npar2 >0 ) {
  Ct <- array(0,c(dim(trans),npar2))
  for (b in 1:npar2) {
    Ct[,,b] <- array(1,dim(nonh))*(nonh==b)
  }
  }
  if (npar3 > 0) {
    C <- array(0,c(dim(covm),npar3))
    for (b in 1:npar3) {
      C[,,,b] <-  array(1,dim(covm))*(covm==b)
    }
  }

  #mapping

  #For each basic parameter, we have a matrix of

  outfun <- function(t,z,x) {
    t <- t - centre_time
    lQ <- Q <- array(0,dim(trans))
    DlQ <- DQ <- array(0,c(dim(trans),npar))
    DlQ[,,1:npar1] <- B
    for (b in 1:npar1) {
      lQ <- lQ + B[,,b]*x[b]
    }
    if (npar2 >0) {
    for (b in 1:npar2) {
      lQ <- lQ + Ct[,,b]*t*x[b+npar1]
      DlQ[,,npar1 + b] <- Ct[,,b]*t
    }
    }
    if (npar3 >0) {
      for (v in 1:length(z)) {
        for (b in 1:npar3) {
          lQ <- lQ + C[,,v,b]*z[v]*x[b+npar1+npar2]
          DlQ[,,npar1+npar2+b] <- DlQ[,,npar1+npar2+b] + C[,,v,b]*z[v]
        }
      }
    }
    Q <- exp(lQ)*(trans>0)
    diag(Q) <- -apply(Q,1,sum)
    for (b in 1:npar) {
      DQ[,,b] <- DlQ[,,b] * Q
      diag(DQ[,,b]) <- -apply(DQ[,,b],1,sum)
    }
    Q<-list(q=Q,qp=DQ)
    return(Q)
  }
  attr(outfun,"parclass") <- parclass
  attr(outfun,"parnames") <- parnames
  attr(outfun,"npar") <- npar #Number of parameters from the qmat model only.
  return(outfun)
}


#Version for Weibull models:

intens_generate_weibull <- function(trans,nonh, covm=NULL,covnames=NULL) {
  #trans: R x R matrix of viable transitions (0= no transition, p = which of the parameter vector)
  #nonh: R x R matrix of viable transitions with a Weibull shape parameter different from 1
  #covm: R x R x Nc array ofviable transitions with a log-linear covariate effect

  gtname <- get_names(trans,nonh,covm,covnames,model="weibull")
  parnames <- gtname$parnames
  parclass <- gtname$parclass

  #approach would be to first establish what the mapping from x to Q is
  npar1 <- max(trans)
  npar2 <- max(nonh)
  if (!is.null(covm)) {
    npar3 <- max(covm)
  }else{
    npar3 <- 0
  }
  npar <- npar1 + npar2 + npar3
  B <- array(0,c(dim(trans),npar1))
  for (b in 1:npar1) {
    B[,,b] <- array(1,dim(trans))*(trans==b)
  }
  Ct <- array(0,c(dim(trans),npar2))
  for (b in 1:npar2) {
    Ct[,,b] <- array(1,dim(nonh))*(nonh==b)
  }
  if (npar3 > 0) {
    C <- array(0,c(dim(covm),npar3))
    for (b in 1:npar3) {
      C[,,,b] <-  array(1,dim(covm))*(covm==b)
    }
  }

  #mapping

  #For each basic parameter, we have a matrix of

  outfun <- function(t,z,x) {
    if (t==0) t<-1e-6 #Avoid possible singularities at zero.
    lQ <- Q <- array(0,dim(trans))
    DlQ <- DQ <- array(0,c(dim(trans),npar))
    DlQ[,,1:npar1] <- B
    A <- L <- array(0,dim(trans))
    for (b in 1:npar1) {
      L <- L + B[,,b]*x[b]
    }
    for (b in 1:npar2) {
      A <- A + Ct[,,b]*x[b + npar1]
    }
    lQ <- A + exp(A)*L + (exp(A)-1)*log(t)
    for (b in 1:npar1) {
      DlQ[,,b] <- exp(A)*B[,,b]
    }
    for (b in 1:npar2) {
      DlQ[,,npar1 + b] <- Ct[,,b]*(1 + exp(A)*(L + log(t)))
    } 
    if (npar3 >0) {
      for (v in 1:length(z)) {
        for (b in 1:npar3) {
          lQ <- lQ + C[,,v,b]*z[v]*x[b+npar1+npar2]
          DlQ[,,npar1+npar2+b] <- DlQ[,,npar1+npar2+b] + C[,,v,b]*z[v]
        }
      }
    }
    Q <- exp(lQ)*(trans>0)
    diag(Q) <- -apply(Q,1,sum)
    for (b in 1:npar) {
      DQ[,,b] <- DlQ[,,b] * Q
      diag(DQ[,,b]) <- -apply(DQ[,,b],1,sum)
    }
    Q<-list(q=Q,qp=DQ)
    return(Q)
  }
  attr(outfun,"parclass") <- parclass
  attr(outfun,"parnames") <- parnames
  attr(outfun,"npar") <- npar #Number of parameters from the qmat model only.
  return(outfun)
}

#Version for B-spline basis models
intens_generate_bspline <- function(trans,nonh, covm=NULL, splinelist=NULL,degrees=NULL,covnames=NULL) {
  #As with the other functions but in addition:
  #splinelist: list of knot locations for spline basis functions including boundary knots
  #NB: If a value above the upper spline point is specified then take the upper knot point.
  #degrees: the degree of each spline. Defaults to 3 if not specified
  if (length(splinelist)!=max(nonh)) stop("Number of non-homogeneous transitions does not match splinelist")
  nspline <- max(nonh)
  if (is.null(degrees)) degrees <- rep(3,nspline)
  #Firstly determine the knots and boundary knots from each of the splines
  maxtimes <- sapply(splinelist,max)
  boundaryknots <- lapply(splinelist,function(x) c(min(x),max(x)))
  internalknots <- lapply(splinelist,function(x) x[-c(which.min(x),which.max(x))])
  #Also establish the required number of parameters
  nparper <- sapply(splinelist,length)+(degrees-2)
  transind <- rep(1:length(nparper),nparper)

  gtname <- get_names(trans,nonh,covm,covnames,model="bspline",nparper=nparper)
  parnames <- gtname$parnames
  parclass <- gtname$parclass

  npar1 <- max(trans)
  npar2 <- sum(nparper)
  if (!is.null(covm)) {
    npar3 <- max(covm)
  }else{
    npar3 <- 0
  }
  npar <- npar1 + npar2 + npar3
  B <- array(0,c(dim(trans),npar1))
  for (b in 1:npar1) {
    B[,,b] <- array(1,dim(trans))*(trans==b)
  }

  #This is the difficult part:
  Ct <- array(0,c(dim(trans),npar2))
  for (b in 1:npar2) {
    Ct[,,b] <- array(1,dim(nonh))*(nonh==transind[b])
  }

  if (npar3 > 0) {
    C <- array(0,c(dim(covm),npar3))
    for (b in 1:npar3) {
      C[,,,b] <-  array(1,dim(covm))*(covm==b)
    }
  }

  outfun <- function(t,z,x) {
    #Establish the spline basis points and put them into a vector length npar2
    u <- unlist(sapply(1:nspline,function(v) bs(x=min(t,maxtimes[v]),knots=internalknots[[v]],Boundary.knots = boundaryknots[[v]],degree=degrees[v])))
    lQ <- Q <- array(0,dim(trans))
    DlQ <- DQ <- array(0,c(dim(trans),npar))
    DlQ[,,1:npar1] <- B
    for (b in 1:npar1) {
      lQ <- lQ + B[,,b]*x[b]
    }
    for (b in 1:npar2) {
      lQ <- lQ + Ct[,,b]*u[b]*x[b+npar1]
      DlQ[,,npar1 + b] <- Ct[,,b]*u[b]
    }
    if (npar3 >0) {
      for (v in 1:length(z)) {
        for (b in 1:npar3) {
          lQ <- lQ + C[,,v,b]*z[v]*x[b+npar1+npar2]
          DlQ[,,npar1+npar2+b] <- DlQ[,,npar1+npar2+b] + C[,,v,b]*z[v]
        }
      }
    }
    Q <- exp(lQ)*(trans>0)
    diag(Q) <- -apply(Q,1,sum)
    for (b in 1:npar) {
      DQ[,,b] <- DlQ[,,b] * Q
      diag(DQ[,,b]) <- -apply(DQ[,,b],1,sum)
    }
    Q<-list(q=Q,qp=DQ)
    return(Q)
  }
  attr(outfun,"parclass") <- parclass
  attr(outfun,"parnames") <- parnames
  attr(outfun,"npar") <- npar #Number of parameters from the qmat model only.
  return(outfun)
}

#Version that will return arrays of times
emat_generate.nhm <- function(emat, ecovm=NULL, censor=NULL, censor.states=NULL, death=NULL, death.states=NULL, intens, nparQ,nparI,covnames=NULL) {
  npar1 <- max(emat)
  if (!is.null(ecovm)) {
    npar2 <- max(ecovm)
  }else{
    npar2 <-0
  }

  gtname <- get_names(emat,nonh=NULL,covm=ecovm,covnames,model="emat")
  parnames <- gtname$parnames
  parclass <- gtname$parclass

  #Main part follows in a similar way to Qmat case:
  B <- array(0,c(dim(emat),npar1))
  for (b in 1:npar1) {
    B[,,b] <- array(1,dim(emat))*(emat==b)
  }
  if (npar2 > 0) {
    C <- array(0,c(dim(ecovm),npar2))
    for (b in 1:npar2) {
      C[,,,b] <-  array(1,dim(ecovm))*(ecovm==b)
    }
  }
  ematI <- 1*(emat!=0)
  diag(ematI) <- 1

  outfun <- function(state,t,z,x,intens,override=FALSE) {
    if (is.vector(z)) z<-matrix(z,c(1,length(z)))
    if (override) death<-FALSE
    e <- array(0,c(dim(emat)[1],length(state)))
    de <- array(0,c(dim(emat)[1],length(state),nparQ + npar1+npar2+nparI))
    EI <- lE <- E <- array(0,c(dim(emat)[1],dim(emat)[2],length(state)))
    DlE <- DE <- array(0,c(dim(emat)[1],dim(emat)[2],nparQ + npar1+npar2+nparI,length(state)))

    E0 <- array(0,c(dim(emat)[1],dim(emat)[2] + length(censor),length(state)))
    DE0 <- array(0,c(dim(emat)[1],dim(emat)[2]+ length(censor),nparQ + npar1+npar2+nparI,length(state)))
    for (k in 1:length(state)) {
      DlE[,,(1:npar1) + nparQ,k] <- B
      for (b in 1:npar1) {
        lE[,,k] <- lE[,,k] + B[,,b]*x[b+nparQ]
      }
      if (npar2 >0) {
        for (v in 1:dim(z)[2]) {
          for (b in 1:npar2) {
            lE[,,k] <- lE[,,k] + C[,,v,b]*z[k,v]*x[b+npar1+nparQ]
            DlE[,,npar1+nparQ+b,k] <- DlE[,,npar1+nparQ+b,k] + C[,,v,b]*z[k,v]
          }
        }
      }
      EI[,,k] <- exp(lE[,,k])*ematI
      Enorm <- array(rep(apply(EI[,,k],1,sum),dim(emat)[1]),dim(emat))
      E[,,k] <- EI[,,k]/Enorm
      for (b in 1:(npar1 + npar2)) {
        Dlnorm <- array(rep(apply(DlE[,,b+nparQ,k]*EI[,,k],1,sum),dim(emat)[1]),dim(emat))
        DE[,,b+nparQ,k] <- E[,,k]*(DlE[,,b+nparQ,k] - Dlnorm/Enorm)
      }
      #Now expand to include death and censoring...

      E0[1:dim(emat)[1], 1:dim(emat)[2],k] <- E[,,k]
      DE0[1:dim(emat)[1], 1:dim(emat)[2],,k] <- DE[,,,k]

      if (death) {
        qintens <- intens(t[k],unlist(z[k,]),x)
        #First the absorbing states
        for (i in death.states) {
          E0[,i,k] <- qintens$q[,i]
          DE0[,i,1:nparQ,k] <- qintens$qp[,i,]
        }
      }
      if (length(censor)>0) {
        for (i in 1:length(censor)) {
          E0[censor.states[[i]],i+dim(emat)[2],k] <- 1
        }
      }
      e[,k] <- E0[,state[k],k]
      de[,k,] <- DE0[,state[k],,k]
    }
    E <- list(e=e,de=de)
    return(E)
  }
  attr(outfun,"parclass") <- parclass
  attr(outfun,"parnames") <- parnames
  attr(outfun,"npar") <- npar1 + npar2
  return(outfun)
}



#Version that will return arrays of times
initp_generate.nhm <- function(initp, initcovm=NULL, nparQ, nparE, covnames=NULL) {
  npar1 <- max(initp) #Same convention as for emat? or force initp[1]=0?
  if (!is.null(initcovm)) {
    npar2 <- max(initcovm)
  }else{
    npar2 <-0
  }
  
  #Need to set up a name generation approach for initp...
  gtname <- get_names(initp,nonh=NULL,covm=initcovm,covnames,model="initp")
  
  parnames <- gtname$parnames
  parclass <- gtname$parclass
  
  #Main part follows in a similar way to Qmat case:
  B <- array(0,c(length(initp),npar1))
  for (b in 1:npar1) {
    B[,b] <- 1*(initp==b)
  }
  if (npar2 > 0) {
    C <- array(0,c(dim(initcovm),npar2))
    for (b in 1:npar2) {
      C[,,b] <-  1*(initcovm==b)
    }
  }
  initpI <- 1*(initp!=0)
  initpI[1] <- 1
  
  outfun <- function(nsub,z,x) {
    if (is.vector(z)) z<-matrix(z,c(1,length(z)))
    EI <- lE <- E <- ip <- array(0,c(length(initp),nsub))
    #For some reason this is the otherway around from elsewhere...
    DlE <- DE <- dip <- array(0,c(length(initp),nparQ + nparE + npar1+npar2,nsub))
    #######
    for (k in 1:nsub) {
      DlE[,(1:npar1) + nparQ+nparE,k] <- B
      for (b in 1:npar1) {
        lE[,k] <- lE[,k] + B[,b]*x[b+nparQ+nparE]
      }
      if (npar2 >0) {
        for (v in 1:dim(z)[2]) {
          for (b in 1:npar2) {
            lE[,k] <- lE[,k] + C[,v,b]*z[k,v]*x[b+npar1+nparQ+nparE]
            DlE[,npar1+nparQ+nparE+b,k] <- DlE[,npar1+nparQ+nparE+b,k] + C[,v,b]*z[k,v]
          }
        }
      }
      EI[,k] <- exp(lE[,k])*initpI
      Enorm <- sum(EI[,k])
      E[,k] <- EI[,k]/Enorm
      for (b in 1:(npar1 + npar2)) {
        Dlnorm <- sum(DlE[,b+nparQ+nparE,k]*EI[,k])
        DE[,b+nparQ+nparE,k] <- E[,k]*(DlE[,b+nparQ+nparE,k] - Dlnorm/Enorm)
      }
    }
    I <- list(initp=E,dinitp=DE)
    return(I)
  }
  attr(outfun,"parclass") <- parclass
  attr(outfun,"parnames") <- parnames
  attr(outfun,"npar") <- npar1 + npar2
  return(outfun)
}




#Come up with a way of generating initial values for a model without misclassification...

generate_inits_nhm <- function(data,trans,censor=NULL, censor.states=NULL,npar) {
  crudevals <- msm::crudeinits.msm(state~time, subject=data$subject, qmatrix=trans,data=data, censor=censor,censor.states=NULL)
  par_ests <- log(tapply(c(crudevals)[c(trans)!=0],c(trans)[c(trans)!=0],mean))
  pars <- rep(0,npar)
  pars[1:max(trans)]<-par_ests
  return(pars)
}
