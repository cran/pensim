create.data <-
  function(nvars=c(100,100,100,100,600),cors=c(0.8,0,0.8,0,0),associations=c(0.5,0.5,0.3,0.3,0),firstonly=c(TRUE,FALSE,TRUE,FALSE,FALSE),nsamples=100,response="timetoevent",basehaz=0.2,logisticintercept=0){
  #check that the input variables are appropriate
  if(class(nvars)!="numeric") stop("nvars must be a numeric vector")
  if(class(cors)!="numeric") stop("cors must be a numeric vector")
  if(class(firstonly)!="logical") stop("firstonly must be a logical vector")
  if(class(associations)!="numeric") stop("associations must be a numeric vector")
  if(length(nvars)!=length(cors)|length(nvars)!=length(firstonly)|length(nvars)!=length(associations)) stop("nvars, cors, firstonly, and associations must all have the same length.")
  library(MASS)
  #create x.out, the matrix of predictor values.
  x.out <- matrix(0.0,ncol=sum(nvars),nrow=nsamples)
  definecors <- data.frame(start=c(1,cumsum(nvars[-length(nvars)])+1),
                           end=cumsum(nvars),
                           cors=cors,
                           associations=associations,
                           num=nvars,
                           firstonly=firstonly,
                           row.names=letters[1:length(nvars)])
  Sigma <- matrix(0.0,ncol=sum(nvars),nrow=sum(nvars))  #covariance matrix
  # wts will be the associations of each variable to the response,
  # defined by the input variable "associations".  Initialize all elements
  # to zero, then only non-zero elements will be set later.
  wts <- rep(0,sum(nvars))
  # This loop defines the wts according to the associations input variable, and creates
  # the predictor variables in x.out
  for (i in 1:nrow(definecors)){
    thisrange <- definecors[i,"start"]:definecors[i,"end"]
    Sigma[thisrange,thisrange] <- definecors[i,"cors"]
    diag(Sigma) <- 1
    # It seems a bit funny using mvrnorm separately for each set of
    # correlated variables, but for large sample sizes or many variables
    # this is faster than using the entire covariance matrix just once
    x.out[,thisrange] <- mvrnorm(n=nsamples,mu=rep(0,nvars[i]),Sigma=Sigma[thisrange,thisrange])
    # Create the weights that define the association between each variable and the outcome.
    # To make the outcome depend on the average of a group of variables, then the weights
    # are just divided by definecors[i,"num"]
    if(definecors[i,"firstonly"]){
     # set only the first in this group of coefficients:
      wts[definecors[i,"start"]] <- definecors[i,"associations"]
    }else{
    # set this entire group of coefficients:
      wts[definecors[i,"start"]:definecors[i,"end"]] <- definecors[i,"associations"]
    }
    varnames <- paste(letters[i],1:nvars[i],sep=".")
    names(wts)[definecors[i,"start"]:definecors[i,"end"]] <- varnames
  }
  names(wts) <- make.unique(names(wts)) #in case anyone ever uses more than 26 types of variables
  dimnames(Sigma) <- list(colnames=names(wts),rownames=names(wts))
  ###got rid of these two line in favour of using mvrnorm within the loop
  ###diag(Sigma) <- 1
  ###x.out <- mvrnorm(n=nsamples,mu=rep(0,sum(nvars)),Sigma=Sigma)
  colnames(x.out) <- names(wts)
  # calculate time to recurrence
  betaX <- x.out%*%wts
  # convert x.out to a dataframe to make adding new columns easy
  x.out <- data.frame(x.out)
  #calculate response
  if(identical(response,"timetoevent")){
    h=basehaz*exp(betaX[,1])  #betaX is a one-column matrix; use the first column
    x.out$time <- rexp(length(h),h)
    x.out$cens <- 1
  }else if(identical(response,"binary")){
    p <- 1/(1+exp(-(betaX+logisticintercept)))
    x.out$outcome <- rbinom(length(p),1,p)
  }else stop("response must be either timetoevent or binary")
  return(list(summary=definecors,associations=wts,covariance=Sigma,data=x.out))
}
