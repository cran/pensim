opt1D <-
function(nsim=50,nprocessors=1,setpen="L1",cl=NULL,...){
  if (!(identical(setpen,"L1")|identical(setpen,"L2"))) stop("setpen must be L1 or L2")
  clusterIsSet <- "cluster" %in% class(cl)
  if(nprocessors>1 | clusterIsSet){
    if(!clusterIsSet){
      nprocessors <- as.integer(round(nprocessors))
      library(snow)
      cl <- makeCluster(nprocessors, type="SOCK")
    }
    myseed=round(2^32*runif(6)) ##rlecuyer wants a vector of six seeds according to the SNOW manual
    tmp <- try(clusterSetupRNG(cl,seed=myseed))
    if(class(tmp) == "try-error") warning("rlecuyer is not properly configured on your system; child nodes may not produce random numbers independently.  Debug using rlecuyer examples if you are concerned about this, or use leave-one-out cross-validation.")
    if(identical(setpen,"L1")){
      thisopt <- parLapply(cl,1:nsim,function(n,...){
        optL1(...)
      },...)
    }else{   ##if(identical(setpen,"L1")){
      thisopt <- parLapply(cl,1:nsim,function(n,...){
        optL2(...)
      },...)
    }
    if(!clusterIsSet){
      stopCluster(cl)
    }
  }else{  ##if(nprocessors>1){
    if(identical(setpen,"L1")){
      thisopt <- lapply(1:nsim,function(n,...) optL1(...),...)
    }else{
      thisopt <- lapply(1:nsim,function(n,...) optL2(...),...)
    }
  }
  output <- sapply(thisopt,function(x){
    coefs <- coefficients(x$fullfit,"all")
    tmp <- c(x$lambda,x$cvl,coefs)
    names(tmp) <- c(setpen,"cvl",names(coefs))
    return(tmp)
  })
  return(t(output))
}

