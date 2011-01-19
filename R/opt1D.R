opt1D <-
function(nsim=50,nprocessors=1,setpen="L1",...){
  if (!(identical(setpen,"L1")|identical(setpen,"L2"))) stop("setpen must be L1 or L2")
  if(nprocessors>1){
    nprocessors <- as.integer(round(nprocessors))
    library(snow)
    library(rlecuyer)
    cl <- makeCluster(nprocessors, type="SOCK")
    myseed=round(2^32*runif(6)) ##rlecuyer wants a vector of six seeds according to the SNOW manual
    clusterSetupRNG(cl,seed=myseed)
    if(identical(setpen,"L1")){
      thisopt <- parLapply(cl,1:nsim,function(n,...){
        library(penalized)
        optL1(...)
      },...)
    }else{   ##if(identical(setpen,"L1")){
      thisopt <- parLapply(cl,1:nsim,function(n,...){
        library(penalized)
        optL2(...)
      },...)
###      thisopt <- parLapply(cl,1:nsim,function(n,...) optL2(...),...)
    }
    stopCluster(cl)
  }else{  ##if(nprocessors>1){
    library(penalized)
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

