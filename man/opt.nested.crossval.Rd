\name{opt.nested.crossval}
\alias{opt.nested.crossval}
\title{Parallelized calculation of cross-validated risk score
  predictions from L1/L2/Elastic Net penalized regression.
}
\description{
  calculates risk score predictions by a nested cross-validation, using
  the optL1 and optL2 functions of the penalized R package for
  regression.  In the outer level of cross-validation, samples are split
  into training and test samples.  Model parameters are tuned by cross-validation within
  training samples only.

  By setting nprocessors > 1, the outer cross-validation is split between multiple processors.

  The functions support z-score scaling of training data, and application
  of these scaling and shifting coefficients to the test data.  It also
  supports repeated tuning of the penalty parameters and selection of the
  model with greatest cross-validated likelihood.
}
\usage{
opt.nested.crossval(outerfold=10,nprocessors=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outerfold}{
    number of folds in outer cross-validation (the level used for validation)
  }
  \item{nprocessors}{
    An integer number of processors to use.  If specified in 
    opt.nested.crossval, iterations of the outer cross-validation are
    sent to different processors.  If specified in opt.splitval,
    repeated starts for the penalty tuning are sent to different processors.
  }
  \item{\dots}{
    optFUN (either "opt1D" or "opt2D"), scaling (TRUE to z-score
    training data then apply the same shift and scale factors to test
    data, FALSE for no scaling) are passed onto the opt.splitval
    function.  Additional arguments are required, to be passed to the
    optL1 or optL2 function of the penalized R package.  See those help
    pages, and it may be desirable to test these arguments directly on
    optL1 or optL2 before using this more CPU-consuming and complex
    function.
  }
}
\details{
  This function calculates cross-validated risk score predictions,
  tuning a penalized regression model using the optL1 or optL2 functions
  of the penalized R package, for each iteration of the
  cross-validation.  Tuning is done by cross-validation in the training
  samples only.  Test samples are scaled using the shift and scale
  factors determined from the training samples.  parameter.  If
  nprocessors > 1, it uses the SNOW package for parallelization,
  dividing the iterations of the outer cross-validation among the
  specified number of processors.

  Some arguments MUST be passed (through the ... arguments) but which
  are documented for the functions in which they are used.  These
  include, from the opt.splitval function:
  
  optFUN="opt1D" for Lasso or Ridge regression, or "opt2D" for Elastic
  Net.  See the help pages for opt1D and opt2D for additional
  arguments associated with these functions.
    
  scaling=TRUE to scale each feature (column) of the training sample to
  z-scores.  These same scaling and shifting factors are applied to the
  test data.  If FALSE, no scaling is done.  Note that only data in the
  penalized argument are scaled, not the optional unpenalized argument
  (see documentation for opt1D, opt2D, or cvl from the penalized package
  for descriptions of the penalized and unpenalized arguments).

  nsim=50 this number specifies the number of times to repeat tuning of
  the penalty parameters on different data foldings for the
  cross-validation.

  setpen="L1" or "L2" : if optFUN="opt1D", this sets regression type to
  LASSO or Ridge, respectively.  See ?opt1D.

  L1range, L2range, dofirst, L1gridsize, L2gridsize: options for Elastic
  Net regression if optFUN="opt2D".  See ?opt2D.
  }
\value{
  Returns a vector of cross-validated continuous risk score predictions.
}
\references{
  Waldron L., Pintilie M., Tsao M.-S., Shepherd F. A., Huttenhower C.*, and Jurisica I.*   Optimized
  application of penalized regression methods to diverse genomic
  data. (2010). Under review.  (*equal contribution)
}
\author{
  Levi Waldron et al.
}
\note{
  Depends on the R packages: penalized, snow, rlecuyer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  opt.splitval
}
\examples{
data(beer.exprs)
data(beer.survival)

gene.quant <- apply(beer.exprs,1,quantile,probs=0.75)
dat.filt <- beer.exprs[gene.quant>log2(100),]
gene.iqr <- apply(dat.filt,1,IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr>0.5,])
dat.filt <- t(dat.filt)

library(survival)
surv.obj <- Surv(beer.survival$os,beer.survival$status)

##Cross-validated calculation and assessment of risk scores
##Ideally nsim would be 50, and outerfold and fold would be
##10, but the values below speed computation 200x compared to these
##recommended values.
set.seed(1)
preds <- opt.nested.crossval(outerfold=5,nprocessors=1,
             optFUN="opt1D",scaling=TRUE,
             setpen="L1",nsim=1,
             response=surv.obj,penalized=dat.filt,fold=5,positive=FALSE,standardize=FALSE,trace=FALSE)

preds.dichot <- preds > median(preds)

coxfit.continuous <- coxph(surv.obj~preds)
coxfit.dichot <- coxph(surv.obj~preds.dichot)
summary(coxfit.continuous)
summary(coxfit.dichot)

if(require(survivalROC)){
  nobs <- length(preds)
  cutoff <- 12
  preds.roc <- survivalROC(Stime=beer.survival$os,status=beer.survival$status,
                         marker=preds,predict.time=cutoff,span = 0.25*nobs^(-0.20))
  plot(preds.roc$FP, preds.roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
       xlab=paste( "FP", "\n", "AUC = ",round(preds.roc$AUC,3)), 
       ylab="TP",main="LASSO predictions\n ROC curve at 12 months")
  abline(0,1)
}
##clean-up
rm(nobs,cutoff,preds.roc,coxfit.dichot,coxfit.continuous,preds,preds.dichot,dat.filt,gene.iqr,gene.quant)
}
\keyword{ regression }
\keyword{ survival }
