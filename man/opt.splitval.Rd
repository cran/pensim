\name{opt.splitval}
\alias{opt.splitval}
\title{Parallelized calculation of split training/test set predictions from L1/L2/Elastic
  Net penalized regression.
}
\description{
  uses a single training/test split to train a penalized
  regression model in the training samples, then use the model to
  calculate values of the linear risk score in the test samples.  This
  function is used by opt.nested.crossval, but can also be used on its own.

  This function support z-score scaling of training data, and application
  of these scaling and shifting coefficients to the test data.  It also
  supports repeated tuning of the penalty parameters and selection of the
  model with greatest cross-validated likelihood.
}
\usage{
opt.splitval(optFUN="opt1D",testset="equal",scaling=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{optFUN}{
    "opt1D" for Lasso or Ridge regression, "opt2D" for Elastic Net.  See
    the help pages for these functions for additional arguments.
  }
  \item{testset}{
    For the opt.splitval function ONLY.
    "equal" for randomly assigned equal training and test sets, or an
    integer vector defining the positions of the test samples in the
    response, penalized, and unpenalized arguments which are passed to
    the optL1, optL2, or cvl functions of the penalized R package.
  }
  \item{scaling}{
    If TRUE, each feature (column) of the training samples (in
    matrix/dataframe specified by the penalized argument) are scaled to
    z-scores, then these scaling and shifting factors are applied to the
    test data.
    If FALSE, no scaling is done.
  }
  \item{\dots}{
    Additional arguments are required, to be passed to the
    optL1 or optL2 function of the penalized R package.  See those help
    pages, and it may be desirable to test these arguments directly on
    optL1 or optL2 before using this more CPU-consuming and complex
    function.
  }
}
\details{
  This function does split sample model training and testing for a
  single split of the data, using the optL1 or optL2 functions
  of the penalized R package, for each
  iteration of the cross-validation.  Scaling of the test samples is
  done independently, using scale factors determined from the training
  samples.  Repeated starts of model training can be parallelized as
  documented in the opt1D and opt2D functions.  This function is used
  for nested cross-validation by the opt.nested.crossval function.
}
\value{
  Returns a vector of cross-validated continuous risk score predictions.
}
\references{
  Waldron L, Pintilie M, Tsao M-S, Shepherd FA, Huttenhower C*, Jurisica
  I*: Optimized application of penalized regression methods to diverse
  genomic data. Bioinformatics 2011, 27:3399-3406.  (*equal contribution)
}
\author{
  Levi Waldron et al.
}
\note{
  Depends on the R packages: penalized, parallel, rlecuyer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  opt1D, opt2D, opt.nested.crossval
}
\examples{
data(beer.exprs)
data(beer.survival)

## select just 250 genes to speed computation:
set.seed(1)
beer.exprs.sample <- beer.exprs[sample(1:nrow(beer.exprs), 250), ]

gene.quant <- apply(beer.exprs.sample, 1, quantile, probs = 0.75)
dat.filt <- beer.exprs.sample[gene.quant > log2(100),]
gene.iqr <- apply(dat.filt, 1, IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr > 0.5,])
dat.filt <- t(dat.filt)

library(survival)
surv.obj <- Surv(beer.survival$os, beer.survival$status)

## Single split training/test evaluation.  Ideally nsim would be 50 and
## fold=10, but this requires 100x more resources.
set.seed(1)
preds50 <- opt.splitval(
  optFUN = "opt1D",
  scaling = TRUE,
  testset = "equal",
  setpen = "L1",
  nsim = 1,
  nprocessors = 1,
  response = surv.obj,
  penalized = dat.filt,
  fold = 5,
  positive = FALSE,
  standardize = FALSE,
  trace = FALSE
)

preds50.dichot <- preds50 > median(preds50)

surv.obj.50 <-
  surv.obj[match(names(preds50), rownames(beer.survival))]
coxfit50.continuous <- coxph(surv.obj.50 ~ preds50)
coxfit50.dichot <- coxph(surv.obj.50 ~ preds50.dichot)
summary(coxfit50.continuous)
summary(coxfit50.dichot)
}
\keyword{ regression }
\keyword{ survival }
