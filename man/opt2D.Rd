\name{opt2D}
\alias{opt2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallelized, two-dimensional tuning of Elastic Net L1/L2 penalties
}
\description{
  This function implements parallelized two-dimensional optimization of Elastic Net
  penalty parameters.  This is accomplished by scanning a regular grid
  of L1/L2 penalties, then using the top five CVL penalty combinations
  from this grid as starting points for the convex optimization problem.
}
\usage{
opt2D(nsim,
      L1range = c(0.001, 100),
      L2range = c(0.001, 100),
      dofirst = "both",
      nprocessors = 1,
      L1gridsize = 10, L2gridsize = 10,
      cl = NULL,
      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{
    Number of times to repeat the simulation (around 50 is suggested)
}
  \item{L1range}{
    numeric vector of length two, giving minimum and maximum constraints
    on the L1 penalty
}
  \item{L2range}{
    numeric vector of length two, giving minimum and maximum constraints
    on the L2 penalty
}
  \item{dofirst}{
    "L1" to optimize L1 followed by L2, "L2" to optimize L2 followed by
    L1, or "both" to optimize both simultaneously in a two-dimensional optimization.
}
  \item{nprocessors}{
        An integer number of processors to use.
}
  \item{L1gridsize}{
    Number of values of the L1 penalty in the regular grid of L1/L2 penalties
}
  \item{L2gridsize}{
    Number of values of the L2 penalty in the regular grid of L1/L2 penalties
}
  \item{cl}{
    Optional cluster object created with the makeCluster() function of
    the parallel package.  If this is not set, pensim calls
    makeCluster(nprocessors, type="SOCK").   Setting this parameter
    can enable parallelization in more diverse scenarios than multi-core
    desktops; see the documentation for the parallel package.  Note that if
    cl is user-defined, this function will not automatically run
    parallel::stopCluster() to shut down the cluster.
}
  \item{\dots}{
    arguments passed on to optL1 and optL2 (dofirst="L1" or "L2"), or
    cvl (dofirst="both") functions of the penalized R package
}
}
\details{
  This function sets up a SNOW (Simple Network of Workstations) "sock"
  cluster to parallelize the task of repeated tunings the Elastic Net
  penalty parameters.  Three methods are implemented, as described by
  Waldron et al. (2011): lambda1 followed by lambda2 (lambda1-lambda2),
  lambda2 followed by lambda1 (lambda2-lambda1), and lambda1 with
  lambda2 simultaneously (lambda1+lambda2).  Tuning of the penalty
  parameters is done by the optL1 or optL2 functions of the penalized R
  package.
}
\value{
  Returns a matrix with the following columns:
  \item{L1}{optimized value of the L1 penalty parameter}
  \item{L2}{optimized value of the L2 penalty parameter}
  \item{cvl}{optimized cross-validated likelihood}
  \item{convergence}{0 if the optimization converged, non-zero otherwise
    (see stats:optim for details)}
  \item{fncalls}{number of calls to cvl function during optimization}
  \item{coef_1, coef_2, ..., coef_n}{argmax coefficients for the model
    with this value of the tuning parameter}

  The matrix contains one row for each repeat of the regression.
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

\seealso{
optL1, optL2, cvl
}
\examples{
data(beer.exprs)
data(beer.survival)

## Select just 100 genes to speed computation:
set.seed(1)
beer.exprs.sample <- beer.exprs[sample(1:nrow(beer.exprs), 100),]

## Apply an unreasonably strict gene filter here to speed computation
## time for the Elastic Net example.
gene.quant <- apply(beer.exprs.sample, 1, quantile, probs = 0.75)
dat.filt <- beer.exprs.sample[gene.quant > log2(150),]
gene.iqr <- apply(dat.filt, 1, IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr > 1,])
dat.filt <- t(dat.filt)

## Define training and test sets
set.seed(9)
trainingset <- sample(rownames(dat.filt), round(nrow(dat.filt) / 2))
testset <-
  rownames(dat.filt)[!rownames(dat.filt) \%in\% trainingset]

dat.training <- data.frame(dat.filt[trainingset,])
pheno.training <- beer.survival[trainingset,]

library(survival)
surv.training <- Surv(pheno.training$os, pheno.training$status)

dat.test <- data.frame(dat.filt[testset,])
all.equal(colnames(dat.training), colnames(dat.test))
pheno.test <- beer.survival[testset,]
surv.test <- Surv(pheno.test$os, pheno.test$status)

set.seed(1)
##ideally set nsim=50, fold=10, but this takes 100x longer.
system.time(
  output <- opt2D(
    nsim = 1,
    L1range = c(0.1, 1),
    L2range = c(20, 1000),
    dofirst = "both",
    nprocessors = 1,
    response = surv.training,
    penalized = dat.training,
    fold = 5,
    positive = FALSE,
    standardize = TRUE
  )
)

cc <- output[which.max(output[, "cvl"]),-1:-5]
output[which.max(output[, "cvl"]), 1:5]  #small L1, large L2
sum(abs(cc) > 0)  #number of non-zero coefficients

preds.training <- as.matrix(dat.training) \%*\% cc
preds.training.median <- median(preds.training)
preds.training.dichot <-
  ifelse(preds.training > preds.training.median, "high risk", "low risk")
preds.training.dichot <-
  factor(preds.training.dichot[, 1], levels = c("low risk", "high risk"))
preds.test <- as.matrix(dat.test) \%*\% cc
preds.test.dichot <-
  ifelse(preds.test > preds.training.median, "high risk", "low risk")
preds.test.dichot <-
  factor(preds.test.dichot[, 1], levels = c("low risk", "high risk"))

coxphfit.training <- coxph(surv.training ~ preds.training.dichot)
survfit.training <- survfit(surv.training ~ preds.training.dichot)
summary(coxphfit.training)
coxphfit.test <- coxph(surv.test ~ preds.test.dichot)
survfit.test <- survfit(surv.test ~ preds.test.dichot)
summary(coxphfit.test)

(p.training <-
    signif(summary(coxphfit.training)$logtest[3], 2))  #likelihood ratio test
(hr.training <- signif(summary(coxphfit.training)$conf.int[1], 2))
(hr.lower.training <- summary(coxphfit.training)$conf.int[3])
(hr.upper.training <- summary(coxphfit.training)$conf.int[4])
par(mfrow = c(1, 2))
plot(
  survfit.training,
  col = c("black", "red"),
  conf.int = FALSE,
  xlab = "Months",
  main = "TRAINING",
  ylab = "Overall survival"
)
xmax <- par("usr")[2] - 50
text(
  x = xmax,
  y = 0.4,
  lab = paste("HR=", hr.training),
  pos = 2
)
text(
  x = xmax,
  y = 0.3,
  lab = paste("p=", p.training, "", sep = ""),
  pos = 2
)
tmp <- summary(preds.training.dichot)
text(
  x = xmax,
  y = c(0.2, 0.1),
  lab = paste(tmp, names(tmp)),
  col = 1:2,
  pos = 2
)
## Now the test set.
## in the test set,  HR=1.7 is not significant - not surprising with the
## overly strict non-specific pre-filter (IQR>1,  75th percentile > log2(150)
(p.test <-
    signif(summary(coxphfit.test)$logtest[3], 2))  #likelihood ratio test
(hr.test <- signif(summary(coxphfit.test)$conf.int[1], 2))
(hr.lower.test <- summary(coxphfit.test)$conf.int[3])
(hr.upper.test <- summary(coxphfit.test)$conf.int[4])
plot(
  survfit.test,
  col = c("black",  "red"),
  conf.int = FALSE,
  xlab = "Months",
  main = "TEST"
)
text(
  x = xmax,
  y = 0.4,
  lab = paste("HR=", hr.test),
  pos = 2
)
text(
  x = xmax,
  y = 0.3,
  lab = paste("p=", p.test, "", sep = ""),
  pos = 2
)
tmp <- summary(preds.test.dichot)
text(
  x = xmax,
  y = c(0.2, 0.1),
  lab = paste(tmp, names(tmp)),
  col = 1:2,
  pos = 2
)
}

\keyword{ regression }
\keyword{ survival }
