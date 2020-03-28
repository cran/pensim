\name{opt1D}
\alias{opt1D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parallelized repeated tuning of Lasso or Ridge penalty parameter
}
\description{
This function is a wrapper to the optL1 and optL2 functions of the
penalized R package, useful for parallelized repeated tuning of the
penalty parameters.
}
\usage{
opt1D(nsim = 50, nprocessors = 1, setpen = "L1", cl = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{
    Number of times to repeat the simulation (around 50 is suggested)
}
  \item{nprocessors}{
    An integer number of processors to use.
}
  \item{setpen}{
    Either "L1" (Lasso) or "L2" (Ridge) penalty
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
    arguments passed on to optL1 or optL2 function of the penalized R package
}
}
\details{
  This function sets up a SNOW (Simple Network of Workstations) "sock"
  cluster to parallelize the task of repeated tunings the L1 or L2 penalty
  parameter.  Tuning of the penalty parameters is done by the optL1 or
  optL2 functions of the penalized R package.
}
\value{
  Returns a matrix with the following columns:
  \item{L1 (or L2)}{optimized value of the penalty parameter}
  \item{cvl}{optimized cross-validated likelihood}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  optL1, optL2
}
\examples{
data(beer.exprs)
data(beer.survival)

##select just 100 genes to speed computation:
set.seed(1)
beer.exprs.sample <- beer.exprs[sample(1:nrow(beer.exprs), 100),]

gene.quant <- apply(beer.exprs.sample, 1, quantile, probs = 0.75)
dat.filt <- beer.exprs.sample[gene.quant > log2(100),]
gene.iqr <- apply(dat.filt, 1, IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr > 0.5,])
dat.filt <- t(dat.filt)

##define training and test sets
set.seed(1)
trainingset <- sample(rownames(dat.filt), round(nrow(dat.filt) / 2))
testset <-
  rownames(dat.filt)[!rownames(dat.filt) \%in\% trainingset]

dat.training <- data.frame(dat.filt[trainingset, ])
pheno.training <- beer.survival[trainingset, ]

library(survival)
surv.training <- Surv(pheno.training$os, pheno.training$status)

dat.test <- data.frame(dat.filt[testset, ])
all.equal(colnames(dat.training), colnames(dat.test))
pheno.test <- beer.survival[testset, ]
surv.test <- Surv(pheno.test$os, pheno.test$status)

##ideally nsim should be on the order of 50,  but this slows computation
##50x without parallelization.
set.seed(1)
output <-
  pensim::opt1D(
    nsim = 1,
    nprocessors = 1,
    setpen = "L2",
    response = surv.training,
    penalized = dat.training,
    fold = 3,
    positive = FALSE,
    standardize = TRUE,
    minlambda2 = 1,
    maxlambda2 = 100
  )

cc <- output[which.max(output[, "cvl"]),-(1:2)]  #coefficients
sum(abs(cc) > 0)  #count non-zero coefficients

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
  x = c(xmax, xmax),
  y = c(0.2, 0.1),
  lab = paste(tmp, names(tmp)),
  col = 1:2,
  pos = 2
)
(p.test <-
    signif(summary(coxphfit.test)$logtest[3], 2))  #likelihood ratio test
(hr.test <- signif(summary(coxphfit.test)$conf.int[1], 2))
(hr.lower.test <- summary(coxphfit.test)$conf.int[3])
(hr.upper.test <- summary(coxphfit.test)$conf.int[4])
plot(
  survfit.test,
  col = c("black", "red"),
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
  x = c(xmax, xmax),
  y = c(0.2, 0.1),
  lab = paste(tmp, names(tmp)),
  col = 1:2,
  pos = 2
)
}
\keyword{ regression }
\keyword{ survival }
