\name{scan.l1l2}
\alias{scan.l1l2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function calculate cross-validated likelihood on a regular grid of L1/L2 penalties
}
\description{
This function generates a grid of values of L1/L2 penalties,
then calculated cross-validated likelihood at each point on the grid.
The grid can be regular (linear progression of the penalty values), or
polynomial (finer grid for small penalty values, and coarser grid for
larger penalty values).
}
\usage{
scan.l1l2(L1range = c(0.1, 100.1), L2range = c(0.1, 100.1), L1.ngrid = 50, L2.ngrid = 50, nprocessors = 1, polydegree = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L1range}{
    numeric vector of length two, giving minimum and maximum constraints
    on the L1 penalty
  }
  \item{L2range}{
    numeric vector of length two, giving minimum and maximum constraints
    on the L2 penalty
  }
  \item{L1.ngrid}{
    Number of values of the L1 penalty in the regular grid of L1/L2 penalties
  }
  \item{L2.ngrid}{
    Number of values of the L2 penalty in the regular grid of L1/L2 penalties
  }
  \item{nprocessors}{
    An integer number of processors to use.
  }
  \item{polydegree}{
    power of the polynomial on which the L1/L2 penalty values are fit.
    ie if polydegree=2, penalty values could be y=x^2, x=1,2,3,..., so y=1,4,9,...
  }
  \item{\dots}{
    arguments passed on to cvl function of the penalized R package
  }
}
\details{
  This function sets up a SNOW (Simple Network of Workstations) "sock"
  cluster to parallelize the task of scanning a grid of penalty values
  to search for suitable starting values for two-dimensional
  optimization of the Elastic Net.
}
\value{
  \item{cvl}{matrix of cvl values along the grid}
  \item{L1range}{L1range,L2range}{L2range}
  \item{xlab}{A text string indicating the range of L1 penalties}
  \item{ylab}{A text string giving the range of L2 penalties}
  \item{zlab}{A text string giving the range of cvl values}
  \item{note}{A note to the user that rows of cvl correspond to values of lambda1, columns to lambda2}
}
\references{
  Waldron L., Pintilie M., Huttenhower C.*, and Jurisica I.*   Optimized
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(beer.exprs)
data(beer.survival)

gene.quant <- apply(beer.exprs,1,quantile,probs=0.75)
dat.filt <- beer.exprs[gene.quant>log2(150),]
gene.iqr <- apply(dat.filt,1,IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr>1,])
dat.final <- scale(t(dat.filt))

##define training and test sets
set.seed(9)
trainingset <- sample(rownames(dat.final),round(nrow(dat.final)/2))
testset <- rownames(dat.final)[!rownames(dat.final)\%in\%trainingset]

dat.training <- data.frame(dat.final[trainingset,])
pheno.training <- beer.survival[trainingset,]

library(survival)
surv.training <- Surv(pheno.training$os,pheno.training$status)

dat.test <- data.frame(dat.final[testset,])
all.equal(colnames(dat.training),colnames(dat.test))
pheno.test <- beer.survival[testset,]
surv.test <- Surv(pheno.test$os,pheno.test$status)

set.seed(9)
system.time(output <- scan.l1l2(L1range=c(0.2,10.2),L2range=c(2,50),
L1.ngrid=20,L2.ngrid=50,polydegree=1,nprocessors=2,
response=surv.training, penalized=dat.training, fold=4,positive=FALSE,standardize=FALSE))

##Note that the cvl surface is not smooth because a different folding of
##the data was used for each cvl calculation
image(x=seq(output$L1range[1],output$L1range[2],length.out=nrow(output$cvl)),
y=seq(output$L2range[1],output$L2range[2],length.out=ncol(output$cvl)),
z=output$cvl,
xlab="lambda1",
ylab="lambda2",
main="red is higher cross-validated likelihood")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ survival }
