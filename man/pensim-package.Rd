\name{pensim-package}
\alias{pensim-package}
\alias{pensim}
\docType{package}
\title{
Functions and data for simulation of high-dimensional data and parallelized repeated penalized regression
}
\description{
Simulation of continuous, correlated high-dimensional data with
time-to-event or binary response, and parallelized functions Lasso,
Ridge, and Elastic Net penalized regression.
}
\details{
\tabular{ll}{
Package: \tab pensim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-01-20\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
This package contains a function for simulation of correlated
high-dimensional data with binary or time-to-event response, and
functions for parallelized repeated optimizations of Lasso, Ridge, and
Elastic Net penalized regression penalties.
}
\author{
Levi Waldron

Maintainer: Levi Waldron <lwaldron.research@gmail.com>
}
\references{
  Waldron L., Pintilie M., Tsao M.-S., Shepherd F. A., Huttenhower C.*, and Jurisica I.*   Optimized
  application of penalized regression methods to diverse genomic
  data. (2010). Under review.  (*equal contribution)
}
\keyword{ package }
\keyword{ survival }
\keyword{ datagen }
\keyword{ regression }
\keyword{ multivariate }

\seealso{
%%\code{\link[snow:makeCluster]{snow package (simple network of workstations)}}
}
\examples{
set.seed(9)
##create some data, with one of a group of five correlated variables
##having an association with the binary outcome:
x <- create.data(nvars=c(15,5),cors=c(0,0.8),associations=c(0,2),firstonly=c(TRUE,TRUE),nsamples=50,response="binary",logisticintercept=0.5)
x$summary
##predictor data frame and binary response vector
pen.data <- x$data[,-match("outcome",colnames(x$data))]
response <- x$data[,match("outcome",colnames(x$data))]
##lasso regression
output <- opt1D(nsim=10,nprocessors=2,penalized=pen.data,response=response)
cc <- output[which.max(output[,"cvl"]),-1:-3]  ##non-zero b.* are true positives
}
