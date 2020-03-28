\name{pensim-package}
\alias{pensim-package}
\alias{pensim}
\docType{package}
\title{
Functions and data for simulation of high-dimensional data and parallelized repeated penalized regression
}
\description{
Simulation of continuous, correlated high-dimensional data with
time-to-event or binary response, and parallelized functions for Lasso,
Ridge, and Elastic Net penalized regression model training and
validation by split-sample or nested cross-validation.  See the help
page for opt.nested.crossval() for the most extensive usage examples.
}
\details{
\tabular{ll}{
Package: \tab pensim\cr
Type: \tab Package\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
Model training and validation by Lasso, Ridge, and Elastic Net penalized
regression.   This package also contains a function for simulation of correlated
high-dimensional data with binary or time-to-event response.
}
\author{
Levi Waldron

Maintainer: Levi Waldron <lwaldron.research@gmail.com>
}
\references{
  Waldron L, Pintilie M, Tsao M-S, Shepherd FA, Huttenhower C*, Jurisica
  I*: Optimized application of penalized regression methods to diverse
  genomic data. Bioinformatics 2011, 27:3399-3406.  (*equal contribution)
}
\keyword{ package }
\keyword{ survival }
\keyword{ datagen }
\keyword{ regression }
\keyword{ multivariate }

\seealso{
penalized-package
}
\examples{
set.seed(9)
## 
## create some data,  with one of a group of five correlated variables
## having an association with the binary outcome:
## 
x <- create.data(
  nvars = c(10, 3),
  cors = c(0, 0.8),
  associations = c(0, 2),
  firstonly = c(TRUE, TRUE),
  nsamples = 50,
  response = "binary",
  logisticintercept = 0.5
)
x$summary
##
##predictor data frame and binary response vector
##
pen.data <- x$data[, -match("outcome", colnames(x$data))]
response <- x$data[, match("outcome", colnames(x$data))]
## lasso regression.  Note that epsilon=1e-2 is passed onto optL1,  and
## reduces the precision of the tuning compared to the default 1e-10.
output <-
  opt1D(
    nsim = 1,
    nprocessors = 1,
    penalized = pen.data,
    response = response,
    epsilon = 1e-2
  )
cc <-
  output[which.max(output[, "cvl"]), -1:-3]  ##non-zero b.* are true positives
}
