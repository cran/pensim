\name{create.data}
\alias{create.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate correlated predictors with time-to-event or binary outcome
}
\description{
  This function creates multiple groups of predictor variables which may
  be correlated within each group, and binary or survival time (without
  censoring) response according to specified weights of the predictors.
}
\usage{
create.data(nvars = c(100, 100, 100, 100, 600),
            cors = c(0.8, 0, 0.8, 0, 0),
            associations = c(0.5, 0.5, 0.3, 0.3, 0),
            firstonly = c(TRUE, FALSE, TRUE, FALSE, FALSE),
            nsamples = 100,
            censoring = "none",
            labelswapprob = 0,
            response = "timetoevent",
            basehaz = 0.2,
            logisticintercept = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvars}{integer vector giving the number of variables of each
    variable type.  The number of variable types is equal to the length
    of this vector.}
  \item{cors}{integer vector of the same length as nvars, giving the population pairwise Pearson
    correlation within each group.}
  \item{associations}{integer vector of the same length as nvars, giving the associations of each
    type with outcome}
  \item{firstonly}{logical vector of the same length as nvars, specifying whether only the first variable
    of each type is associated with outcome (TRUE) or all variables of that
    type (FALSE)}
  \item{nsamples}{an integer giving the number of observations}
  \item{censoring}{"none" for no censoring, or a vector of length two c(a,b)
    for uniform U(a,b) censoring.}
  \item{labelswapprob}{This provides an option to add uncertainty to
    binary outcomes by randomly switching labels with probability
    labelswapprob.  The probability of a label being swapped is
    independent for each observation.  The value is ignored if
    response is "timetoevent"}
  \item{response}{either "timetoevent" or "binary"}
  \item{basehaz}{baseline hazard, used for "timetoevent"}
  \item{logisticintercept}{intercept which is added to X\%*\%Beta for "binary"}
}
\details{
  This function simulates "predictor" variables in one or more groups, which are standard
  normally distributed.  The user can specify the population correlation
  within each variable group, the association of each variable group to
  outcome, and whether the first or all variables of that type should be
  associated with outcome.  The simulated response variable can be time to event
  with an exponential distribution, or binary survival with a logistic distribution.
}
\value{
  Returns a list with items:
  \item{summary}{a summary of the variable
    types produced}
  \item{associations}{weights of each variable in
    computing the outcome}
  \item{covariance}{covariance matrix used for
  generating potentially correlated random predictors}
  \item{data}{dataframe containing the predictors and response.
  Response is the last column for binary outcome ("outcome"), and the
  last two columns for timetoevent outcome ("time" and "cens")}
}
\references{
Waldron L., Pintilie M., Tsao M.-S., Shepherd F. A., Huttenhower C.*, and Jurisica I.*   Optimized
application of penalized regression methods to diverse genomic
data. (2010). Under review.  (*equal contribution)
}
\author{Levi Waldron et al.}
\note{Depends on the MASS package for correlated random number generation}
\examples{
##binary outcome example
set.seed(9)
x <-
  create.data(
    nvars = c(15, 5),
    cors = c(0, 0.8),
    associations = c(0, 2),
    firstonly = c(TRUE, TRUE),
    nsamples = 50,
    response = "binary",
    logisticintercept = 0.5
  )
summary(x)
x$summary
model <- glm(outcome ~ .,  data = x$data, family = binomial)
summary(model)
dat <- t(as.matrix(x$data[, -match("outcome", colnames(x$data))]))
heatmap(dat, ColSideColors = ifelse(x$data$outcome == 0, "black", "white"))

##censored survival outcome example:
set.seed(1)
x <- create.data(
  nvars = c(15, 5),
  cors = c(0, 0.8),
  associations = c(0, 2),
  firstonly = c(TRUE, TRUE),
  nsamples = 50,
  censoring = c(2, 10),
  response = "timetoevent"
)
sum(x$data$cens == 0) / nrow(x$data)  #34 percent censoring

library(survival)
surv.obj <- Surv(x$data$time, x$data$cens)
plot(survfit(surv.obj ~ 1), ylab = "Survival probability", xlab = "time")
}
\keyword{ datagen }
\keyword{ survival }
