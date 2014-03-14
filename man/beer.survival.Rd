\name{beer.survival}
\alias{beer.survival}
\docType{data}
\title{
Survival data for Beer et al. (2002) lung adenocarcinoma study
}
\description{
  Overall survival time for 86 lung adenocarcinoma patients, with 62 of
  the 86 events being censored.  
}
\usage{data(beer.survival)}
\format{
  A data frame with 86 observations on the following 2 variables.
  \describe{
    \item{\code{status}}{a numeric vector}
    \item{\code{os}}{a numeric vector}
  }
}
\source{
Beer DG, Kardia SL, Huang C, Giordano TJ, Levin AM, Misek DE, Lin L,
Chen G, Gharib TG, Thomas DG, Lizyness ML, Kuick R, Hayasaka S, Taylor
JM, Iannettoni MD, Orringer MB, Hanash S: Gene-expression pr
ofiles predict survival of patients with lung adenocarcinoma. Nat Med 2002, 8:816-824.
}
\examples{
data(beer.survival)
library(survival)
surv.obj <- with(beer.survival, Surv(os, status))
surv.obj.rev <- with(beer.survival, Surv(os, 1-status))
survfit(surv.obj.rev~1)  #reverse KM estimate of follow-up time (months)
(my.survfit <- survfit(surv.obj~1))  ##KM estimate of survival
plot(my.survfit, xlab="Time (months)", 
     ylab="KM estimate of overall survival")
legend("bottomright", lty=c(1, 2), pch=-1,
       legend=c("KM estimate", "95 percent confidence interval"))
}
\keyword{datasets}
