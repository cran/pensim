\name{beer.exprs}
\alias{beer.exprs}
\docType{data}
\title{
  Lung adenocarcinoma microarray expression data of Beer et al. (2002)
}
\description{
  Lung adenocarcinomas were profiled by Beer et al. (2002) using Affymetrix hu6800 microarrays.  The data here were normalized from raw .CEL files by RMAExpress (v0.3).  The expression matrix contains expression data for 86 patients with 7,129 probe sets.
}
\usage{
data(beer.exprs)
}
\format{
  A data frame with 7129 probe sets (rows) for 86 patients (columns)
}
\source{
  Beer DG, Kardia SL, Huang C, Giordano TJ, Levin AM, Misek DE, Lin L, Chen G, Gharib TG, Thomas DG, Lizyness ML, Kuick R, Hayasaka S, Taylor JM, Iannettoni MD, Orringer MB, Hanash S: Gene-expression profiles predict survival of patients with lung adenocarcinoma. Nat Med 2002, 8:816-824. 
}
\references{
  Irizarry, R.A., et al. (2003) Summaries of Affymetrix GeneChip probe level data, Nucl. Acids Res., 31, e15+-e15+.
}
\examples{
data(beer.exprs)
mysd <- apply(beer.exprs,1,sd)
%%##100 most changing probesets
beer.subset <- as.matrix(beer.exprs[rank(-mysd) <= 100,])  
heatmap(beer.subset)
}
\keyword{datasets}
