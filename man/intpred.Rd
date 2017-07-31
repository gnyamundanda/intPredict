\name{intpred}
\alias{intpred}
\title{
intpred
}
\description{
A pipeline for different gene selection and class prediction methods. 
}
\usage{
intpred(edata, labs, seqp, fac, seed = 166)
}
\arguments{
  \item{edata}{
Expression data, samples by genes.
}
  \item{labs}{
Sample known classes.
}
  \item{seqp}{
Sequence of gene selection.
}
  \item{fac}{
Ratio of splitting learning and testing data.
}
  \item{seed}{
Fix randomness.
}
}
\details{
A pipeline for different gene selection and class prediction methods. 
}
\value{
Misclassification error rate for each method and the selected genes are in your working directory. 
}
\examples{
#require(devtools)
#install_github("syspremed/intPredict")
data(edata) 
data(labs)
seqp=c(seq(2,44,2),47)
fac=2/3
intpred(edata,labs,seqp,fac,seed=166)
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
