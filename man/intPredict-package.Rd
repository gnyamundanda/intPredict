\name{intPredict-package}
\alias{intPredict-package}
\alias{intPredict}
\docType{package}
\title{
\packageTitle{intPredict}
}
\description{
\packageDescription{intPredict}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{intPredict}
\packageIndices{intPredict}
A pipeline for different gene selection and class prediction methods.
}
\author{
\packageAuthor{intPredict}

Maintainer: \packageMaintainer{intPredict}
}
\keyword{ package }
\examples{
#require(devtools)
#install_github("syspremed/intPredict")
data(edata) 
data(labs)
seqp=c(seq(2,44,2),47)
fac=2/3
intpred(edata,labs,seqp,fac,seed=166)
}
