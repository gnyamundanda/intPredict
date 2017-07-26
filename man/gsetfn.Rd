\name{gsetfn}
\alias{gsetfn}
\title{
Gene selection function}
\description{
Pipeline for gene selection.}
\usage{
gsetfn(response, array, fac, narray, np, p, nrmeth, ngenes, nsubs, test.idx)
}
\arguments{
  \item{response}{
Sample known classes.
}
  \item{array}{
Expression data, genes by samples
}
  \item{fac}{
Ratio of splitting learning and testing data.
}
  \item{narray}{
Number of samples}
  \item{np}{
Number of gene sets}
  \item{p}{
Size of each gene set}
  \item{nrmeth}{
Number of gene selection methods}
  \item{ngenes}{
Total number of genes}
  \item{nsubs}{
Number of subtypes}
  \item{test.idx}{
Index of samples selected as test}
}
\details{
Pipeline for gene selection}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
