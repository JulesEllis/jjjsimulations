\name{CBP}
\alias{CBP}
\title{Applies the conditionalized Bonferroni procedure}
\usage{
CBP(pvalues, lambda)
}
\arguments{
  \item{pvalues}{Array with p-values (e.g., ravenp1, ravenp2, ravenp3, ravenp4, or natalcarep$p)}
  \item{lambda}{Conditionalization threshold}
}
\value{The number of p-values less than or equal to lambda; the adjusted critical value; the adjusted p-values}
\description{
Used in Examples 1, 2 and 3 of the paper.
}
\examples{
CBP(c(0.01, 0.02, 0.6, 0.7, 0.8, 0.9), 0.5)
}
