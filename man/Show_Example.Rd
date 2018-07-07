\name{Show_Example}
\alias{Show_Example}
\title{Computes the outcomes of the empirical examples in the paper}
\usage{
Show_Example(example, lambda)
}
\arguments{
  \item{example}{Array with p-values (ravenp1, ravenp2, ravenp3, ravenp4, or natalcarep)}
  \item{lambda}{Conditionalization threshold}
}
\value{The number of p-values less than or equal to lambda; the adjusted critical value; the adjusted p-values}
\description{
Used in Examples 1 and 3 of the paper.
}
\examples{
Show_Example(c(0.01, 0.02, 0.6, 0.7, 0.8, 0.9), 0.5)
}
