\name{StudyMiniMax}
\alias{StudyMiniMax}
\title{Compute the minimax value of lambda within a range}
\usage{
StudyMiniMax(npi,minpi,maxpi,ndelta,maxdelta,nepsilon,maxepsilon,nlambda,minlambda)
}
\arguments{
  \item{npi}{number of pi0-values to be studied}
#  minpi = 0.5
#minimum value of pi0 (this value is included in the studied interval)
#  maxpi = 1.0
#maximum value of pi0 (this value is NOT included in the studied interval)
#  ndelta = 100
#number of delta-values to be studied
#  maxdelta = 10
#maximum value of delta (this value is included in the studied interval)
#  nepsilon = 100
#number of epsilon-values to be studied
#  maxepsilon = 10
#maximum value of epsilon (this value is included in the studied interval)
#  nlambda = 500
#number of lambda-values to be studied
#  minlambda = 0.5
#minimum value of lambda (this value is included in the studied interval)
  \item{minpi}{Minimum value of pi0 (this value is included in the studied interval)}
  \item{maxpi}{Maximum value of pi0 (this value is NOT included in the studied interval)}
  \item{ndelta}{Number of delta-values to be studied}
  \item{maxdelta}{Maximum value of delta (this value is included in the studied interval)}
  \item{nepsilon}{Number of epsilon-values to be studied}
  \item{maxepsilon}{Maximum value of epsilon (this value is included in the studied interval)}
  \item{nlambda}{Number of lambda-values to be studied}
  \item{minlambda}{Minimum value of lambda (this value is included in the studied interval) }
}
\value{The minimax solution of lambda, the minimax risk, pi with minimax risk, delta with minimax risk, epsilon with minimax risk}
\description{
Compute the minimax value of lambda within a range. The grid is defined in the arguments.
}
\examples{
StudyMiniMax(10, 0.5, 1.0, 10, 10, 10, 10, 50, 0.5)
}
