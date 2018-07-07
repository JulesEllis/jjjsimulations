\name{Relative_Power}
\alias{Relative_Power}
\title{Study the expected ratio of critical values of Bonferroni and conditionalized Bonferroni}
\usage{
Relative_Power(lambda, pi0, delta, epsilon)
}
\arguments{
  \item{lambda}{Conditionalization threshold, between 0 and 1}
  \item{pi0}{Fraction of true hypotheses, between 0 and 1}
  \item{delta}{Noncentrality parameter of true hypotheses, greater or equal to 0}
  \item{epsilon}{Noncentrality parameter of false hypotheses, greater or equal to 0}
}
\value{The expected ratio (rho) of the critical values of the Bonferroni and conditionalized Bonferroni with conditionalization threshold lambda}
\description{
Used in table of paper
}
\examples{
Relative_Power(0.9, 0.8, 2, 2)
}
