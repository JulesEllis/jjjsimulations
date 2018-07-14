\name{Simulate_Power2}
\alias{Simulate_Power2}
\title{Study the power of unconditional and conditional Bonferroni and FGS}
\usage{Simulate_Power2(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)}
\arguments{
  \item{nsim}{Number of simulations}
  \item{ncp0}{Noncentrality parameter of test statistic of true hypotheses, greater or equal to 0}
  \item{ncp1}{Noncentrality parameter of test statistic of false hypotheses, greater or equal to 0}
  \item{sigma0}{Standard error of test statistic of true hypotheses, greater or equal to 0}
  \item{sigma1}{Standard error of test statistic of false hypotheses, greater or equal to 0}
  \item{n0}{Number of true hypotheses}
  \item{n1}{Number of false hypotheses}
  \item{lambda}{Conditionalization threshold, between 0 and 1}
  \item{alpha}{Nominal significance level}
}
\value{List with rejection rates of these multiple testing procedures:
Bonferroni, FGS, Conditional Bonferroni, Conditional FGS}
\description{
Used in the Supplementary Material. Works like Simulate_Power, but evaluaties only the Bonferroni and FGS rules in order to save time.
}
\examples{
Simulate_Power2(100, 2, 2, 1, 1,15,5,0.5,0.05)
}
