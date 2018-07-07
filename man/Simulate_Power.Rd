\name{Simulate_Power}
\alias{Simulate_Power}
\title{Study the power of unconditional and conditional tests}
\usage{Simulate_Power(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)}
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
Bonferroni, FGS, Fisher, Tippet, Iplus, LR,
Conditional Bonferroni, Conditional FGS, Conditional Fisher,
Conditional Tippet, Conditional Iplus, Conditional LR}
\description{
Used in the paper
}
\examples{
Simulate_Power(100, 2, 2, 1, 1,15,5,0.5,0.05)
}
