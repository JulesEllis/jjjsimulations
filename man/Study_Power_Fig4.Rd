\name{Study_Power_Fig4}
\alias{Study_Power_Fig4}
\title{Study the power in a range of situations, as displayed in figure 4 of the Supplementary Material}
\usage{
Study_Power_Fig4(nsim,effects, sizes, lambda, alpha, fraction)
}
\arguments{
  \item{nsim}{Number of simulations}
  \item{effects}{List with noncentrality parameters}
  \item{sizes}{List with numbers of hypotheses}
  \item{lambda}{Conditionalization threshold, between 0 and 1}
  \item{alpha}{Nominal significance level}
  \item{fraction}{Fraction of true hypotheses}
}
\value{Array with values of minimum power for Bonferroni, FGS, Conditionalized Bonferroni, and Conditionalized FGS}
\description{
Computes the minimum power of various multiple testing procedures, as defined in Figure 4 of the Supplementary Material.
}
\examples{
Study_Power_Fig4(100, c(0.5, 1.0, 1.5), c(200, 500), 0.5, 0.05, 0.75)
}
