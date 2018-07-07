\name{Study_FWER}
\alias{Study_FWER}
\title{Study the FWER of models P1 to N3}
\usage{
Study_FWER(model, sizes, nsim, nsim2)
}
\arguments{
  \item{model}{Name of model, can be "P1", "P2", "P3", "N1", "N2", "N3"}
  \item{sizes}{List with sizes that needs to be studied}
  \item{nsim}{Number of correlation matrices that is drawn per size}
  \item{nsim2}{Number of simulations per correlation matrix}
}
\value{Array with nsim * sim2 FWERs and the minimum, maximum and mean correlation}
\description{
Used in the paper
}
\examples{
Study_FWER("P3", c(2, 3, 5, 10), 10, 100)
}

