\name{PSM-package}
\alias{PSM}
\docType{package}
\title{
  Population Stochastic Modelling
}
\description{
  Mixed-effects models using Stochastic Differential Equations
  
  This package provides functions for estimation and simulation of
  multivariate linear mixed-effects models using stochastic
  differential equations. The setup allows for multidimensional input,
  specification of dosage regimen and is able to return smoothed
  estimates of model states.
}

\details{
  \tabular{ll}{
    Package: \tab PSM\cr
    Name: \tab Population Stochastic Modelling\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2007-11-23\cr
    License: \tab GPL (>=2)\cr
  }
  Function overview:
  
  \code{PSM.estimate}\cr
  Estimate population parameters for any linear model.
  
  \code{PSM.smooth}\cr
  Optimal estimates of model states based on estimated parameters.
  
  \code{PSM.simulate}\cr
  Simulate data for multiple individuals.

  This package depends on the mexp-package by David Firth for
  calculation of the matrix exponential. It is currently not available
  at CRAN but can be download from
  \url{http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic/firth/software/mexp/} 
}
\author{
  Stig B. Mortensen and Søren Klim
  
  Maintainer: Søren Klim <skl@imm.dtu.dk>
}
\references{
  Stig B. Mortensen, Søren Klim, Bernd Dammann, Niels R. Kristensen,
  Henrik Madsen, Rune V. Overgaard. A matlab framework for estimation of
  NLME models using stochastic differential equations: Application for
  estimation of insulin secretion rates. J Pharmacokinet Pharmacodyn
  (2007) 34:623-642.
}
\keyword{ package }
\seealso{
  \code{\link{PSM.estimate}}, \code{\link{PSM.smooth}}, \code{\link{PSM.simulate}} 
}
%\examples{
%print("There are no examples.")
%}
