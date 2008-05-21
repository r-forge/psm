\name{APL.KF.individualloop}
\alias{APL.KF.individualloop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Inner loop of APL.KF }
\description{
  This in an internal function. Contains the innner loop over
  individuals for APL.KF.
}
\usage{
APL.KF.individualloop(theta, OMEGA, Model, Data, GUIFlag = 0, fast=TRUE,Linear)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ Fixed effects parameters  }
  \item{OMEGA}{ Covariance matrix for random effects}
  \item{Model}{ Model list}
  \item{Data}{ Data list}
  \item{GUIFlag}{ Display tracing information }
  \item{fast}{ Boolean. Use compiled Fortran code }
  \item{Linear}{Boolean. Is the model Linear or Non-Linear}
}
\value{
  Returns list containing:
  \item{LiPart_i}{Value of individual neg. log-likelihood.}
  \item{eta_i}{Random effect at optimum of
    indiv. neg. log-likelihood.}
  \item{optimStat_i}{Vector of selected output from optim containg the
      elements convergence, value and counts[1].}
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
