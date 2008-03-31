\name{APL.KF}
\alias{APL.KF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{<INTERNAL> Evaluate population likelihood function}
\description{
  This is an internal function. It evaluates the population likelihood function.
}
\usage{
APL.KF(THETA, Model, Pop.Data, LB = NULL, UB = NULL, GUIFlag = 0, longOutput = FALSE, fast=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{THETA}{ Parameters to be optimized}
  \item{Model}{ Model list }
  \item{Pop.Data}{ Data list }
  \item{LB}{ Lower bounds }
  \item{UB}{ Upper bounds }
  \item{GUIFlag}{ Display tracing information}
  \item{longOutput}{ Output negLogLike, etaList, optimStat}
  \item{fast}{ Boolean. Use compiled Fortran code }
}
\value{
  Returns value of population likelihood functiion.
}
\author{ Stig B. Mortensen and Søren Klim}
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

