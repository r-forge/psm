\name{APL.KF.gr}
\alias{APL.KF.gr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Gradient of APL.KF }
\description{
  This is an internal function. Evaluates the gradient of APL.KF.
}
\usage{
APL.KF.gr(THETA, Model, Pop.Data, LB = NULL, UB = NULL, GradSTEP = 1e-04, GUIFlag = 0, fast=TRUE,Linear=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{THETA}{ Parameters to be optimized}
  \item{Model}{ Model list }
  \item{Pop.Data}{ Data list }
  \item{LB}{ Lower bounds }
  \item{UB}{ Upper bounds }
  \item{GradSTEP}{Steplength for gradient}
  \item{GUIFlag}{ Display tracing information}
  \item{fast}{ Boolean. Use compiled Fortran code }
  \item{Linear}{Boolean. Is the model Linear or Non-Linear}
}
\value{
  Vector with gradient of APL.KF.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
