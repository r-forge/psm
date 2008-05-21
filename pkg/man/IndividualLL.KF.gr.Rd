\name{IndividualLL.KF.gr}
\alias{IndividualLL.KF.gr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Gradient of indivdual neg. log-likelihood function. }
\description{
  This in an internal function. Evaluates the gradient of the indivdual
  neg. log-likelihood function. 
}
\usage{
IndividualLL.KF.gr(eta, theta, OMEGA, Model, Data, GradSTEP = 1e-04, GUIFlag = 0, fast=TRUE,Linear=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ Random effects}
  \item{theta}{ Fixed effects}
  \item{OMEGA}{ Covariance matrix for random effects}
  \item{Model}{ Model list}
  \item{Data}{ Data list}
  \item{GradSTEP}{ Step length for gradient }
  \item{GUIFlag}{ Display tracing information }
  \item{fast}{ Boolean. Use compiled Fortran code }
  \item{Linear}{Boolean. Is the model Linear or Non-Linear}
}

\value{
  Returns vector with gradient.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
