\name{IndividualLL.KF}
\alias{IndividualLL.KF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Indivdual neg. log-likelihood function. }
\description{
  This in an internal function. Evaluates the indivdual neg. log-likelihood function.
}
\usage{
IndividualLL.KF(eta, theta, OMEGA, Model, Data, fast=TRUE,Linear=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ Random effects}
  \item{theta}{ Fixed effects}
  \item{OMEGA}{ Covariance matrix for random effects}
  \item{Model}{ Model list}
  \item{Data}{ Data list}
  \item{fast}{ Boolean. Use compiled Fortran code }
  \item{Linear}{Boolean. Is the model Linear or Non-Linear}
}
\value{
  Returns individual neg. log-likelihood.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
