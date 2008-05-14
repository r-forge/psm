\name{ExtKalmanFilter}
\alias{ExtKalmanFilter}

\title{ <INTERNAL> Extended Kalmal filter }
\description{ 
  This in an internal function. Performs a Extended Kalman filtering.
}
\usage{
ExtKalmanFilter(phi, Model, Data, outputInternals = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{ Individual parameters }
  \item{Model}{ Model list}
  \item{Data}{ Data list}
  \item{outputInternals}{ If true, output is a list containing
    negLogLike, X, Xf, Yp, KfGain, Pf, Pp, R}
}
\value{
  Returns neg. log-likelihood.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
