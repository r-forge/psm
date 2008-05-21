\name{ExtKalmanSmoother}
\alias{ExtKalmanSmoother}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Non-linear Kalman smoother }
\description{
  This in an internal function. Performs a non-linear Kalman smoothing.
}
\usage{
ExtKalmanSmoother(phi, Model, Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{ Individual parameters }
  \item{Model}{ Model list}
  \item{Data}{ Data list}
}

\value{
  Returns a list containing:
  
  TIME, Xs, Ps, negLogLike, Xp, Xf, Yp, KfGain, Pf, Pp, R
}

\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
