\name{PSM.simulate}
\alias{PSM.simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create simulation data for multiple individuals}
\description{
  Simulates data for multiple individuals in a mixed effects model based
  on stochastic differential equations.
}
\usage{
PSM.simulate(Model, THETA, dt, Ulist = NA, Tlist, individuals = 1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{ ~~Describe \code{Model} here~~ }
  \item{THETA}{ ~~Describe \code{THETA} here~~ }
  \item{dt}{ ~~Describe \code{dt} here~~ }
  \item{Ulist}{ ~~Describe \code{Ulist} here~~ }
  \item{Tlist}{ ~~Describe \code{Tlist} here~~ }
  \item{individuals}{ ~~Describe \code{individuals} here~~ }
}
\details{
  If necessary, more details than the description above
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\references{ Please refer to the help page for \code{\link{PSM}}. }

\author{ Stig B. Mortensen and Søren Klim }

\seealso{  \code{\link{PSM}}, \code{\link{PSM.estimate}},
  \code{\link{PSM.smooth}}
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
