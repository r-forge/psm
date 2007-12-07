\name{PSM.smooth}
\alias{PSM.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing of model states based on estimated population parameters.}
\description{
  Gives estimates of model states and random effects \eqn{\eta}. The
  function is intended to be used based on population parameters found
  using \code{\link{PSM.estimate}} or to check initial values before
  parameter estimation.
}

\usage{
PSM.smooth(Model, Data, THETA, subsample = 0, trace = 0, etaList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
    Model list.*
  }
  \item{Data}{
    Data list.*
  }
  \item{THETA}{
    Vector of population parameters used for the state estimation.
  }
  \item{subsample}{
    Number of points to estimate states in between
    measurements. The extra points are linearly spaced.
  }
  \item{trace}{
    Non-negative integer. If positive, tracing
    information on the progress of the optimization is produced. Higher
    values produces more tracing information.
  }
  \item{etaList}{
    Matrix where each column contains an etimate of
    \eqn{\eta_i}. \code{etaList} has the same format as the output of 
    \code{\link{PSM.estimate}}. If ommitted, the function will evalutate
    the population likehood function to find estimates of \eqn{\eta_i}{eta}
    for all individuals.
  }
  * See description in \link{PSM.estimate}.
}
\details{
  The function produces three types of estimates.

  \describe{
    \item{Predicted}{Only past measurements are used for the state
      estimate at time t. }
    
    \item{Filtered}{ Only past and the current measurements are used for
      the state estimate at time t.}
    
    \item{Smoothed}{ All measurements (both past and future) are used to
      form the state estimate at time t. This is usually the prefered type
      of state estimate.}
  }
  If \code{subsample}>0 then the data is automatically subsampled to
  provide estimated of the model states between observation time points.
  
}

\value{
  \item{negLogL}{Value of the negative log-likelihood function at \code{THETA}.}
  \item{etaList}{Matrix containing estimated eta's for each individual.}
  \item{smooth}{
    An unnamed list with one element for each individual. Each element
    contains the following elements:
    \describe{
      \item{\code{negLogLike}}{Value of neg. log-likelihood function }
      \item{\code{TIME}}{Possibly subsampled time-vector corresponding
	to the estimated states }
      \item{\code{Xs}, \code{Ps}}{Smoothed state and state co-variance
	estimate} 
      \item{\code{Xf}, \code{Pf}}{Filtered state and state co-variance
	estimate} 
      \item{\code{Xp}, \code{Pp}}{Predicted state and state co-variance
	estimate} 
      \item{\code{Yp}, \code{R}}{Predicted observations and observation
	variances}      
    }
  }
}

\note{
  For further details please also read the package vignette pdf-document
  by writing vignette("PSM",package="PSM") in R.
}

\references{Please refer to the help page for \code{\link{PSM}}. }

\author{ Stig B. Mortensen and Søren Klim }

\seealso{   \code{\link{PSM}}, \code{\link{PSM.estimate}},
  \code{\link{PSM.simulate}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
