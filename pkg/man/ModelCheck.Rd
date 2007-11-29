\name{ModelCheck}
\alias{ModelCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Check Model and Data structure }
\description{
  The Model used in PSM can easily be misspecified. This function checks
  for dimensionalities and model objects. Furthermore it tests the Model
  objects and the dimensions in the Data set.
}

\usage{
ModelCheck(Model, Data, Par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\code{Model}}{
    Model list*
  }
  \item{\code{Data}}{
    Single subject data - ie. only one element of the \code{Data} for
    \code{PSM.estimate}
  }
  \item{\code{Par}}{
    List of parameters for the Model*.
    \code{Par$Init} is used to test Model functions.
    If present \code{Par$LB} and \code{Par$UB} are tested for validity.
  }

  * See description \link{PSM.estimate}
}

\details{
  
  The ModelCheck checks for:

  \itemize{
    \item Presence of correctly named functions in the Model
    \item Ability to evaluate the functions in the Model
    \item Dimensions of the Model objects vs. Data objects
    \item Positive semidefinit properties of S and SIG
    \item Optionally Parameter bounds
    
  }  
}

\value{
  Boolean indicating whether the Model passed the requirements
}

\references{ Please refer to the help page for \code{\link{PSM}}. }

\author{ Stig B. Mortensen and Søren Klim}

%\note{
%  (Further notes, possibly other sections using section{Warning }{....})
%}

\seealso{
  \code{\link{PSM}}
} 

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
