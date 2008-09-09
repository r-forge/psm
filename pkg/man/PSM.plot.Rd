\name{PSM.plot}
\alias{PSM.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Basic plots of output}
\description{
  Create basic plots for objects from PSM.smooth. 
}

\usage{
  PSM.plot(Smooth, Data, indiv = NULL, type = c('Xs','YsObs')) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Smooth}{
    Object from \code{\link{PSM.smooth}}.
  }
  \item{Data}{
    Data list.*
  }
  \item{indiv}{
    A vector of integers with which individuals to include.
  }
  \item{type}{
    A vector of strings listing the types of plots to create. The
    possibilities are:     Xp, Xf, Xs, Yp, Ys, YpObs, YsObs, res.p, acf.p, res.s, acf.s
  }
  * See description in \link{PSM.estimate}.
}

\value{
  None (invisible \code{NULL}).
}

\note{
  For further details please also read the package vignette pdf-document
  by writing vignette("PSM") in R.
}

\references{Please visit \url{http://www.imm.dtu.dk/psm} or refer to the
  help page for \code{\link{PSM}}. } 

\author{ Stig B. Mortensen and Søren Klim }

\seealso{   \code{\link{PSM}}, \code{\link{PSM.estimate}}, \code{\link{PSM.smooth}},
  \code{\link{PSM.simulate}}, \code{\link{PSM.template}}
}
\examples{
cat("\nExamples are included in the package vignette.\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
