\name{PSM.estimate}
\alias{PSM.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate population parameters}
\description{
  Estimates population parameters in a mixed effects model based on
  stochastic differential equations by use of maximum likelihood and the
  Kalman filter.
}
\usage{
PSM.estimate(Model, Data, Par, CI = FALSE, trace = 0, optimizer = "optim",controllist=NULL,fast=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
    A list containing the following elements:

    \describe{
      \item{\code{Matrices = function(phi)}}{
	Defines the matrices \eqn{A}, \eqn{B}, \eqn{C} and \eqn{D} in
	the model equation. Must return a list of matrices named
	\code{matA}, \code{matB}, \code{matC} and \code{matD}.
	If there is no input, \code{matB} and \code{matD} may be omitted
	by setting them to \code{NULL}.
      }
      \item{\code{X0 = function(Time, phi, U)}}{
	Defines the model state at \code{Time[1]} before
	update. \code{Time[1]} and \code{U[,1]} can be used in the
	evaluation of X0. Must return a column matrix.
      }
      \item{\code{SIG = function(phi)}}{
	Defines the matrix \eqn{\sigma}{SIG} for the diffusion
	term. Returns a square matrix. 
      }
      \item{\code{S = function(phi)}}{
	Defines a covariance matrix for the observation noise. Returns a
	square matrix. 
      }
      \item{\code{h = function(eta,theta,covar)}}{
	Second stage model. Defines how random effects (\code{eta}) and
	covariates (\code{covar}) affects the fixed effects parameters (\code{theta}).
	Must return a list \code{phi} of individual parameters. 
      }
      \item{\code{ModelPar = function(THETA)}}{
	Defines the population parameters to be optimized. Returns a
	list containing 2 elements, named: 
	\describe{
	  \item{\code{theta}}{A list of fixed effects parameters
	    \eqn{\theta} which are used as input to the function
	    \code{h} listed above.} 
	  \item{\code{OMEGA}}{A square covariance matrix
	    \eqn{\Omega}{OMEGA} for the random effects. If
	    \code{OMEGA} is missing or \code{NULL} then no 2nd stage
	    model is used and \code{phi} = \code{theta} }
	}
      }
      \item{\code{Dose}}{
	A list containing the 3 elements listed below. If the element
	\code{Dose} is missing or \code{NULL}, no dose is assumed. 
	\describe{
	  \item{\code{Time}}{ A vector of timepoints for the dosing. Each must
	    coinside with a measurement time. Remember to insert a
	    missing measurement in \code{Data$Y} if a corresponding
	    timepoint is not present. Dose is considered added to the system
	    just after the measurement.}
	  \item{\code{State}}{ A vector with indexes of the state for dosing.}
	  \item{\code{Amount}}{ A vector of amounts to be added.}
	}
      }
    }
  }
  \item{Data}{An unnamed list where each element contains
    data for one individual. Each element in \code{Data} is a list
    containing: 

    \describe{
      \item{\code{Time}}{A vector of timepoints for measurements}

      \item{\code{Y}}{A matrix of multivariate observations for each timepoint, where
	each column is a multivariate measurement. Y may contain \code{NA} for missing
	observations and a column may consist of both some or only
	\code{NA}s. The latter is useful if a dose is given when no measurement
	is taken. }

      \item{\code{U}}{A matrix of multivariate input to the model for each
	timepoint. \code{U} is assumed constant between measurements and may not
	contain any \code{NA}. If \code{U} is ommitted, the model is 
	assumed to have no input and \code{matB} and \code{matD} need no
	to be specified.}
    }
  }
  \item{Par}{A list containing the following elements:

    \describe{
      \item{\code{Init}}{A vector with initial estimates for \code{THETA}, vector of
	population parameters to be optimized.}

      \item{\code{LB}, \code{UB}}{: Two vectors with lower and upper bounds for
	parameters. If ommitted, the program performs unconstrained
	optimization. It is highly recommended to specify bounds to ensure
	robust optimization. }
    }
  }
  \item{CI}{ Boolean. If true, the program estimates 95\% confidence
    intervals based on the Hessian of the likelihood function. The
    Hessian is estimated by the optimizer.
  }
  \item{trace}{ Non-negative integer. If positive, tracing
    information on the progress of the optimization is produced. Higher
    values produces more tracing information.
  }
  \item{optimizer}{ May be either \code{"optim"} or
    \code{"nlm"}. Chooses between the native optimizers in R.}

  \item{controllist}{ Controllist for the optimizer. Default in package
    is:
    \describe{
      \item{\code{maxit}}{100}
      \item{\code{abstol}}{1e-5}
      \item{\code{trace}}{Set to trace as in argument}
      \item{\code{REPORT}}{1}
      }
  }
  \item{fast}{  Boolean. Use compiled Fortran code for faster
    estimation.}

  \item{...}{Further arguments to be passed to the chosen optimzer.}
}

\details{
  The first stage model describing intra-individual variations is
  defined as: 

  \deqn{ dX = (AX + BU)dt + \sigma d\omega}{ dX = (A*X + B*U)dt + SIG*dw}
  \deqn{ Y = CX + DU + e}{Y = C*X + D*U + e}

  where \eqn{e \sim N(0,S)}{e ~ N(0,S)} and \eqn{\omega}{w} is a
  standard Brownian motion. The matrices A, B, C, D, \eqn{\sigma}{SIG}
  and S may depend on the individual parameters \eqn{\phi_i}{phi}. 

  The second stage model describing inter-individual variations is
  defined as:

  \deqn{\phi_i = h(\eta_i,\theta,Z_i)}{phi = h(eta,theta,Z)}

  where \eqn{\eta_i \sim N(0,\Omega)}{eta ~ N(0,OMEGA)}, \eqn{\theta}
  are the fixed effect parameters and \eqn{Z_i}{Z} are covariates for individual i.
  
}
\value{
  A list containing the following elements:
  
  \item{NegLogL}{Value of the negative log-likelihood function at optimum.}
  \item{THETA}{Population parameters at optimum}
  \item{CI}{95\% confidence interval for the estimated parameters}
  \item{opt}{Raw output from the optimizer}
  \item{sec}{Time for the estimation in seconds}
  
}

\note{
  For further details please also read the package vignette pdf-document
  by writing vignette("PSM",package="PSM") in R.
}

\references{ Please refer to the help page for \code{\link{PSM}}. }

\author{ Stig B. Mortensen and Søren Klim}

\seealso{
  \code{\link{PSM}}, \code{\link{PSM.smooth}},
  \code{\link{PSM.simulate}}
} 

\examples{
cat("\nExamples are included in the package vignette.\n")
vignette("PSM",package="PSM")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series
