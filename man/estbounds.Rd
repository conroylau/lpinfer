% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estbounds.R
\name{estbounds}
\alias{estbounds}
\title{Estimate bounds with shape restrictions}
\usage{
estbounds(
  data = NULL,
  lpmodel,
  kappa = 0,
  norm = 2,
  estimate = TRUE,
  solver = NULL
)
}
\arguments{
\item{data}{An \code{data frame} or a \code{matrix}.}

\item{lpmodel}{The \code{lpmodel} object.}

\item{kappa}{The tuning parameter used in the second step of the two-step
procedure for obtaining the bounds subject to the shape constraints.
It can be any nonnegative number.}

\item{norm}{The norm used in the optimization problem. It can be either a
1-norm or a 2-norm. See the details section of
\code{\link[lpinfer]{estbounds}} for a list of accepted strings.}

\item{estimate}{A boolean variable to indicate whether the bounds should be
estimated or not.}

\item{solver}{The name of the linear and quadratic programming solver that
is used to obtain the solution to linear and quadratic programs.
The solvers supported by this package are \code{cplexAPI}, \code{gurobi},
\code{limSolve} and \code{Rcplex}.}
}
\value{
Returns the bounds subject to the shape constraints.
  \item{ub}{The upper bound.}
  \item{lb}{The lower bound.}
  \item{mincriterion}{The objective value of the first-stage problem (i.e.
    \code{\link[lpinfer]{mincriterion}}).}
  \item{est}{The indicator of whether estimation is involved in the
    estimation.}
  \item{call}{The function that has been called.}
  \item{norm}{The norm used.}
  \item{ub.status}{The status of the optimization problem for the upper
    bound.}
  \item{lb.status}{The status of the optimization problem for the lower
    bound.}
  \item{solver}{The name of the solver used.}
}
\description{
This function estimates the bounds subject to the shape
  constraints. The estimation is done via a two-step procedure with either
  the 1-norm or the 2-norm.
}
\details{
The following input for \code{norm} will be interpreted as the
  1-norm:
  \itemize{
    \item{\code{1} (\code{numeric})}
    \item{\code{"1"} (\code{string})}
    \item{\code{"L1"}}
    \item{\code{"one"}}
    \item{\code{"o"}}
    \item{\code{"taxicab"}}
  }
  The following input for \code{norm} will be interpreted as the 2-norm:
  \itemize{
    \item{\code{2} (\code{numeric})}
    \item{\code{"2"} (\code{string})}
    \item{\code{"L2"}}
    \item{\code{"two"}}
    \item{\code{"t"}}
    \item{\code{"e"}}
    \item{\code{"euclidean"}}
  }
}
\examples{
## ========================================================================= ##
##
##  Example for the estbounds function
##
##  This followings illustrate how the function can be used to estimate the
##  bounds subject to shape constraints using the missing data problem. 
##  Currently, only Gurobi is supported for the estbounds function with a
##  2-norm.
##
## ========================================================================= ##
rm(list = ls())

# ---------------- #
# Part 1: Load packages
# ---------------- #
library(lpinfer)
library(future)

# ---------------- #
# Part 2: Data and lpmodel preparation
# ---------------- #
source("./inst/example/dgp_missingdata.R")
J <- 5
N <- 1000
data <- missingdata_draw(J = J, n = N, seed = 1, prob.obs = .5)
lpmodel.full <- missingdata_lpm(J = J, info = "full", data = data)
lpmodel.twom <- missingdata_lpm(J = J, info = "mean", data = data)

tau <- sqrt(log(N)/N)
beta.tgt <- .2
kappa <- 1e-5

# ---------------- #
# Step 3: Run estbounds
# ---------------- #
# Example 1 - Compute the true bounds by setting estimate = FALSE
estb1 <- estbounds(data = data,
                   lpmodel = lpmodel.full,
                   kappa = kappa,
                   norm = 1,
                   solver = "gurobi",
                   estimate = FALSE)

# Example 2 - Estimated bounds with full-information approach and 1-norm
estb2 <- estbounds(data = data,
                   lpmodel = lpmodel.full,
                   kappa = kappa,
                   norm = 1,
                   solver = "gurobi",
                   estimate = TRUE)

# Example 3 - Estimated bounds with two-moments approach and 2-norm
estb3 <- estbounds(data = data,
                   lpmodel = lpmodel.twom,
                   kappa = kappa,
                   norm = 2,
                   solver = "gurobi",
                   estimate = TRUE)
}
