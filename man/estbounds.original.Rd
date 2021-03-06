% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estbounds.R
\name{estbounds.original}
\alias{estbounds.original}
\title{Computes the true bounds subjected to shape constraints}
\usage{
estbounds.original(data, lpmodel, original.sense, solver)
}
\arguments{
\item{data}{An \code{data frame} or a \code{matrix}.}

\item{lpmodel}{The \code{lpmodel} object.}

\item{original.sense}{The sense of the constraints for the true bounds.}

\item{solver}{The name of the linear and quadratic programming solver that
is used to obtain the solution to linear and quadratic programs.
The solvers supported by this package are \code{cplexAPI}, \code{gurobi},
\code{limSolve} and \code{Rcplex}.}
}
\value{
Returns the solution to the linear program.
 \item{objval}{The optimal value.}
 \item{status}{The status of the linear program.}
}
\description{
This function computes the true bounds subjected to the shape
   constraints without approximation.
}
