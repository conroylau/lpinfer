% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{check.betatgt}
\alias{check.betatgt}
\title{Check function: check if \code{beta.tgt} is within the logical bound}
\usage{
check.betatgt(data, lpmodel, beta.tgt, solver)
}
\arguments{
\item{data}{An \code{data frame} or a \code{matrix}.}

\item{lpmodel}{An \code{lpmodel} object.}

\item{beta.tgt}{The value to be tested.}

\item{solver}{A linear or quadratic programming solver. The exact solver
that is supported depends on the test chosen.}
}
\value{
Returns the following objects:
  \item{inout}{The variable that equals 1 if it is within the logical bound.
  This equals 0 otherwise.}
  \item{lb}{The logical lower bound.}
  \item{ub}{The logical upper bound.}
}
\description{
This function checks whether the parameter \code{beta.tgt}
  is within the logical bounds. If it is not within the logical bound,
  then reject immediately.
}
