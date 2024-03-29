% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim.R
\name{lpsolveapi.optim}
\alias{lpsolveapi.optim}
\title{LP solver by \code{lpSolveAPI}}
\usage{
lpsolveapi.optim(Af, bf, nf, A, rhs, sense, modelsense, lb, weight = NULL, ...)
}
\arguments{
\item{Af}{The matrix that is involved in the objective function.}

\item{bf}{The vector that is involved in the objective function.}

\item{nf}{The number of observations in the data frame.}

\item{A}{The constraint matrix.}

\item{rhs}{The rhs vector for the linear constraints.}

\item{sense}{The sense of the linear constraints.}

\item{modelsense}{The indicator of whether the model is to max or min an
objective function.}

\item{lb}{The lower bound vector.}

\item{weight}{The weighting matrix.}

\item{...}{List of options to be passed to the Gurobi solver. This part is
optional.}
}
\value{
Returns the optimal point, the optimal value and the status of the
 solution.
 \item{objval}{The optimal value.}
 \item{x}{The optimal point.}
 \item{status}{The status of the optimization problem.}
}
\description{
This function computes the solution to the linear program
   using the \code{lpSolveAPI} package.
}
\details{
The package \code{lpSolveAPI} cannot be used to solve quadratic
  programs.
}
