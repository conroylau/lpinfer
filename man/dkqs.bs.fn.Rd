% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dkqs.R
\name{dkqs.bs.fn}
\alias{dkqs.bs.fn}
\title{Carries out one bootstrap replication for the \code{\link[lpinfer]{dkqs}}
procedure}
\usage{
dkqs.bs.fn(
  x,
  data,
  lpmodel,
  beta.obs.hat,
  beta.tgt,
  s.star.list,
  tau.list,
  solver,
  n,
  pbar,
  eval.count,
  n.bs,
  any.list,
  progress
)
}
\arguments{
\item{x}{This is either the list of indices that represent the bootstrap
replications, or the list of bootstrap components of the \code{lpmodel}
object passed from the user.}

\item{data}{An \code{data frame} or a \code{matrix}.}

\item{lpmodel}{An \code{lpmodel} object.}

\item{beta.obs.hat}{The value of sample \eqn{\hat{\bm{\beta}}_{\mathrm{obs}}}
from the \code{lpmodel} object.}

\item{beta.tgt}{The value to be tested.}

\item{s.star.list}{The list of values of
\eqn{\hat{\bm{s}}^\star \equiv \bm{A}_{\mathrm{obs}}\hat{\bm{x}}_n^\star}
in the cone-tightening procedure for each \eqn{\tau}.}

\item{tau.list}{The list of feasible parameters \eqn{\tau}.}

\item{solver}{The name of the linear and quadratic programming solver that
is used to obtain the solution to linear and quadratic programs.
The solvers supported by this package are \code{cplexAPI}, \code{gurobi},
\code{limSolve} and \code{Rcplex}.}

\item{n}{The sample size. This is only required if \code{data} is omitted
in the input.}

\item{pbar}{The progress bar object.}

\item{eval.count}{The count for the number of times the \code{future_map}
function has been called. If this object is zero, it means that the
\code{future_map} function is being called for the first time in this
subprocedure. Otherwise, it means that the \code{future_map} function
has been called for more than once. This situation typically refers to the
situations where there are some errors in the first time of the
replications.}

\item{n.bs}{The total number of replications to be conducted in this
procedure.}

\item{any.list}{This object contains information on whether any components
in the \code{lpmodel} object is a \code{list}.}

\item{progress}{The boolean variable for whether the progress bars should
be displayed. If it is set as \code{TRUE}, the progress bars will be
displayed while the code is running.}
}
\value{
Returns a list of output that are obtained from the DKQS
  procedure:
  \item{Ts}{A bootstrap test statistic.}
  \item{beta}{A bootstrap estimator.}
  \item{param}{A list of problematic parameters in the DKQS test.}
  \item{msg}{An error message (if applicable).}
}
\description{
This function carries out the one bootstrap replication of the
  DKQS procedure This function is used in the \code{\link[lpinfer]{dkqs.bs}}
  function via the \code{future_map} function.
}
