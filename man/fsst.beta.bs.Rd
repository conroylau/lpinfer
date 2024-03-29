% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsst.R
\name{fsst.beta.bs}
\alias{fsst.beta.bs}
\title{Computing the bootstrap estimates of \code{beta.obs}}
\usage{
fsst.beta.bs(
  n,
  data,
  beta.obs.hat,
  lpmodel,
  R,
  maxR,
  progress,
  df.error,
  iseq,
  eval.count
)
}
\arguments{
\item{n}{The sample size. This is only required if \code{data} is omitted
in the input.}

\item{data}{An \code{data frame} or a \code{matrix}.}

\item{beta.obs.hat}{The sample estimator
\eqn{\widehat{\bm{\beta}}_{\mathrm{obs}, n}} based on the given
information in \code{lpmodel} (and \code{data} if applicable).}

\item{lpmodel}{An \code{lpmodel} object.}

\item{R}{The number of bootstrap replications.}

\item{maxR}{The maximum number of bootstrap replications to be considered in
case there are some errors.}

\item{progress}{The boolean variable for whether the progress bars should
be displayed. If it is set as \code{TRUE}, the progress bars will be
displayed while the code is running.}

\item{df.error}{A table showing the id of the bootstrap replication(s)
with error(s) and the corresponding error message(s).}

\item{iseq}{The list of indices or betas to iterate over.}

\item{eval.count}{The count for the number of times the \code{future_map}
function has been called. If this object is zero, it means that the
\code{future_map} function is being called for the first time in this
subprocedure. Otherwise, it means that the \code{future_map} function
has been called for more than once. This situation typically refers to the
situations where there are some errors in the first time of the
replications.}
}
\value{
Returns the bootstrap estimators.
  \item{beta.obs.bs}{A list of bootstrap estimators
  \eqn{\{\hat{\beta}_{{\rm obs}, n, b}\}^B_{b=1}}.}
  \item{df.error}{An updated table showing the id of the bootstrap
    replication(s) with error(s) and the corresponding error message(s).}
  \item{R.eval}{The number of bootstrap replications that has been
    conducted.}
  \item{R.succ}{The number of successful bootstrap replications.}
}
\description{
This function computes the bootstrap estimates of
  \eqn{\hat{\beta}_{{\rm obs}, n}}.
}
