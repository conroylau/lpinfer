% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{bs.index}
\alias{bs.index}
\title{Auxiliary function to return the indices for bootstrap
replications}
\usage{
bs.index(R, R.eval, R.succ, maxR)
}
\arguments{
\item{R}{The number of bootstrap replications.}

\item{R.eval}{The number of bootstrap replications that has been evaluated.}

\item{R.succ}{The number of successful bootstrap replications.}

\item{maxR}{The maximum number of bootstrap replications.}
}
\value{
Returns the following three objects.
  \item{i0}{The starting index to be evaluated.}
  \item{i1}{The last index to be evaluated.}
}
\description{
This function is used to return the starting and ending
  indices for bootstrap replications.
}
