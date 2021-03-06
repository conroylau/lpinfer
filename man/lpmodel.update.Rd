% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpmodel.R
\name{lpmodel.update}
\alias{lpmodel.update}
\title{Combines deterministic components and one bootstrap estimate in
\code{lpmodel}}
\usage{
lpmodel.update(lpm.de, lpm.st)
}
\arguments{
\item{lpm.de}{The deterministic components of the \code{lpmodel} object.}

\item{lpm.st}{An \code{lpmodel} object that only contains one bootstrap
replication of the stochastic component(s). The deterministic component
is set as \code{NULL}.}
}
\value{
Returns an \code{lpmodel} object that combines the deterministic
  and stochastic component.
  \item{lpm.de}{An updated \code{lpmodel} object.}
}
\description{
This function is used in the bootstrap replications to combine
  the deterministic components in \code{lpmodel} to the stochastic component.
}
