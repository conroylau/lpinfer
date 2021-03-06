% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsst.R
\name{sigma.summation}
\alias{sigma.summation}
\title{Computes the asymptotic variance estimator}
\usage{
\method{sigma}{summation}(n, beta.bs.list, progress, eval.count)
}
\arguments{
\item{n}{Sample size.}

\item{beta.bs.list}{A list of bootstrap estimators
\eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.}
}
\value{
Returns the estimator of the asymptotic variance.
    \item{sigma.mat}{The estimator of the asymptotic variance.}
}
\description{
Based on the bootstrap estimates
  \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}, this function computes the
  asymptotic variance estimator of the bootstrap estimator, i.e.
  \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b -
  \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b -
  \widehat{\bm{\beta}}\right)'.}
 This function supports parallel programming via the \code{future.apply}
 package.
}
