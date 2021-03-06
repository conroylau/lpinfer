% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chorussell.R
\name{chorussell.simp}
\alias{chorussell.simp}
\title{Simplifies the candidates to be considered in
\code{\link[lpinfer]{chorussell}}}
\usage{
chorussell.simp(lb.can1, lb.can2, ub.can1, ub.can2, progress)
}
\arguments{
\item{lb.can1}{The vector of values that corresponds to
\eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm lb} - \hat{\theta}_{\rm lb}\right)}.}

\item{lb.can2}{The vector of values that corresponds to
\eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm lb} - \hat{\theta}_{\rm lb} -
\Delta \right)}.}

\item{ub.can1}{The vector of values that corresponds to
\eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm ub} - \hat{\theta}_{\rm ub}\right)}.}

\item{ub.can2}{The vector of values that corresponds to
\eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm ub} - \hat{\theta}_{\rm ub} +
\Delta \right)}.}

\item{progress}{The boolean variable for whether the progress bars should
be displayed. If it is set as \code{TRUE}, the progress bars will be
displayed while the code is running.}
}
\value{
Returns the list of updated candidates of lower bounds and upper
  bounds.
  \item{lb}{The updated list of candidates of lower bounds.}
  \item{ub}{The updated list of candidates of upper bounds.}
}
\description{
This function simplifies the list of candidates to be
  considered in the optimization problem in the
  \code{\link[lpinfer]{chorussell}} function. In particular, because
  \deqn{\mathbf{1}\left[\sqrt{n}\left(\hat{\theta}_{\rm lb}^b -
  \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb}\right]
  \geq
  \mathbf{1}\left[\sqrt{n}\left(\hat{\theta}^b_{\rm lb} -
  \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb} \quad \mathrm{ and } \quad
  -c_{\rm ub} \leq \sqrt{n} \left(\hat{\theta}^b_{\rm ub} -
  \hat{\theta}_{\rm ub} + \Delta\right)
  \right],}
  the values of \eqn{c_{\rm lb}} that satisfy
  \deqn{\frac{1}{B}\sum^B_{b=1} \mathbf{1}
  \left[\sqrt{n}\left(\hat{\theta}_{\rm lb}^b -
  \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb}\right] < 1 - \alpha}
  will be removed from the set of the values under consideration. Similarly,
  the list of \eqn{c_{\rm ub}} that satisfy
  \deqn{\frac{1}{B}\sum^B_{b=1} \mathbf{1}
  \left[- \sqrt{n}\left(\hat{\theta}_{\rm ub}^b -
  \hat{\theta}_{\rm ub}\right) \leq c_{\rm lb}\right] < 1 - \alpha}
  will be removed from the set of the values under consideration.
}
