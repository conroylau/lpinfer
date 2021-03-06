% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{check.norm}
\alias{check.norm}
\title{Check function: norm}
\usage{
check.norm(x, name.var)
}
\arguments{
\item{x}{The variable to be checked.}

\item{name.var}{The name of the variable.}
}
\value{
Nothing is returned.
}
\description{
This function checks whether the the norm used in the
  problem is a 1-norm or a 2-norm. If not, an error message is displayed.
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
  Capitalization is not an issue here as the text will be brought
  to the lower case.
}
