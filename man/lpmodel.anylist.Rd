% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpmodel.R
\name{lpmodel.anylist}
\alias{lpmodel.anylist}
\title{Check if there is any list in the \code{lpmodel} object}
\usage{
lpmodel.anylist(lpmodel)
}
\arguments{
\item{lpmodel}{A \code{lpmodel} object.}
}
\value{
Returns the following objects:
  \item{list}{A boolean variable that indicates whether there is any object
     with the class \code{list} in the \code{lpmodel} object.}
  \item{name}{The names of the components with class \code{list} in the
    \code{lpmodel} object.}
  \item{consol}{A consolidated and updated \code{lpmodel} object that lists
    the \code{lpmodel} by observation instead of by component (if
    applicable).}
  \item{len}{The length of the component in the \code{lpmodel} object that
    is a list (if applicable).}
}
\description{
This function checks if there is any components in the
  \code{lpmodel} object with class \code{list}.
}
