% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpsf.R
\name{standard.lpmodel}
\alias{standard.lpmodel}
\title{Obtains standard form of linear program for constraints in \code{lpmodel}}
\usage{
standard.lpmodel(lpm.natural)
}
\arguments{
\item{lpm.natural}{A \code{lpmodel.natural} object.}
}
\value{
Returns a \code{lpmodel} object.
  \item{lpmodel}{A \code{lpmodel} object.}
}
\description{
This function is uses the \code{standard.form} function to
   convert a \code{lpmodel.natural} object into a \code{lpmodel} object.
}
\examples{
## ========================================================================= ##
##
##  Example for the standard.lpmodel function
##
##  The following illustrates how to convert an lpmodel.natural object into an
##  lpmodel object.
##  
## ========================================================================= ##

### Step 1: Create an object in the `lpmodel.natural` class
# Obs
Aobs0 <- matrix(c(1, 2), nrow = 1)
bobs0 <- c(10)

# Shp
Ashp0 <- matrix(c(3, 4, 5, 6), nrow = 2, byrow = TRUE)
bshp0 <- matrix(c(15, 100))
sshp0 <- matrix(c(">=", "<="))

# Tgt
Atgt0 <- matrix(c(1, 1), nrow = 1)

# Upper bounds
xub0 <- c(200, 200)

# Lower bounds
xlb0 <- c(0.1, 0.1)

# Formulate the `lpmodel.natural` object
lpmn0 <- lpmodel.natural(A.obs = Aobs0,
                         A.shp = Ashp0,
                         A.tgt = Atgt0,
                         beta.obs = bobs0,
                         beta.shp = bshp0,
                         sense.shp = sshp0,
                         x.ub = xub0,
                         x.lb = xlb0)

### Step 2: Apply the `standard.lpmodel` function
lpm1 <- standard.lpmodel(lpmn0)
}
