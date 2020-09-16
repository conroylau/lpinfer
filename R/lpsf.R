#' Obtain standard form of linear program
#'
#' @description This function is mainly a wrapper of the \code{simplex}
#'    function from the \code{boot} package. It returns the standard form
#'    of the constraints matrix of a linear program.
#'
#' @import boot
#'
#' @param A Matrix that represents the constraints.
#' @param b RHS vector that represents the constraints.
#' @param sense Sense of the constraints.
#' @param lb Lower bound of the variables.
#' @param obj Objective function.
#'
#' @return Returns the following four objects:
#'   \item{A}{Matrix representing the constraints in standard form.}
#'   \item{b}{RHS vector representing the constraints in standard form.}
#'   \item{lb}{Lower bound of the variables.}
#'   \item{obj}{Objective function.}
#'
#' @export
#'
standard.form <- function(A, b, sense, lb = NULL, obj = NULL) {
  # Obtain the matrices by the type of the inequality or equality
  A1 <- A[sense == "<=", ]
  A2 <- A[sense == ">=", ]
  A3 <- A[sense == "=", ]

  # Obtain the RHS vectors by the type of the inequality or equality
  b1 <- b[sense == "<="]
  b2 <- b[sense == ">="]
  b3 <- b[sense == "="]
  # Assign objective function
  a <- rep(1, ncol(A))

  # Obtain the constraints in standard form
  A.sf <- boot::simplex(a, A1, b1, A2, b2, A3, b3,
                        maxi = FALSE,
                        n.iter = 0)$A[, 1:(ncol(A) + nrow(rbind(A1, A2)))]

  # Obtain the new RHS vector
  b.sf <- c(b1, b2, b3)

  if (!is.null(lb)) {
    # Append the zero lower bound for the slack and surplus variables
    lb.sf <- c(lb, rep(0, nrow(rbind(A1, A2))))
  } else {
    lb.sf <- NULL
  }

  if (!is.null(obj)) {
    # Return the updated objective value
    obj.sf <- c(obj, rep(0, nrow(rbind(A1, A2))))
  } else {
    obj.sf <- NULL
  }

  return(list(A = A.sf,
              b = b.sf,
              lb = lb.sf,
              obj = obj.sf))
}

#' Obtains standard form of linear program for constraints in \code{lpmodel}
#'
#' @description This function is uses the \code{standard.form} function to
#'    convert a \code{lpmodel.natural} object into a \code{lpmodel} object.
#'
#' @param lpm.natural A \code{lpmodel.natural} object.
#'
#' @return Returns a \code{lpmodel} object.
#'   \item{lpmodel}{A \code{lpmodel} object.}
#'
#' @example ./inst/example/standard.lpmodel_example.R
#'
#' @export
#'
standard.lpmodel <- function(lpm.natural) {
  # ---------------- #
  # Step 1: Update the sense constraints and the inequality constraints
  # ---------------- #
  # Extract the objects from lpm.natural
  A.obs <- lpm.natural$A.obs
  A.shp <- lpm.natural$A.shp
  A.tgt <- lpm.natural$A.tgt
  beta.obs <- lpm.natural$beta.obs
  beta.shp <- lpm.natural$beta.shp
  sense.shp <- lpm.natural$sense.shp

  # Update the upper bounds
  if (!is.null(lpm.natural$x.ub)) {
    ub.temp <- list()
    ub.temp$A <- diag(length(lpm.natural$x.ub))
    ub.temp$b <- c(lpm.natural$x.ub)
    ub.temp$sense <- rep("<=", length(lpm.natural$x.ub))

    # Attach to the shape matrices
    A.shp <- rbind(A.shp, ub.temp$A)
    beta.shp <- c(beta.shp, ub.temp$b)
    sense.shp <- c(sense.shp, ub.temp$sense)
  }

  # Update the lower bounds
  if (!is.null(lpm.natural$x.lb)) {
    lb.temp <- list()
    lb.temp$A <- diag(length(lpm.natural$x.lb))
    lb.temp$b <- c(lpm.natural$x.lb)
    lb.temp$sense <- rep(">=", length(lpm.natural$x.lb))

    # Attach to the shape matrices
    A.shp <- rbind(A.shp, lb.temp$A)
    beta.shp <- c(beta.shp, lb.temp$b)
    sense.shp <- c(sense.shp, lb.temp$sense)
  }

  # ---------------- #
  # Step 2: Prepare models to be passed to the standard.form function
  # ---------------- #
  # Make the shape constraints in standard form
  if (length(sense.shp) > 0) {
    lpm.temp <- standard.form(A = A.shp,
                              b = matrix(beta.shp, ncol = 1, byrow = TRUE),
                              sense = matrix(sense.shp, ncol = 1, byrow = TRUE))
    lpm.new <- list()
    lpm.new$A.shp <- lpm.temp$A
    lpm.new$beta.shp <- lpm.temp$b

    # Add the zeros to the equality constraints
    k <- ncol(lpm.new$A.shp) - ncol(A.obs)
    lpm.new$A.obs <- cbind(A.obs, matrix(rep(0, k * nrow(A.obs)), ncol = k))
    lpm.new$A.tgt <- cbind(A.tgt, matrix(rep(0, k * nrow(A.tgt)), ncol = k))
  }

  # ---------------- #
  # Step 3: Create the new lpmodel object
  # ---------------- #
  lpm <- lpmodel(A.obs = lpm.new$A.obs,
                 A.tgt = lpm.new$A.tgt,
                 A.shp = lpm.new$A.shp,
                 beta.obs = matrix(c(beta.obs), ncol = 1),
                 beta.shp = matrix(c(lpm.new$beta.shp), ncol = 1))
  return(lpm)
}
