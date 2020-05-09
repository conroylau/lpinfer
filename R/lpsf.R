#' Obtain standard form of linear program
#'
#' @description This function is mainly a wrapper of the \code{simplex}
#'    function from the \code{boot} package. It returns the standard form
#'    of the constraints matrix in a linear program.
#'
#' @import boot
#'
#' @param A Matrix that represents the constraints
#' @param b RHS vector that represents the constraints
#' @param sense Sense of the constraints
#' @param lb Lower bound of the vectors
#' @param obj Objective function
#'
#' @return Returns the bounds subject to the shape constraints.
#'   \item{A}{Matrix representing the constraints in standard form}
#'   \item{b}{RHS vector representing the constraints in standard form}
#'   \item{lb}{Lower bound of the variables}
#'   \item{obj}{Objective function}
#'
#' @export
#'
standard.form <- function(A, b, sense, lb = NULL, obj = NULL){
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

  if (!is.null(lb)){
    # Append the zero lower bound for the slack and surplus variables
    lb.sf <- c(lb, rep(0, nrow(rbind(A1, A2))))
  } else {
    lb.sf <- NULL
  }

  if (is.null(obj)){
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
