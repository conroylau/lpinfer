#' Check function: data frame
#'
#' @description This function coerce the data provided from the user as a
#'   \code{data.frame}.
#'
#' @param data A \code{data} object.
#'
#' @return Returns the object as a \code{data.frame}.
#'
#' @export
#'
check.dataframe <- function(data) {
  # Check data
  data <- as.data.frame(data)

  # Return updated data
  return(data)
}

#' Check function: passing data to function
#'
#' @description This function checks whether there will be any error if an
#'   object of class \code{data.frame} is passed to the user-defined function.
#'
#' @param data The data object in the class \code{data.frame}.
#' @param f A user-defined function in \code{lpmodel}.
#' @param lpmodel.comp The name of the \code{lpmodel} object.
#'
#' @details An error message will be displayed.
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.datafunction <- function(data, f, lpmodel.comp) {
  # General warning message
  display.msg <- sprintf(paste0("The function defined in the '%s' component ",
                                "for 'lpmodel' needs to accept data of class ",
                                "'data.frame'."), lpmodel.comp)

  # try-catch
  tryCatch(
    expr = {
      f(data)
    },
    error = function(e) {
      stop(display.msg)
    },
    warning = function(w) {
      stop(display.msg)
    },
    finally = {
    }
  )
}

#' Check function: positive integer
#'
#' @description This function checks whether the class of the variable
#'   is \code{numeric}, has length 1 and is a positive integer. If not, an
#'   error message is displayed.
#'
#' @param x The variable to be checked.
#' @param name.var The name of the variable.
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.positiveinteger <- function(x, name.var) {
  if (!is.numeric(x)) {
    # Call the general error message function
    check.errormsg(name.var, "a positive integer")
  } else if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0 & x %% 1 == 0)
             == FALSE) {
    # Call the general error message function
    check.errormsg(name.var, "a positive integer")
  }

  return(x)
}

#' Check function: nonnegative number
#'
#' @description This function checks whether variable satisfies the following
#'   requirements:
#'   \itemize{
#'     \item{The class of the variable is \code{numeric}.}
#'     \item{The length of the variable is 1.}
#'     \item{The variable is nonnegative.}
#'   }
#'   If not, an error message is displayed.
#'
#' @param x The variable to be checked.
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.nonnegative <- function(x, name.var) {
  if (!is.numeric(x)) {
    # Call the general error message function
    check.errormsg(name.var, "a nonnegative number")
  } else if ((is.numeric(x) == TRUE & length(x) == 1 & x >= 0) == FALSE) {
    # Call the general error message function
    check.errormsg(name.var, "a nonnegative number")
  }

  return(x)
}

#' Check function: positive number
#'
#' @description This function checks whether variable satisfies the following
#'   requirements:
#'   \itemize{
#'     \item{The class of the variable is \code{numeric}.}
#'     \item{The length of the variable is 1.}
#'     \item{The variable is strictly positive.}
#'   }
#'   If not, an error message is displayed.
#'
#' @param x The variable to be checked.
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.positive <- function(x, name.var) {
  if (!is.numeric(x)) {
    # Call the general error message function
    check.errormsg(name.var, "a strictly positive number")
  } else if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0) == FALSE) {
    # Call the general error message function
    check.errormsg(name.var, "a strictly positive number")
  }

  return(x)
}

#' Check function: numeric
#'
#' @description This function checks whether the class of the variable
#'   is \code{numeric} and has length 1. If not, an error message is displayed.
#'
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.numeric <- function(x, name.var) {
  if ((is.numeric(x) == TRUE & length(x) == 1) == FALSE) {
    stop(sprintf("The class of the variable '%s' has to be numeric.",
                 name.var),
         call. = FALSE)
  }
}

#' Check function: boolean variable
#'
#' @description This function checks whether the variable is boolean. If not,
#'   an error message is displayed.
#'
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.boolean <- function(x, name.var) {
  if (!(x == TRUE | x == FALSE)) {
    # Call the general error message function
    check.errormsg(name.var, "a Boolean expression")
  }
}

#' Check function: range of a variable
#'
#' @description This function checks whether the variable is within a
#'   certain interval. If not, an error message is displayed.
#'
#' @param left.type The type of the left interval (\code{open} or
#'   \code{closed}).
#' @param left The value of lower bound.
#' @param right.type The type of the right interval (\code{open} or
#'   \code{closed}).
#' @param right the value of the upper bound.
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.numrange <- function(x, name.var, left.type, left, right.type, right) {
  # ---------------- #
  # Step 1: Check if the number is numeric
  # ---------------- #
  check.numeric(x, name.var)

  # ---------------- #
  # Step 2: Check if the number is within the range
  # ---------------- #
  # General message
  msg.interval <- paste0("The variable '%s' has to be inside the interval ",
                         "%s%s, %s%s.")

  # Start the check
  if (left.type == "open" & right.type == "open") {
    if ((x > left & x < right) == FALSE) {
      stop(sprintf(msg.interval, name.var, "(", left, right, ")"),
           call. = FALSE)
    }
  } else if (left.type == "open" & right.type == "closed") {
    if ((x > left & x <= right) == FALSE) {
      stop(sprintf(msg.interval, name.var, "(", left, right, "]"),
           call. = FALSE)
    }
  } else if (left.type == "closed" & right.type == "open") {
    if ((x >= left & x < right) == FALSE) {
      stop(sprintf(msg.interval, name.var, "[", left, right, ")"),
           call. = FALSE)
    }
  } else if (left.type == "closed" & right.type == "closed") {
    if ((x >= left & x <= right) == FALSE) {
      stop(sprintf(msg.interval, name.var, "[", left, right, "]"),
           call. = FALSE)
    }
  }
}

#' Check function: norm
#'
#' @description This function checks whether the the norm used in the
#'   problem is a 1-norm or a 2-norm. If not, an error message is displayed.
#'
#' @details The following input for \code{norm} will be interpreted as the
#'   1-norm:
#'   \itemize{
#'     \item{\code{1} (\code{numeric})}
#'     \item{\code{"1"} (\code{string})}
#'     \item{\code{"L1"}}
#'     \item{\code{"one"}}
#'     \item{\code{"o"}}
#'     \item{\code{"taxicab"}}
#'   }
#'   The following input for \code{norm} will be interpreted as the 2-norm:
#'   \itemize{
#'     \item{\code{2} (\code{numeric})}
#'     \item{\code{"2"} (\code{string})}
#'     \item{\code{"L2"}}
#'     \item{\code{"two"}}
#'     \item{\code{"t"}}
#'     \item{\code{"e"}}
#'     \item{\code{"euclidean"}}
#'   }
#'   Capitalization is not an issue here as the text will be brought
#'   to the lower case.
#'
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.norm <- function(x, name.var) {
  # Bring the variable to lower case
  x <- tolower(x)

  if (x %in% c(1, "1", "l1", "one", "o", "taxicab")) {
    # Case 1: user provided an input that corresponds to 1-norm
    norm <- 1
  } else if (x %in% c(2, "2", "l2", "two", "t", "e", "euclidean")) {
    # Case 2: user provided an input that corresponds to 2-norm
    norm <- 2
  } else {
    stop(gsub("\\s+", " ",
              paste0("Only 1-norm and 2-norm are supported in this function. ",
                     "For 1-norm, please use one of the followings: 1, '1', ",
                     "'l1', 'one', 'o' or 'taxicab'. For 2-norm, please use ",
                     "one of the followings: 2, '2', 'l2', 'two', 't', 'e', ",
                     "'euclidean'.")),
         call. = FALSE)
  }

  # Return the updated norm
  return(norm)
}

#' Check function: function
#'
#' @description This function checks whether the class of the function
#'   provided by the user is \code{function}, and whether the output
#'   produced by the function is \code{numeric} with the correct dimension.
#'   If not, an error message is displayed.
#'
#' @param f The function to be tested.
#' @param A A matrix that contains the information about the dimension for the
#'    output of f.
#' @param data A data frame.
#' @param name.A The name of the matrix \code{A}.
#' @param mat.type The type of the matrix to be checked.
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Returns the output of the function in the format of \code{numeric}
#'   and the correct dimension.
#'   \itemize{
#'      \item If \code{mat.type} is \code{col}, then a column vector is
#'      returned.
#'      \item If \code{mat.type} is \code{matrix}, then a square matrix is
#'      returned.}
#'
#' @export
#'
check.func <- function(f, A, data, name.var, name.A, mat.type) {
  # ---------------- #
  # Step 1: General messages
  # ---------------- #
  msg.class <- paste0("The %s of '%s' has to be %s.")

  # ---------------- #
  # Step 2: Start the check
  # ---------------- #
  # Check the class of the function
  if (!inherits(f, "function")) {
    stop(sprintf(msg.class, "input", name.var, "a function"),
         call. = FALSE)
  } else {
    out <- f(data)
    out <- asmat(out)
    # Check if the output is numeric
    if (is.numeric(out[, 1]) == FALSE) {
      stop(sprintf(msg.class, "output", name.var, "numeric"),
           call. = FALSE)
    } else{
      if (mat.type == "col") {
        # If the output has to be a column vector
        if (dim(out)[2] != 1) {
          stop(sprintf(msg.class, "output", name.var, "a column vector"),
               call. = FALSE)
        } else if (dim(out)[1] != dim(A)[1]) {
          stop(sprintf(paste0("The number of rows in the output of '%s' need ",
                              "to be the same as the number of rows in the ",
                              "matrix '%s'."), name.var, name.A),
               call. = FALSE)
        }
      } else if (mat.type == "square") {
        # If the output has to be a square matrix
        if (nrow(out) != ncol(out)) {
          stop(sprintf(msg.class, "output", name.var, "a square matrix"),
               call. = FALSE)
        }
        if ((nrow(out) != nrow(A)) | (ncol(out) != nrow(A))) {
          stop(sprintf("The number of rows and columns for the output ",
                       "of '%s' need to be equal to the number of rows in ",
                       "the matrix '%s'.", name.var, name.A),
               call. = FALSE)
        }
      }
    }
  }

  return(out)
}

#' Check function: constraint matrix and the corresponding rhs vector
#'
#' @description This function checks the constraint matrix \eqn{\bm{A}} and
#'    the corresponding vector \eqn{\bm{\beta}}.
#'
#' @param A The constraint matrix.
#' @param b The corresponding rhs vector.
#' @param Aname The name of the constraint matrix.
#' @param bname The name of the corresponding rhs vector.
#'
#' @return Returns the updated matrix and rhs vector or prints a
#'    stop message.
#'    \item{A}{The updated constraint matrix.}
#'    \item{b}{The updated rhs vector.}
#'
#' @export
#'
check.Ab <- function(A, b, Aname, bname) {
  if (is.null(A) + is.null(b) == 1) {
    # ---------------- #
    # Step 1: Check that if A and b must be both NULL or both non-NULL
    # ---------------- #
    msg.temp <- paste0("'%s' is NULL but '%s' is not NULL. Please ",
                       "ensure that either they are both NULL or they are ",
                       "both valid matrices.")
    if (is.null(A) == TRUE) {
      stop(sprintf(msg.temp, Aname, bname))
    } else {
      stop(sprintf(msg.temp, bname, Aname))
    }
  } else if (is.null(A) + is.null(b) == 2) {
    # ---------------- #
    # Step 2: Both A and b are NULL. Nothing else to do
    # ---------------- #
  } else {
    # ---------------- #
    # Step 3: Checks for the case where both A and b are non-NULL
    # ---------------- #
    matrix.names <- c(Aname, bname)
    matrix.list <- list(A, b)
    for (i in 1:2) {
      ## Part 1: Check the format of the matrices
      if (!(inherits(matrix.list[[i]], "data.frame") |
            inherits(matrix.list[[i]], "matrix") |
            inherits(matrix.list[[i]], "numeric"))) {
        stop(gsub("\\s+", " ",
                  paste0("The argument '", matrix.names[i],
                  "' must either be a data.frame, data.table, or matrix.")),
             call. = FALSE)
      } else {
        # Ensure the variable is in matrix form
        matrix.list[[i]] <- asmat(matrix.list[[i]])

        ## Part 2: Check whether the matrices are numeric
        if (is.numeric(matrix.list[[i]]) == FALSE) {
          stop(paste0("The argument '", matrix.names[i],
                      "' has to be numeric."), call. = FALSE)
        }
      }
    }

    ## Part 3: Ensure that beta is a column vector. This also captures the
    ## case of b being a scalar because it is turned into a matrix in the
    ## condition.
    if (dim(as.matrix(b))[2] != 1) {
      stop(sprintf("The argument '%s' has to be a column vector.", bname),
           call. = FALSE)
    }

    ## Part 4: Ensure that the number of rows of A and b are identical
    if (nrow(A) != length(b)) {
      stop(sprintf(paste0("The number of rows of '%s' has to be equal to ",
                          "the number of rows of '%s."), Aname, bname),
           call. = FALSE)
    }

    # ---------------- #
    # Step 4: Update class
    # ---------------- #
    # Ensure that both A and b to ensure that they are both matrices
    A <- matrix.list[[1]]
    b <- matrix.list[[2]]
  }

  # ---------------- #
  # Step 5: Return results
  # ---------------- #
  return(list(A = A,
              b = b))
}

#' Check function: solvers
#'
#' @description This function checks the solver used is supported and is
#'   appropriate for the testing function. If not, an error message is
#'   displayed.
#'
#' @param norm The norm used (if applicable).
#' @param qc An indicator of whether the problem consists of a quadratic
#'   constraint.
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Returns the function name that corresponds to the solver for the
#'   problem.
#'    \item{solver}{The function for the solver used.}
#'    \item{solver.name}{The name of the solver used in lower case.}
#'
#' @export
#'
check.solver <- function(x, name.var, norm = 2 , qc = FALSE) {
  # ---------------- #
  # Step 1: Preparation and hard-coded information
  # ---------------- #
  # Only change x to lower case if it is non-null. Otherwise, keep it as null
  if (!is.null(x)) {
    x <- tolower(x)
  }

  # Package recommendation messages
  gurobi.msg <- "'gurobi' (version 8.1-1 or later)"
  cplexapi.msg <- "'cplexAPI' (version 1.3.3 or later)"
  rcplex.msg <- "'Rcplex' (version 0.3-3 or later)"
  limsolve.msg <- "'limSolve' (version 1.5.6 or later)"
  lpsolveapi.msg <- "lpSolveAPI (version 5.5.2.0 or later)"

  # General error message
  msg.notcompat <- paste0("This function with a %s-norm in the estimation ",
                          "is not compatible with '%s'.")
  msg.install <- "Please install one of the following packages: %s"

  # Message for installing packages
  ## 1-norm
  msg.l1packages <- sprintf(msg.install,
                            paste0(gurobi.msg, "; ",
                                   cplexapi.msg, "; ",
                                   rcplex.msg, "; ",
                                   limsolve.msg, "; ",
                                   lpsolveapi.msg, "."))
  ## 2-norm with non-quadratically constrained QP
  msg.l2packages <- sprintf(msg.install,
                            paste0(gurobi.msg, "; ",
                                   cplexapi.msg, "; ",
                                   rcplex.msg, "; ",
                                   limsolve.msg, "."))
  ## 2-norm with quadratically constrained QP
  msg.l2norm.qc <- paste0("This function with a 2-norm in the estimation ",
                          "procedure and a quadratically constrained ",
                          "quadratic program is only compatible with  ",
                          "'gurobi'. Please install ", gurobi.msg, ".")

  # ---------------- #
  # Step 2a: If no solver name is provided by the user
  # ---------------- #
  if (is.null(x) == TRUE) {
    # If 'gurobi' is installed, the 'gurobi' solver will be used for 1-norm &
    # 2-norm
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE) {
      solver = gurobi.optim
      x <- "gurobi"
    } else if (norm == 1) {
      # If 1-norm is used, other solvers will be checked
      if (requireNamespace("limSolve", quietly = TRUE) == TRUE) {
        solver <- limsolve.optim
        x <- "limSolve"
      } else if (requireNamespace("Rcplex", quietly = TRUE) == TRUE) {
        solver <- rcplex.optim
        x <- "Rcplex"
      } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE) {
        solver <- cplexapi.optim
        x <- "cplexAPI"
      } else if (requireNamespace("lpSolveAPI", quietly = TRUE) == TRUE) {
        solver <- lpsolveapi.optim
        x <- "lpSolveAPI"
      }
    } else {
      if (norm == 1) {
        stop(gsub("\\s+", " ", msg.l1packages), call. = FALSE)
      } else if (norm == 2 & isFALSE(qc)) {
        stop(gsub("\\s+", " ", msg.l2packages), call. = FALSE)
      } else if (norm == 2 & isTRUE(qc)) {
        stop(gsub("\\s+", " ", msg.l2norm.qc), call. = FALSE)
      }
    }
  } else if (x == "gurobi") {
    # ---------------- #
    # Step 2b: If a solver name is provided by the user
    # ---------------- #
    ## Case 1: If the user specified the solver as 'gurobi'
    solver <- gurobi.optim
    x <- "gurobi"
  } else if (x == "limsolve" & !(norm == 2 & isTRUE(qc))) {
    ## Case 2: If the user specified the solver as 'limSolve'
    solver <- limsolve.optim
    x <- "limSolve"
  } else if (x == "rcplex" & !(norm == 2 & isTRUE(qc))) {
    ## Case 3: If the user specified the solver as 'rcplex'
    solver <- rcplex.optim
    x <- "Rcplex"
  } else if (x == "cplexapi" & !(norm == 2 & isTRUE(qc))) {
    ## Case 4: If the user specified the solver as 'cplexapi'
    solver <- cplexapi.optim
    x <- "cplexAPI"
  } else if (x == "lpsolveapi" & norm == 1) {
    ## Case 5: If the user specified the solver as 'lpsolveapi'
    solver <- lpsolveapi.optim
    x <- "lpSolveAPI"
  } else {
    ## Case 6: If the user specified a solver that is not compatible
    if (norm == 1) {
      stop(gsub("\\s+", " ", paste(sprintf(msg.notcompat, 1, x),
                                   msg.l1packages)), call. = FALSE)
    } else if (norm == 2 & isTRUE(qc)) {
      stop(gsub("\\s+", " ", msg.l2norm.qc), call. = FALSE)
    } else if (norm == 2 & isFALSE(qc)) {
      stop(gsub("\\s+", " ", paste(sprintf(msg.notcompat, 2, x),
                                   msg.l2packages)), call. = FALSE)
    }
  }

  # ---------------- #
  # Step 3: Returns the function name that corresponds to the solver
  # ---------------- #
  return(list(solver = solver,
              solver.name = x))
}

#' Check function: \code{lpmodel}
#'
#' @description This function checks if the object \code{lpmodel} is in the
#'    correct format.
#'
#' @param lpmodel An \code{lpmodel} object.
#' @param name.var The name of the \code{lpmodel} object.
#' @param A.tgt.cat The category of the \code{A.tgt} object.
#' @param A.obs.cat The category of the \code{A.obs} object.
#' @param A.shp.cat The category of the \code{A.shp} object.
#' @param beta.obs.cat The category of the \code{beta.obs} object.
#' @param beta.shp.cat The category of the \code{beta.shp} object.
#' @param is.estbounds A boolean variable that indicates whether the test
#'   function being called is \code{\link[lpinfer]{estbounds}}. If this
#'   function is being called, then it means that when \code{data} is
#'   \code{NULL}, then each component of the \code{lpmodel} object is not
#'   required to be a \code{list}. Otherwise, the \code{lpmodel} object needs
#'   to contain the bootstrap estimates if \code{data} is \code{NULL}.
#' @inheritParams dkqs
#'
#' @details In each of the testing procedures, there are six possible
#' categories for the each of the objects:
#' \itemize{
#'   \item{\code{not_used}: This refers to the case where the object is not
#'    used in the function.}
#'   \item{\code{matrix}: This refers to the case where the object has to be
#'    a matrix.}
#'  \item{\code{function_mat}: This refers to the case where
#'    the object has to be a function that produces a matrix.}
#'  \item{\code{list}: This refers to the case where the object is a list.}
#'  \item{\code{function_obs_var}: This refers to the case where
#'    the object is a function that produces a list that contains a matrix
#'    and a vector. This is typically the case for \code{beta.obs} when
#'    the testing procedure requires both the observed value of
#'    \code{beta.obs} and the estimator of the asymptotic variance.}
#'  \item{\code{function_obs_var_bs}: This is essentially the same as the
#'    \code{function_obs_var} category. However, when the variance matrix is
#'    not provided, it will be estimated by standard nonparametric bootstrap
#'    so there is no problem in case the variance matrix is not provided.}
#' }
#'  Each object can belong to one of more categories.
#'
#' @return Returns the updated \code{lpmodel} object.
#'    \item{lpmodel}{An updated \code{lpmodel} object.}
#'
#' @export
#'
check.lpmodel <- function(data, lpmodel, name.var, A.tgt.cat, A.obs.cat,
                          A.shp.cat, beta.obs.cat, beta.shp.cat, R,
                          is.estbounds = FALSE) {
  # ---------------- #
  # Step 1: Check if lpmodel is a list
  # ---------------- #
  if (is.null(lpmodel)) {
    stop("The 'lpmodel' object is required.")
  }

  if (!inherits(lpmodel, "lpmodel")) {
    # Call the general error message function
    check.errormsg(name.var, "an object in the 'lpmodel' class")
  }

  # ---------------- #
  # Step 2: Check if lpmodel includes the bootstrap estimates if data is NULL
  # (only check if lpmodel is deterministic if the function is estbounds)
  # ---------------- #
  # If data is NULL, then at least one of the components in lpmodel
  # has to be a list to represent the bootstrap replications except for
  # 'estbounds' where only deterministic components are required.
  error.msg.det <- paste0("When 'data' is 'NULL', each component in the ",
                          "'lpmodel' object has to be deterministic.")
  if (isFALSE(is.estbounds)) {
    if (is.null(data)) {
      notlist <- 0
      for (i in seq_along(lpmodel)) {
        if (inherits(lpmodel[[i]], "function")) {
          stop(error.msg.det)
        } else if (!inherits(lpmodel[[i]], "list")) {
          notlist <- notlist + 1
        }
      }

      # If none of the components inside the 'lpmodel' object is a list,
      # return an error message
      if (notlist == 0) {
        stop(paste0("When 'data' is 'NULL', the 'lpmodel' object needs ",
                    "to contain the bootstrap estimates."))
      }
    }
  } else {
    for (i in seq_along(lpmodel)) {
      if (is.null(data)) {
        if (inherits(lpmodel[[i]], "function")) {
          stop(error.msg.det)
        }
      }
    }
  }

  # ---------------- #
  # Step 3: Check each of the objects (treat them as matrices)
  # and check if beta.obs and beta.shp are matrices
  # ---------------- #
  if (!("not_used" %in% A.tgt.cat)) {
    A.tgt.return <- check.lpobjects(data, lpmodel$A.tgt, "A.tgt", A.tgt.cat, R)
  }
  if (!("not_used" %in% A.obs.cat)) {
    A.obs.return <- check.lpobjects(data, lpmodel$A.obs, "A.obs", A.obs.cat, R)
  }
  if (!("not_used" %in% A.shp.cat)) {
    A.shp.return <- check.lpobjects(data, lpmodel$A.shp, "A.shp", A.shp.cat, R)
  }
  if (!("not_used" %in% beta.obs.cat)) {
    beta.obs.return <- check.lpobjects(data, lpmodel$beta.obs, "beta.obs",
                                       beta.obs.cat, R)
    check.vector(beta.obs.return$sample, "beta.obs", FALSE)
  }
  if (!("not_used" %in% beta.shp.cat)) {
    beta.shp.return <- check.lpobjects(data, lpmodel$beta.shp, "beta.shp",
                                       beta.shp.cat, R)
    check.vector(beta.shp.return$sample, "beta.shp", FALSE)
  }

  # ---------------- #
  # Step 4: Check whether the dimension matches
  # ---------------- #
  # General message telling the user to provide the "A" matrix and the "beta"
  # vector with the same number of rows
  msg.row <- paste0("The objects '%s' and '%s' in 'lpmodel' need ",
                    "to have the same number of rows.")

  # Start the check
  if (!("not_used" %in% A.obs.cat)) {
    if (nrow(A.obs.return$sample) != NROW(beta.obs.return$sample)) {
      stop(sprintf(msg.row, "A.obs", "beta.obs"))
    }
  }
  if (!("not_used" %in% A.shp.cat)) {
    if (nrow(A.shp.return$sample) != NROW(beta.shp.return$sample)) {
      stop(sprintf(msg.row, "A.shp", "beta.shp"))
    }
  }

  return(lpmodel)
}

#' Check function: matrices and vectors in \code{lpmodel}
#'
#' @description This function checks if the matrix objects in \code{lpmodel}
#'    are in the correct format.
#'
#' @inheritParams dkqs
#' @param mat The matrix object in \code{lpmodel}.
#' @param mat.name The name of the matrix object in \code{lpmodel}.
#' @param mat.cat The category of the matrix object.
#'
#' @details See the details section for the function
#'   \code{\link[lpinfer]{check.lpmodel}} for more details of each category.
#'
#' @return Returns two objects:
#'    \item{mat}{The updated object in \code{lpmodel}.}
#'    \item{sample}{A sample object that is used to compare the dimension.}
#'
#' @export
#'
check.lpobjects <- function(data, mat, mat.name, mat.cat, R) {
  if (!is.null(mat)) {
    # ---------------- #
    # Step 1: Check error indicator
    # ---------------- #
    err.ind <- NULL

    # Category 1: Check if it is a single matrix or a sparseMatrix
    if ("matrix" %in% mat.cat) {
      mat.return <- check.matrix(mat, mat.name, mat.cat, FALSE)
      if (mat.return$err.ind != 1) {
        return(list(mat = mat.return$mat.update,
                    sample = mat.return$mat.update))
      } else {
        err.ind <- c(err.ind, 1)
      }
    }

    # Category 2: Check if it is a function
    if ("function_mat" %in% mat.cat) {
      if (inherits(mat, "function")) {
        sample.return <- lpmodel.beta.eval(data, mat, 1)

        # Check whether the function can accept 'data' in the 'data.frame'
        # format
        check.datafunction(data, mat, mat.name)

        return(list(mat = mat,
                    sample = sample.return$beta.obs))
      } else if (length(mat.cat) == 1) {
        stop(sprintf(paste0("The object '%s' in 'lpmodel' needs to be a ",
                            "function."),
                     mat.name),
             call. = FALSE)
      } else {
        err.ind <- c(err.ind, 2)
      }
    }

    # Category 3: Check if it is a list
    if ("list" %in% mat.cat) {
      if (inherits(mat, "list")) {
        # Check if the length of the list is R+1
        if (length(mat) != (R + 1)) {
          stop(sprintf(paste0("The object '%s' in 'lpmodel' needs to have ",
                              "exactly %s elements"),
                              mat.name, R + 1),
               call. = FALSE)
        }

        # Loop through each element of the list to check if they are matrices
        df.dim <- data.frame(matrix(vector(), nrow = R + 1))
        for (i in 1:R) {
          mat.return <- check.matrix(mat[[i + 1]], mat.name, mat.cat, TRUE)
          # Check if all objects inside the list has the same dimension
          df.dim[i, 1] <- mat.return$dim[1]
          df.dim[i, 2] <- mat.return$dim[2]
          if (i > 1) {
            if ((df.dim[i, 1] != df.dim[i - 1, 1]) |
                (df.dim[i, 1] != df.dim[i - 1, 1])) {
              stop(sprintf(paste0("The dimension of the objects inside the ",
                                  "list '%s' in 'lpmodel' need to have the ",
                                  "same dimension."),
                           mat.name),
                   call. = FALSE)
            }
          }
        }
        return(list(mat = mat,
                    sample = mat.return[[1]]))
      } else if (length(mat.cat) == 1) {
        stop(sprintf("The object '%s' in 'lpmodel' needs to be a list.",
                     mat.name),
             call. = FALSE)
      } else {
        err.ind <- c(err.ind, 3)
      }
    }

    # Category 4: Check if it is a function that produces two matrices
    # Category 5: If it is a list that contains two elements - one matrix
    # and one vector
    if ("function_obs_var_bs" %in% mat.cat |
        "function_obs_var" %in% mat.cat |
        "list_vector" %in% mat.cat) {
      if (inherits(mat, "function") | "list_vector" %in% mat.cat) {
        # Check whether the function can accept 'data' in the 'data.frame'
        # format
        if (inherits(mat, "function")) {
          check.datafunction(data, mat, mat.name)
          func.output <- mat(data)
        } else if (inherits(mat, "list")) {
          func.output <- mat
        }

        msg.vecmat <- sprintf(paste0("The output of '%s' in 'lpmodel' needs ",
                                     "to be a list of two objects (one vector ",
                                     "and one matrix)."),
                              mat.name)

        # Check if there are two outputs of the function
        if (length(func.output) != 2) {
          if (!("function_obs_var_bs" %in% mat.cat)) {
            stop(msg.vecmat, call. = FALSE)
          } else {
            return(list(sample = func.output))
          }
        } else {
          for (i in 1:2) {
            if (inherits(func.output[[i]], "numeric")) {
              func.output[[i]] <- matrix(func.output[[i]])
            }
          }
          out1 <- func.output[[1]]
          out2 <- func.output[[2]]

          # Check if one of out1 and out2 if a vector and the remaining one is
          # a matrix
          if ((ncol(out1) != 1 & ncol(out2) != 1)) {
            if (!("function_obs_var_bs" %in% mat.cat)) {
              stop(msg.vecmat, call. = FALSE)
            } else {
              return(list(sample = func.output))
            }
          } else if (ncol(out1) == 1) {
            sample <- out1
            mat <- out2
          } else if (ncol(out2) == 1) {
            sample <- out2
            mat <- out1
          }
        }
        return(list(mat = mat,
                    sample = sample))
      } else if (length(mat.cat) == 1) {
        if (!("function_obs_var_bs" %in% mat.cat)) {
          stop(msg.vecmat, call. = FALSE)
        } else {
          return(list(sample = func.output))
        }
      } else {
        err.ind <- c(err.ind, 4)
      }
    }

    # Display an overall error message
    err.msg.comb <- NULL
    if (length(err.ind) == length(mat.cat)) {
      # Pick the messages to display
      for (i in 1:length(mat.cat)) {
        if (i != 1) {
          err.msg.comb <- paste0(err.msg.comb, ", ")
        }
        if (mat.cat[i] == 1) {
          err.msg.comb <- paste0(err.msg.comb, "a matrix")
        } else if (mat.cat[i] == 2) {
          err.msg.comb <- paste0(err.msg.comb,
                                 "a function that returns ",
                                 "a matrix")
        } else if (mat.cat[i] == 3) {
          err.msg.comb <- paste0(err.msg.comb, "a list")
        } else if (mat.cat[i] == 3) {
          err.msg.comb <- paste0(err.msg.comb, paste0("a function that ",
                                                      "returns a matrix and ",
                                                      "a vector"))
        }
      }

      # Display the error message
      stop(sprintf(paste0("The object '%s' in 'lpmodel' has to be one of ",
                          "the followings: ", err.msg.comb),
                   mat.name),
           call. = FALSE)
    }
  } else {
    stop(sprintf("The component '%s' is required in the 'lpmodel' object.",
                 mat.name),
         call. = FALSE)
  }
}

#' Check function: matrix
#'
#' @importFrom methods is
#'
#' @description This function checks if the matrix objects in the
#'    \code{lpmodel} object are in the correct format.
#'
#' @inheritParams check.lpobjects
#' @param mat The matrix object in \code{lpmodel}.
#' @param mat.name The name of the matrix object in \code{lpmodel}.
#' @param inside.list An indicator variables of whether the object that is
#'    being checked is inside a list or not.
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.matrix <- function(mat, mat.name, mat.cat, inside.list) {
  # General message indicating the object has to be either a data.frame, matrix,
  # numeric, or sparseMatrix
  msg.matdf <- paste0("The objects inside the list object '%s' of ",
                      "'lpmodel' has to be one of the followings: data.frame,",
                      "matrix, numeric, or sparseMatrix.")

  if (inherits(mat, "data.frame") |
      inherits(mat, "matrix") |
      inherits(mat, "numeric")) {
    if (is.null(dim(mat))) {
      if (mat.name %in% c("A.obs", "A.shp", "A.tgt")) {
        mat.update <- matrix(mat, nrow = 1)
      } else if (mat.name %in% c("beta.obs", "beta.shp")) {
        mat.update <- matrix(mat, ncol = 1)
      }
    } else {
      mat.update <- as.matrix(mat)
    }
    return(list(mat.update = mat.update,
                err.ind = 0,
                dim = dim(mat.update)))
  } else if (isTRUE(methods::is(mat, "sparseMatrix"))) {
    return(list(mat.update = mat,
                err.ind = 0,
                dim = dim(mat)))
  } else if (length(mat.cat) == 1) {
    if (inside.list == FALSE) {
      stop(sprintf(paste0("The class of the object '%s' in 'lpmodel' has to ",
                          "be one of the followings: data.frame, matrix, ",
                          "numeric, or sparseMatrix."),
                   mat.name),
           call. = FALSE)
    } else {
      stop(sprintf(msg.matdf, mat.name), call. = FALSE)
    }
  } else {
    if (inside.list == FALSE) {
      mat.update <- NULL
      return(list(mat.update = mat.update,
                  err.ind = 1))
    } else if (inherits(mat, "list")) {
      mat.update <- NULL
      return(list(mat.update = mat.update,
                  err.ind = 1))
    } else {
      stop(sprintf(msg.matdf, mat.name), call. = FALSE)
    }
  }
}

#' Check function: vector
#'
#' @importFrom methods is
#'
#' @description This function checks if the matrix objects in \code{lpmodel}
#'    are in the correct format. If not, an error message is displayed.
#'
#' @param vec The vector object in \code{lpmodel}.
#' @param vec.name The name of the vector object in \code{lpmodel}.
#' @param inside.list An indicator variables of whether the object that is
#'    being checked is inside a list or not.
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.vector <- function(vec, vec.name, inside.list) {
  msg.vector <- paste0("The object '%s' in 'lpmodel' has to be a %s.")

  # Turn it into a matrix if it is not a list
  if (!is.list(vec) | !methods::is(vec, "sparseMatrix")) {
    vec <- as.matrix(vec)
  }
  if (nrow(vec) != 1 & ncol(vec) != 1) {
    if (inside.list == FALSE) {
      stop(sprintf(msg.vector, vec.name, "vector"), call. = FALSE)
    } else {
      stop(sprintf(msg.vector, vec.name, "list of vectors"), call. = FALSE)
    }
  }
}

#' Check function: check the number of cores
#'
#' @description This function checks if the number of cores specified is a
#'    positive integer. If not, set it as one.
#'
#' @param x Number of cores.
#'
#' @return Returns the updated number of cores.
#'   \item{x}{The updated number of cores}
#'
#' @export
#'
check.cores <- function(x) {
  if ((is.numeric(x) == TRUE &
       length(x) == 1 & x > 0 & x %% 1 == 0) == FALSE) {
    x <- 1
  }

  return(x)
}

#' Check function: check if \code{beta.tgt} is within the logical bound
#'
#' @description This function checks whether the parameter \code{beta.tgt}
#'   is within the logical bounds. If it is not within the logical bound,
#'   then reject immediately.
#'
#' @inheritParams dkqs
#' @param lpmodel An \code{lpmodel} object.
#' @param solver A linear or quadratic programming solver. The exact solver
#'   that is supported depends on the test chosen.
#'
#' @return Returns the following objects:
#'   \item{inout}{The variable that equals 1 if it is within the logical bound.
#'   This equals 0 otherwise.}
#'   \item{lb}{The logical lower bound.}
#'   \item{ub}{The logical upper bound.}
#'
#' @export
#'
check.betatgt <- function(data, lpmodel, beta.tgt, solver) {
  # ---------------- #
  # Step 1: Compute the logical upper and lower bound
  # ---------------- #
  ub <- check.betatgt.lp(data, lpmodel, "max", solver)
  lb <- check.betatgt.lp(data, lpmodel, "min", solver)

  # ---------------- #
  # Step 2: Return the indicator
  # ---------------- #
  if (beta.tgt <= ub & beta.tgt >= lb) {
    inout <- 1
  } else {
    inout <- 0
  }
  return(list(inout = inout,
              lb = lb,
              ub = ub))
}

#' Construct the linear program for the function
#'   \code{\link[lpinfer]{check.betatgt}}
#'
#' @description This function solves the linear program that is used to find
#'   the logical bounds in the \code{\link[lpinfer]{check.betatgt}} function.
#'
#' @inheritParams check.betatgt
#' @param modelsense A string that indicates whether the program is a
#'   maximization or minimization problem. If it is "max", then it is referring
#'   to a maximization problem. Otherwise, it is referring to a minimization
#'   problem.
#'
#' @return Returns the solution to the linear program.
#'   \item{objval}{The optimal value of the linear program. It is either the
#'   logical upper bound or the logical lower bound.}
#'
#' @details The linear program used to find the logical lower bound is
#' \deqn{
#'   \min_{x \in \mathbf{R}^d} \,\,
#'   A_{\rm tgt} x
#'   \quad \mathrm{s.t.} \quad
#'   A_{\rm shp} x = \beta_{\rm shp}
#'   \quad \mathrm{ and } \quad
#'   x \geq 0.
#' }
#' The linear program used to find the logical upper bound replaces the min
#' operator by the max operator.
#'
#' @export
#'
check.betatgt.lp <- function(data, lpmodel, modelsense, solver) {
  # ---------------- #
  # Step 1: Obtain the deterministic components of the objects in lpmodel
  # ---------------- #
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)

  # ---------------- #
  # Step 2: Solve the linear program
  # ---------------- #
  # Set up the arguments for optimizer
  optim.arg <- list(Af = NULL,
                    bf = A.tgt.hat,
                    nf = 1,
                    A = A.shp.hat,
                    rhs = beta.shp.hat,
                    sense = "=",
                    modelsense = modelsense,
                    lb = rep(0, ncol(A.shp.hat)))

  # Solve the linear program
  ans <- do.call(solver, optim.arg)

  if (ans$status == "INF_OR_UNBD") {
    ans$objval <- ifelse(modelsense == "max", Inf, -Inf)
  }
  # ---------------- #
  # Step 3: Return result
  # ---------------- #
  return(ans$objval)
}

#' General message for infeasible \code{beta.tgt}
#'
#' @description This function prints the message to inform the user that
#'   the \eqn{p}-value is directly set to 0 because they have specified a
#'   \code{beta.tgt} parameter that is infeasible, i.e. outside the logical
#'   bounds of the program.
#'
#' @return Returns a string of the message.
#'   \item{msg.explain}{The message to explain why the \eqn{p}-value is zero.}
#'   \item{msg.pval}{The message indicating that the \eqn{p}-value is zero.}
#'
#' @export
#'
infeasible.msg.betatgt <- function() {
  msg.explain <- paste0("Computation is skipped because the parameter ",
                        "'beta.tgt' is outside the logical bound.")
  msg.pval <- "p-value: 0"

  return(list(msg.explain = msg.explain,
              msg.pval = msg.pval))
}

#' Wrapper for \code{infeasible.msg.betatgt}
#'
#' @description This function is a wrapper for the
#'   \code{infeasible.msg.betatgt} function to print that \eqn{p}-value is
#'   zero in the \code{print} or \code{summary} messages.
#'
#' @return Nothing is returned.
#'
#' @export
#'
infeasible.pval.msg <- function() {
  msg <- infeasible.msg.betatgt()
  cat(msg$msg.pval)
}

#' Display warning message for infeasible \code{beta.tgt}
#'
#' @description This function displays the warning message for infeasible
#'   \code{beta.tgt}.
#'
#' @return Nothing is returned.
#'
#' @export
#'
infeasible.betatgt.warning <- function() {
  msg <- infeasible.msg.betatgt()$msg.explain
  warning(msg, call. = FALSE, immediate. = TRUE)
}

#' General error for checking the objects
#'
#' @description This is a general function to produce the error messages in
#'   the check functions.
#'
#' @param name.var The name of the object.
#' @param needs.to.be What the object is supposed to be.
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.errormsg <- function(name.var, needs.to.be) {
  general.msg <- "The object '%s' has to be %s."
  stop(sprintf(general.msg, name.var, needs.to.be), call. = FALSE)
}

#' Check function: sample size \code{n} if \code{data} is \code{NULL}
#'
#' @description This function checks the sample size \code{n} if \code{data}
#'   is \code{NULL}. This function is only used if \code{data} is \code{NULL}.
#'   Typically, when \code{data} is \code{NULL}, it refers to the case where
#'   the bootstrap replications are all passed to the \code{lpmodel} object.
#'
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.samplesize <- function(x, name.var) {
  # General message
  samplesize.msg <- sprintf(paste0("When 'data' is 'NULL', the object '%s' ",
                                   "has to be a positive integer."),
                           name.var)

  # Return error if x is not numeric or if it is not a positive integer
  if (!is.numeric(x)) {
    stop(samplesize.msg)
  } else if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0 & x %% 1 == 0)
             == FALSE) {
    # Call the general error message function
    stop(samplesize.msg)
  }
}

#' Check function: check brackets
#'
#' @description This function checks the initial brackets provided by the
#'   user for the \code{\link[lpinfer]{invertci}} function. It can either
#'   be \code{NULL} or a \code{numeric} of length 2. The initial brackets have
#'   to be finite.
#'
#' @param init.b An initial bracket.
#' @param name.var The name of the initial bracket.
#' @param bd A string that indicates whether the initial bracket for the lower
#'   bound (\code{lb}) or the upper bound (\code{ub}) is being evaluated.
#'
#' @return Returns the following intervals.
#'   \item{lb}{The lower bound of the initial bracket.}
#'   \item{ub}{The upper bound of the initial bracket.}
#'
#' @export
#'
check.initb <- function(init.b, name.var, bd) {
  # Check whether the length is 0, 1, 2 and whether it is numeric
  if (is.null(init.b)) {
    lb <- NULL
    ub <- NULL
  } else if ((length(init.b) == 2) & is.numeric(init.b)) {
    lb <- init.b[1]
    ub <- init.b[2]
    # Return an error if the upper bound is not larger than the lower bound
    if (lb >= ub) {
      stop(sprintf(paste0("When the initial bracket '%s' is a numeric of ",
                          "length 2, the first element has to be smaller ",
                          "then than second element."),
                   name.var))
    }

    # Return an error if the bounds are infinite
    if (is.infinite(lb) | is.infinite(ub)) {
      stop(sprintf("The initial bracket '%s' has to be finite. ", name.var))
    }
  } else {
    stop(sprintf(paste0("The initial bracket '%s' can either be ",
                        "NULL or a numeric of length 2."),
                 name.var))
  }

  return(list(lb = lb,
              ub = ub))
}
