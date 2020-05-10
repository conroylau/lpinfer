#' Check function: data frame
#'
#' @description This function force the data input as a \code{matrix}.
#'
#' @param data Data frame to be checked.
#'
#' @return Returns the object as a \code{matrix}.
#'
#' @export
#'
check.dataframe <- function(data){
  # Check data
  data <- as.matrix(data)

  # Return updated data
  return(data)
}

#' Check function: positive integer
#'
#' @description This function checks whether the class of the variable
#'   is \code{numeric}, has length 1 and is a positive integer. If not, an
#'   error message is displayed.
#'
#' @param x Variable to be checked.
#' @inheritParams check.dataframe
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.positiveinteger <- function(x, name.var){
  if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0 & x%%1 == 0) == FALSE){
    stop(sprintf("The variable '%s' has to be a positive integer.", x),
         call. = FALSE)
  }
  return(x)
}

#' Check function: numeric
#'
#' @description This function checks whether the class of the variable
#'   is \code{numeric} and has length 1.
#'
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.numeric <- function(x, name.var){
  if ((is.numeric(x) == TRUE & length(x) == 1) == FALSE){
    stop(sprintf("The class of the variable '%s' has to be numeric.", x),
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
check.boolean <- function(x, name.var){
  if (!(x == TRUE | x == FALSE)){
    stop(sprintf(paste0("The variable '%s' has to be boolean variable."),
                 name.var), call. = FALSE)
  }
}

#' Check function: number within a range
#'
#' @description This function checks whether the variable is within a
#'   certain interval. If not, an error message is displayed.
#'
#' @param left.type Type of the left interval (\code{open} or \code{closed})
#' @param left Value of lower bound
#' @param right.type Type of the right interval (\code{open} or \code{closed})
#' @param right Value of the upper bound
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.numrange <- function(x, name.var, left.type, left, right.type, right){
  # ---------------- #
  # Step 1: Check if the number is numeric
  # ---------------- #
  check.numeric(x, name.var)

  # ---------------- #
  # Step 2: Check if the number is within the range
  # ---------------- #
  if (left.type == "open" & right.type == "open"){
    if ((x > left & x < right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "(%s, %s)."), name.var, left, right),
           call. = FALSE)
    }
  } else if (left.type == "open" & right.type == "closed"){
    if ((x > left & x <= right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "(%s, %s]."), name.var, left, right),
           call. = FALSE)
    }
  } else if (left.type == "closed" & right.type == "open"){
    if ((x >= left & x < right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "[%s, %s)."), name.var, left, right),
           call. = FALSE)
    }
  } else if (left.type == "closed" & right.type == "closed"){
    if ((x >= left & x <= right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "[%s, %s]."), name.var, left, right),
           call. = FALSE)
    }
  }
}

#' Check function: norm
#'
#' @description This function checks whether the the norm used in the
#'   problem is an L1 or L2 norm. If not, an error message is displayed.
#'
#' @details The following input for \code{norm} will be interpreted as the L1-
#'   norm:
#'   \itemize{
#'     \item{\code{1} (numeric)}
#'     \item{\code{"1"} (string)}
#'     \item{\code{"L1"}}
#'     \item{\code{"one"}}
#'     \item{\code{"o"}}
#'     \item{\code{"taxicab"}}
#'   }
#'   The following input for \code{norm} will be interpreted as the L2-norm:
#'   \itemize{
#'     \item{\code{2} (numeric)}
#'     \item{\code{"2"} (string)}
#'     \item{\code{"L2"}}
#'     \item{\code{"two"}}
#'     \item{\code{"t"}}
#'     \item{\code{"e"}}
#'     \item{\code{"euclidean"}}
#'   }
#'   Note that capitalization is not an issue here as the text will be brought
#'   to the lower case.
#'
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Nothing is returned.
#'
#' @export
#'
check.norm <- function(x, name.var){
  # Bring the variable to lower case
  x <- tolower(x)

  if (x %in% c(1, "1", "l1", "one", "o", "taxicab")){
    # Case 1: user provied an input that corresponds to L1-norm
    norm <- 1
  } else if (x %in% c(2, "2", "l2", "two", "t", "e", "euclidean")){
    # Case 1: user provied an input that corresponds to L2-norm
    norm <- 2
  } else {
    stop(gsub("\\s+", " ",
              paste0("Only 1-norm and 2-norm are supported in this function. ",
                     "For 1-norm, please enter one of the followings: 1, '1', ",
                     "'l1', 'one', 'o' or 'taxicab'. For 2-norm, please enter ",
                     "one of the followings: 2, '2', 'l2', 'two', 't', 'e', ",
                     "'euclidean'.")),
         call. = FALSE)
  }

  # Return the updated norm
  return(norm)
}

#' Check function: function
#'
#' @description This function checks the class of the function provided by
#'   the user is \code{function}, and whether the output produced by the
#'   function is \code{numeric} and has the correct dimension. If not, an
#'   error message is displayed.
#'
#' @param f Function to be tested
#' @param A Matrix that will be contains dimensional information for the
#'    output of f
#' @param data Data frame that will be passed to the function
#' @param name.A Name of matrix
#' @param mat.type Type of the matrix to be checked
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Returns the output of the function in the format of \code{numeric}
#'   and the correct dimension. If \code{mat.type} is \code{col}, then a
#'   column vector is returned. If \code{mat.type} is \code{matrix}, then a
#'   square matrix is returned.
#'
#' @export
#'
check.func <- function(f, A, data, name.var, name.A, mat.type){
  # Check the class of the function
  if (class(f) != "function"){
    stop(sprintf("The input of '%s' has to be a function.", name.var),
         call. = FALSE)
  } else{
    out <- f(data)
    out <- as.matrix(out)
    # Check if the output is numeric
    if (is.numeric(out[,1]) == FALSE){
      stop(sprintf("The output of '%s' has to be numeric.", name.var),
           call. = FALSE)
    } else{
      if (mat.type == "col"){
        # If the output has to be a column vector
        if (dim(out)[2] != 1){
          stop(sprintf("The output of '%s' has to be a column vector",
                       name.var), call. = FALSE)
        } else if (dim(out)[1] != dim(A)[1]){
          stop(sprintf(paste0("The number of rows in the output of '%s' has ",
                              "to be the same as the number of rows in the ",
                              "matrix '%s'."), name.var, name.A),
               call. = FALSE)
        }
      } else if (mat.type == "square"){
        # If the output has to be a square matrix
        if (nrow(out) != ncol(out)){
          stop(sprintf("The output of '%s' has to be a square matrix",
                       name.var), call. = FALSE)
        }
        if ((nrow(out) != nrow(A)) | (ncol(out) != nrow(A))){
          stop(sprintf("The number of rows and columns of the output ",
                       "for '%s' has to be equal to the number of rows in ",
                       "matrix '%s'.", name.var, name.A),
               call. = FALSE)
        }
      }
    }
  }

  return(out)
}

#' Check function: Constraint matrix and the corresponding RHS vector
#'
#' @description This function checks the constraint matrix \eqn{\bm{A}} and
#'    the corresponding vector \eqn{\bm{\beta}}.
#'
#' @param A Constraint matrix.
#' @param b Corresponding vector.
#' @param Aname Variable name of constraint matrix.
#' @param bname Varaible name of the corresponding vector.
#'
#' @return Returns the updated matrix and vector or print stop message.
#'    \item{A_updated}{Updated constraint matrix.}
#'    \item{b_updated}{Updated rhs vector.}
#'
#' @export
#'
check.Ab <- function(A, b, Aname, bname){
  if (is.null(A) + is.null(b) == 1){
    # ---------------- #
    # Step 1: Check that if A and b must be both NULL or both non-NULL
    # ---------------- #
    if (is.null(A) == TRUE){
      stop(sprintf("'%s' is NULL but '%s' is not NULL. Please ensure
                   that either they are both NULL or they are both
                   valid matrices.", Aname, bname))
    } else {
      stop(sprintf("'%s' is NULL but '%s' is not NULL. Please ensure
                   that either they are both NULL or they are both
                   valid matrices.", bname, Aname))
    }
  } else if (is.null(A) + is.null(b) == 2){
    # ---------------- #
    # Step 2: Both A and b are NULL. Nothing else to do
    # ---------------- #
  } else {
    # ---------------- #
    # Step 3: Checks for the case where both A and b are non-NULL
    # ---------------- #
    matrix.names <- c(Aname, bname)
    matrix.list <- list(A, b)
    for (i in 1:2){
      ## Part 1: Check the format of the matrices
      if (class(matrix.list[[i]]) %in%
          c("data.frame", "matrix", "numeric") == FALSE){
        stop(gsub("\\s+", " ",
                  paste0("The argument '", matrix.names[i], "' must either be
                a data.frame, data.table, or matrix.")), call. = FALSE)
      } else{
        # Ensure the variable is in matrix form
        matrix.list[[i]] <- as.matrix(matrix.list[[i]])

        ## Part 2: Check whether the matrices are numeric
        if (is.numeric(matrix.list[[i]]) == FALSE){
          stop(paste0("The argument '", matrix.names[i],
                      "' has to be numeric."), call. = FALSE)
        }
      }
    }

    ## Part 3: Ensure that beta is a column vector. This also captures the
    ## case of b being a scalar because it is turned into a matrix in the
    ## condition.
    if (dim(as.matrix(b))[2] != 1){
      stop(sprintf("The argument '%s' has to be a column vector", bname),
           call. = FALSE)
    }

    ## Part 4: Ensure that the number of rows of A and b are identical
    if (nrow(A) != length(b)){
      stop(sprintf("The number of rows of '%s' has to be equal to the number
                  of rows of '%s.", Aname, bname), call. = FALSE)
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

#' Check function: check solvers
#'
#' @description This function checks the solver used is supported and is
#'   appropriate for the problem considered.
#'
#' @param norm Norm used in the problem
#' @inheritParams check.dataframe
#' @inheritParams check.positiveinteger
#'
#' @return Returns the function name that corresponds to the solver for the
#'   problem.
#'    \item{solver}{Function for the solver used.}
#'    \item{solver_name}{Name of the solver used in lower case.}
#'
#' @export
#'
check.solver <- function(x, name.var, norm = 2){
  # ---------------- #
  # Step 1: Preparation and hard-coded information
  # ---------------- #
  x <- tolower(x)
  # Package recommendation messages
  gurobi.msg <- "'gurobi' (version 8.1-1 or later)"
  cplexapi.msg <- "'cplexAPI' (version 1.3.3 or later)"
  rcplex.msg <- "'Rcplex' (version 0.3-3 or later)"
  limsolve.msg <- "'limSolve' (version 1.5.6 or later)"
  lpsolveapi.msg <- "lpSolveAPI (version 5.5.2.0 or later)"

  # ---------------- #
  # Step 2a: If no solver name is provided by the user
  # ---------------- #
  if (is.null(x) == TRUE){
    # If gurobi is installed, the gurobi solver will be used for L1- & L2-norm
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
      solver = gurobi.optim
    } else if (norm == 1) {
      # If L1-norm is used, other solvers will be checked
      if (requireNamespace("limSolve", quietly = TRUE) == TRUE){
        solver <- limsolve.optim
      } else if (requireNamespace("Rcplex", quietly = TRUE) == TRUE){
        solver <- rcplex.optim
      } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE){
        solver <- cplexapi.optim
      } else if (requireNamespace("lpsolveAPI", quietly = TRUE) == TRUE){
        solver <- lpsolveapi.optim
      }
    } else {
      if (norm == 1){
        stop(gsub("\\s+", " ",
                  paste0("This function is incompatible with '", x,
                         "' when L1-norm is chosen in the estimation procedure.
                         Please install one of the following packages to solve
                         the linear and quadratic programs: ",
                         gurobi.msg, "; ",
                         cplexapi.msg, "; ",
                         rcplex.msg, "; ",
                         limsolve.msg, ";",
                         lpsolveapi.msg, ".")),
             call. = FALSE)
      } else if (norm == 2){
        stop(gsub("\\s+", " ",
                  paste0("This function with L2-norm in the estimation
                         procedure is only incompatible with 'gurobi'. ",
                         "Please install ", gurobi.msg, " to obtain the
                       bounds of the problem subject to shape restriction.")),
             call. = FALSE)
      }
    }
    if (progress == TRUE){
      cat(paste("No solver solver is suggested by the user. The solver",
                x, "is chosen.\n", sep = ""))
    }
  } else if (x == "gurobi"){
    # ---------------- #
    # Step 2b: If a solver name is provided by the user
    # ---------------- #
    ## Case 1: If the user specified the solver as 'gurobi'
    solver <- gurobi.optim
  } else if (x == "limsolve"){
    ## Case 2: If the user specified the solver as 'limSolve'
    solver <- limsolve.optim
  } else if (x == "rcplex"){
    ## Case 3: If the user specified the solver as 'rcplex'
    solver <- rcplex.optim
  } else if (x == "cplexapi"){
    ## Case 4: If the user specified the solver as 'cplexapi'
    solver <- cplexapi.optim
  } else if (x == "lpsolveapi" & norm == 1){
    ## Case 5: If the user specified the solver as 'lpsolveapi'
    solver <- lpsolveapi.optim
  } else {
    ## Case 6: If the user specified a solver that is not compatible
    if (norm == 1){
      stop(gsub("\\s+", " ",
                paste0("This function is incompatible with '", x,
                       "' when L1-norm is chosen in the estimation procedure.
                       Please install one of the following packages to solve
                       the linear and quadratic programs: ",
                       gurobi.msg, "; ",
                       cplexapi.msg, "; ",
                       rcplex.msg, "; ",
                       limsolve.msg, ".")),
           call. = FALSE)
    } else if (norm == 2){
      stop(gsub("\\s+", " ",
                paste0("This function with L2-norm in the estimation procedure
                       is only incompatible with 'gurobi'. ",
                       "Please install ", gurobi.msg, " to obtain the
                       bounds of the problem subject to shape restriction.")),
           call. = FALSE)
    }
  }

  # ---------------- #
  # Step 3: Returns the function name that corresponds to the solver
  # ---------------- #
  return(list(solver = solver,
              solver.name = x))
}

#' Check function: check lpmodel
#'
#' @description This function checks if the object \code{lpmodel} is in the
#'    correct format.
#'
#' @param lpmodel An \code{lpmodel} object.
#' @param name.var Name of the \code{lpmodel} object.
#' @param A.tgt.cat Category of the \code{A.tgt} object.
#' @param A.obs.cat Category of the \code{A.obs} object.
#' @param A.shp.cat Category of the \code{A.shp} object.
#' @param beta.obs.cat Category of the \code{beta.obs} object.
#' @param beta.shp.cat Category of the \code{beta.shp} object.
#' @inheritParams dkqs
#'
#' @details In each of the testing procedures, there are five possible
#' categories of the each of the object:
#'   \itemize{
#'     \item{Category '\code{not_used}': This refers to the case where the 
#'       object is not used in the function}
#'     \item{Category '\code{matrix}': This refers to the case where the
#'       object has to be a matrix}
#'     \item{Category '\code{function_mat}': This refers to the case where 
#'       the object has to be a function that produces a matrix.}
#'     \item{Category '\code{list}': This refers to the case where the
#'       object is a list}
#'     \item{Category '\code{function_obs_var}': This refers to the case 
#'       where the object is a function that produces a list that contains 
#'       a vector and a matrix. This is typically the case for \code{beta.obs}
#'       that the function produces an estimator and the asymptotic variance.}
#'   }
#'   Each object can belong to one of more categories.
#'
#' @return Returns the updated \code{lpmodel} object.
#'    \item{lpmodel}{Updated \code{lpmodel} object.}
#'
#' @export
#'
check.lpmodel <- function(data, lpmodel, name.var, A.tgt.cat, A.obs.cat,
                          A.shp.cat, beta.obs.cat, beta.shp.cat, R){
  # ---------------- #
  # Step 1: Check if lpmodel is a list
  # ---------------- #
  if (class(lpmodel) != "lpmodel"){
    stop(sprintf("The object '%s' has to be a list.", name.var),
         call. = FALSE)
  }

  # ---------------- #
  # Step 2: Check each of the objects (treat them as matrices)
  # and check if beta.obs and beta.shp are matrices
  # ---------------- #
  if (!("not_used" %in% A.tgt.cat)){
    A.tgt.return <- check.lpobjects(data, lpmodel$A.tgt, "A.tgt", A.tgt.cat, R)
  }
  if (!("not_used" %in% A.obs.cat)){
    A.obs.return <- check.lpobjects(data, lpmodel$A.obs, "A.obs", A.obs.cat, R)
  }
  if (!("not_used" %in% A.shp.cat)){
    A.shp.return <- check.lpobjects(data, lpmodel$A.shp, "A.shp", A.shp.cat, R)
  }
  if (!("not_used" %in% beta.obs.cat)){
    beta.obs.return <- check.lpobjects(data, lpmodel$beta.obs, "beta.obs",
                                       beta.obs.cat, R)
    check.vector(beta.obs.return$sample, "beta.obs", FALSE)
  }
  if (!("not_used" %in% beta.shp.cat)){
    beta.shp.return <- check.lpobjects(data, lpmodel$beta.shp, "beta.shp",
                                       beta.shp.cat, R)
    check.vector(beta.shp.return$sample, "beta.shp", FALSE)
  }

  # ---------------- #
  # Step 3: Check whether the dimension matches
  # ---------------- #
  if (!("not_used" %in% A.obs.cat)){
    if (nrow(A.obs.return$sample) != nrow(beta.obs.return$sample)){
      stop(paste0("The objects 'A.obs' and 'beta.obs' in 'lpmodel' has ",
                  "to have the same number of rows."))
    }
  }
  if (!("not_used" %in% A.shp.cat)){
    if (nrow(A.shp.return$sample) != nrow(beta.shp.return$sample)){
      stop(paste0("The objects 'A.shp' and 'beta.shp' in 'lpmodel' has ",
                  "to have the same number of rows."))
    }
  }
  return(lpmodel)
}

#' Check function: check matrices and vectors in lpmodel
#'
#' @description This function checks if the matrix objects in \code{lpmodel}
#'    are in the correct format.
#'
#' @param mat The matrix object in \code{lpmodel}.
#' @param mat.name Name of the matrix object in \code{lpmodel}.
#' @param mat.cat Category of the matrix object.
#' @inheritParams dkqs
#' 
#' @details See the details for the function `\code{check.lpmodel}` for the 
#'   details on the strings for each category.
#'
#' @return Returns two objects:
#'    \item{mat}{The updated object in \code{lpmodel}.}
#'    \item{sample}{A sample object that is used to compare the dimension.}
#'
#' @export
#'
check.lpobjects <- function(data, mat, mat.name, mat.cat, R){
  # Ignore if mat is null
  if (!is.null(mat)){
    # ---------------- #
    # Step 1: Check error indicator
    # ---------------- #
    err.ind <- NULL

    # Category 1: Check if it is a single matrix
    if ("matrix" %in% mat.cat) {
      mat.return <- check.matrix(mat, mat.name, mat.cat, FALSE)
      if (mat.return$err.ind != 1){
        return(list(mat = mat.return$mat.update,
                    sample = mat.return$mat.update))
      } else {
        err.ind <- c(err.ind, 1)
      }
    }

    # Category 2: Check if it is a function
    if ("function_mat" %in% mat.cat) {
      if (class(mat) == "function"){
        return(list(mat = mat,
                    sample = mat(data)))
      } else if (length(mat.cat) == 1){
        stop(sprintf(paste0("The object '%s' in 'lpmodel' has to be a ",
                            "function."),
                     mat.name),
             call. = FALSE)
      } else {
        err.ind <- c(err.ind, 2)
      }
    }

    # Category 3: Check if it is a list
    if ("list" %in% mat.cat) {
      if (class(mat) == "list"){
        # Check if the length of the list is R+1
        if (length(mat) != (R+1)){
          stop(sprintf(paste0("The object '%s' in 'lpmodel' needs to have ",
                              "exactly %s elements",
                              mat.name, R+1)),
               call. = FALSE)
        }

        # Loop through each element of the list to check if they are matrices
        df.dim <- data.frame(matrix(vector(), nrow = R+1))
        for (i in 1:R){
          mat.return <- check.matrix(mat[[i+1]], mat.name, mat.cat, TRUE)
          # Check if all objects inside the list has the same dimension
          df.dim[i,1] <- mat.return$dim[1]
          df.dim[i,2] <- mat.return$dim[2]
          if (i > 1){
            if ((df.dim[i,1] != df.dim[i-1,1]) |
                (df.dim[i,1] != df.dim[i-1,1])){
              stop(sprintf(paste0("The dimension of the objects inside the list ",
                                  "'%s' in 'lpmodel' need to have the same",
                                  "dimension.",
                                  mat.name)),
                   call. = FALSE)
            }
          }
        }
        return(list(mat = mat,
                    sample = mat.return[[1]]))
      } else if (length(mat.cat) == 1){
        stop(sprintf("The object '%s' in 'lpmodel' has to be a list.",
                     mat.name),
             call. = FALSE)
      } else {
        err.ind <- c(err.ind, 3)
      }
    }

    # Category 4: Check if it is a function that produces two matrices
    if ("function_obs_var" %in% mat.cat) {
      if (class(mat) == "function"){
        # Check if there are two outputs of the function
        func.output <- mat(data)
        if (length(func.output) != 2){
          stop(sprintf(paste0("The output of '%s' in 'lpmodel' has to be a ",
                              "list of two objects (one vector and one ",
                              "matrix)."),
                       mat.name),
               call. = FALSE)
        } else {
          out1 <- matrix(func.output[[1]])
          out2 <- matrix(func.output[[2]])

          # Check if one of out1 and out2 if a vector and the remaining one is
          # a matrix
          if ((ncol(out1) != 1 & ncol(out2) != 1)){
            stop(sprintf(paste0("The output of '%s' in 'lpmodel' has to be a ",
                                "list of two objects (one vector and one ",
                                "matrix)."),
                         mat.name),
                 call. = FALSE)
          } else if (ncol(out1) == 1){
            sample <- out1
          } else if (ncol(out1) == 2){
            sample <- out2
          }
        }
        return(list(mat = mat,
                    sample = sample))
      } else if (length(mat.cat) == 1){
        stop(sprintf(paste0("The object '%s' in 'lpmodel' has to be a ",
                            "function that produces a list of two objects ",
                            "(one vector and one matrix)."),
                     mat.name),
             call. = FALSE)
      } else {
        err.ind <- c(err.ind, 4)
      }
    }

    # Display an overall error message
    err.msg.comb <- NULL
    if (err.ind == length(mat.cat)){
      # Pick the messages to display
      for (i in 1:length(mat.cat)){
        if (i != 1){
          err.msg.comb <- paste0(err.msg.comb, ", ")
        }
        if (mat.cat[i] == 1){
          err.msg.comb <- paste0(err.msg.comb, "a matrix")
        } else if (mat.cat[i] == 2){
          err.msg.comb <- paste0(err.msg.comb,
                                 "a function that returns ",
                                 "a matrix")
        } else if (mat.cat[i] == 3){
          err.msg.comb <- paste0(err.msg.comb, "a list")
        } else if (mat.cat[i] == 3){
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
  }
}

#' Check function: check matrix
#'
#' @description This function checks if the matrix objects in \code{lpmodel}
#'    are in the correct format.
#'
#' @param mat The matrix object in \code{lpmodel}.
#' @param mat.name Name of the matrix object in \code{lpmodel}.
#' @param inside.list Indicator variables of whether the object that is
#'    being checked is inside a list or not.
#'
#' @export
#'
check.matrix <- function(mat, mat.name, mat.cat, inside.list){
  if (class(mat) == "data.frame" | class(mat) == "matrix" |
      class(mat) == "numeric"){
    mat.update <- as.matrix(mat)
    return(list(mat.update = mat.update,
                err.ind = 0,
                dim = dim(mat.update)))
  } else if (length(mat.cat) == 1){
    if (inside.list == FALSE){
      stop(sprintf(paste0("The class of teh object '%s' in 'lpmodel' has to ",
                          "be one of the followings: data.frame, matrix, or ",
                          "numeric."),
                   mat.name),
           call. = FALSE)
    } else {
      stop(sprintf(paste0("The objects inside the list object '%s' of ",
                          "'lpmodel' has to be either a matrix or a ",
                          "data.frame."),
                   mat.name),
           call. = FALSE)
    }
  } else {
    if (inside.list == FALSE) {
      mat.update <- NULL
      return(list(mat.update = mat.update,
                  err.ind = 1))
    } else if (class(mat) == "list"){
      mat.update <- NULL
      return(list(mat.update = mat.update,
                  err.ind = 1))
    } else {
      stop(sprintf(paste0("The objects inside the list object '%s' of ",
                          "'lpmodel' has to be either a matrix or a ",
                          "data.frame."),
                   mat.name),
           call. = FALSE)
    }
  }
}

#' Check function: check vector
#'
#' @description This function checks if the matrix objects in \code{lpmodel}
#'    are in the correct format.
#'
#' @param vec The vector object in \code{lpmodel}.
#' @param vec.name Name of the vector object in \code{lpmodel}.
#' @param inside.list Indicator variables of whether the object that is
#'    being checked is inside a list or not.
#'
#' @export
#'
check.vector <- function(vec, vec.name, inside.list){
  vec <- as.matrix(vec)

  if (nrow(vec) != 1 & ncol(vec) != 1){
    if (inside.list == FALSE){
      stop(sprintf(paste0("The object '%s' of 'lpmodel' has to be a vector."),
                   vec.name),
           call. = FALSE)
    } else {
      stop(sprintf(paste0("The object '%s' of 'lpmodel' has to be a list ",
                          "of vectors"),
                   vec.name),
           call. = FALSE)
    }
  }
}

#' Check function: check the number of cores
#'
#' @description This function checks if the number of cores specified is a
#'    positive integer. If not, set it as one.
#'
#' @param cores Number of cores.
#'
#' @return Returns the updated number of cores.
#'   \item{cores}{Updated number of cores}
#'
#' @export
#'
check.cores <- function(x){
  if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0 & x%%1 == 0) == FALSE){
    cores <- 1
  }
  return(cores)
}
