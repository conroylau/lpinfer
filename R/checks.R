#' Check function: data frame
#' 
#' @description This function checks whether the class of the variable 
#'   is \code{data.frame} or \code{matrix}. If yes, the object will be 
#'   returned as an object in the \code{data.frame} class. If not, an
#'   error message is displayed.
#' 
#' @param data Data frame to be checked. 
#' @param name_var Name of the variable.
#' 
#' @return Returns the object in the format of \code{data.frame}.
#'   
#' @export
#' 
check_dataframe <- function(data, name_var){
  # Check data
  if (class(data) %in% c("data.frame", "matrix") == TRUE){
    data = as.data.frame(data)  
  } else {
    stop(sprintf(paste0("The data povided '%s' must either be a data.frame, ", 
    "a data.table, or a matrix."), name_var), call. = FALSE)    
  }
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
#' @inheritParams check_dataframe
#' 
#' @return Nothing is returned.
#'   
#' @export
#' 
check_positiveinteger <- function(x, name_var){
  if ((is.numeric(x) == TRUE & length(x) == 1 & x > 0 & x%%1 == 0) == FALSE){
    stop(sprintf("The variable '%s' has to be a positive integer.", x),
         call. = FALSE)
  }
}

#' Check function: numeric
#' 
#' @description This function checks whether the class of the variable 
#'   is \code{numeric} and has length 1.
#' 
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Nothing is returned.
#'   
#' @export
#' 
check_numeric <- function(x, name_var){
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
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Nothing is returned.
#'   
#' @export
#' 
check_boolean <- function(x, name_var){
  if (!(x == TRUE | x == FALSE)){
    stop(sprintf(paste0("The variable '%s' has to be boolean variable."), 
                 name_var), call. = FALSE)  
  }
}

#' Check function: number within a range
#' 
#' @description This function checks whether the variable is within a 
#'   certain interval. If not, an error message is displayed.
#' 
#' @param left_type Type of the left interval (\code{open} or \code{closed})
#' @param left Value of lower bound
#' @param right_type Type of the right interval (\code{open} or \code{closed})
#' @param right Value of the upper bound
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Nothing is returned.
#' 
#' @export
#' 
check_numrange <- function(x, name_var, left_type, left, right_type, right){
  # = = = = = = 
  # Step 1: Check if the number is numeric
  # = = = = = = 
  check_numeric(x, name_var)
  
  # = = = = = = 
  # Step 2: Check if the number is within the range
  # = = = = = = 
  if (right_type == "open" & right_type == "open"){
    if ((x > left & x < right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "(%s, %s)."), name_var, left, right), 
           call. = FALSE)  
    }
  } else if (right_type == "open" & right_type == "closed"){
    if ((x > left & x <= right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "(%s, %s]."), name_var, left, right), 
           call. = FALSE)  
    }   
  } else if (right_type == "closed" & right_type == "open"){
    if ((x >= left & x < right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "[%s, %s)."), name_var, left, right), 
           call. = FALSE)  
    }
  } else if (right_type == "closed" & right_type == "closed"){
    if ((x >= left & x <= right) == FALSE){
      stop(sprintf(paste0("The variable '%s' has to be inside the interval ",
                          "[%s, %s]."), name_var, left, right), 
           call. = FALSE)  
    }
  }
}

#' Check function: norm
#' 
#' @description This function checks whether the the norm used in the
#'   problem is an L1 or L2 norm. If not, an error message is displayed.
#' 
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Nothing is returned.
#' 
#' @export
#' 
check_norm <- function(x, name_var){
  if (x != 1 & x != 2){
    stop("Only 1-norm and 2-norm are supported in this function.", 
         call. = FALSE)
  }
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
#' @param name_A Name of matrix
#' @param mat_type Type of the matrix to be checked
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Returns the output of the function in the format of \code{numeric}
#'   and the correct dimension. If \code{mat_type} is \code{col}, then a 
#'   column vector is returned. If \code{mat_type} is \code{matrix}, then a
#'   square matrix is returned.
#' 
#' @export
#' 
check_func <- function(f, A, data, name_var, name_A, mat_type){
  # Check the class of the function
  if (class(f) != "function"){
    stop(sprintf("The input of '%s' has to be a function.", name_var), 
         call. = FALSE)
  } else{
    out = f(data)
    out = as.matrix(out)
    # Check if the output is numeric
    if (is.numeric(out[,1]) == FALSE){
      stop(sprintf("The output of '%s' has to be numeric.", name_var),
           call. = FALSE)
    } else{
      if (mat_type == "col"){
        # If the output has to be a column vector
        if (dim(out)[2] != 1){
          stop(sprintf("The output of '%s' has to be a column vector", 
                       name_var), call. = FALSE)
        } else if (dim(out)[1] != dim(A)[1]){
          stop(sprintf(paste0("The number of rows in the output of '%s' has ",
                              "to be the same as the number of rows in the ",
                              "matrix '%s'."), name_var, name_A), 
               call. = FALSE)
        } 
      } else if (mat_type == "square"){
        # If the output has to be a square matrix
        if (nrow(out) != ncol(out)){
          stop(sprintf("The output of '%s' has to be a square matrix", 
                       name_var), call. = FALSE)
        }
        if ((nrow(out) != nrow(A)) | (ncol(out) != nrow(A))){
          stop(sprintf("The number of rows and columns of the output ",
                       "for '%s' has to be equal to the number of rows in ", 
                       "matrix '%s'.", name_var, name_A),
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
check_Ab <- function(A, b, Aname, bname){
  if (is.null(A) + is.null(b) == 1){
    # = = = = = = 
    # Step 1: Check that if A and b must be both NULL or both non-NULL
    # = = = = = = 
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
    # = = = = = = 
    # Step 2: Both A and b are NULL. Nothing else to do
    # = = = = = = 
    #### Step 2: Both A and b are NULL. Nothing else to do
  } else {
    # = = = = = = 
    # Step 3: Checks for the case where both A and b are non-NULL
    # = = = = = = 
    matrix_names = c(Aname, bname)
    matrix_list = list(A, b)
    for (i in 1:2){
      ## Part 1: Check the format of the matrices
      if (class(matrix_list[[i]]) %in% 
          c("data.frame", "matrix", "numeric") == FALSE){
        stop(gsub("\\s+", " ",
                  paste0("The argument '", matrix_names[i], "' must either be 
                a data.frame, data.table, or matrix.")), call. = FALSE)   
      } else{
        # Ensure the variable is in matrix form
        matrix_list[[i]] = as.matrix(matrix_list[[i]])
        
        ## Part 2: Check whether the matrices are numeric
        if (is.numeric(matrix_list[[i]]) == FALSE){
          stop(paste0("The argument '", matrix_names[i], 
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
    
    # = = = = = = 
    # Step 4: Update class
    # = = = = = = 
    # Ensure that both A and b to ensure that they are both matrices
    A = matrix_list[[1]]
    b = matrix_list[[2]]
  }
  
  # = = = = = = 
  # Step 5: Return results
  # = = = = = = 
  return(list(A = A,
              b = b))
}

#' Check function: check solvers
#' 
#' @description This function checks the solver used is supported and is
#'   appropriate for the problem considered.
#' 
#' @param norm Norm used in the problem
#' @inheritParams check_dataframe
#' @inheritParams check_positiveinteger
#' 
#' @return Returns the function name that corresponds to the solver for the 
#'   problem.
#'    \item{solver}{Function for the solver used.}
#'    \item{solver_name}{Name of the solver used in lower case.}
#'      
#' @export
#' 
check_solver <- function(x, name_var, norm = 2){
  # = = = = = = 
  # Step 1: Preparation and hard-coded information
  # = = = = = = 
  x = tolower(x)
  # Package recommendation messages
  gurobi_msg = "'gurobi' (version 8.1-1 or later)"
  cplexapi_msg = "'cplexAPI' (version 1.3.3 or later)"
  rcplex_msg = "'Rcplex' (version 0.3-3 or later)"
  limsolve_msg = "'limSolve' (version 1.5.6 or later)"
  lpsolveapi_msg = "lpSolveAPI (version 5.5.2.0 or later)"
  
  # = = = = = = 
  # Step 2a: If no solver name is provided by the user
  # = = = = = = 
  if (is.null(x) == TRUE){
    # If gurobi is installed, the gurobi solver will be used for L1- & L2-norm
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
      solver = gurobi_optim
    } else if (norm == 1) {
      # If L1-norm is used, other solvers will be checked
      if (requireNamespace("limSolve", quietly = TRUE) == TRUE){
        solver = limsolve_optim
      } else if (requireNamespace("Rcplex", quietly = TRUE) == TRUE){
        solver = rcplex_optim
      } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE){
        solver = cplexapi_optim
      } else if (requireNamespace("lpsolveAPI", quietly = TRUE) == TRUE){
        solver = lpsolveapi_optim
      }      
    } else {
      if (norm == 1){
        stop(gsub("\\s+", " ",
                  paste0("This function is incompatible with '", x, 
                         "' when L1-norm is chosen in the estimation procedure. 
                         Please install one of the following packages to solve 
                         the linear and quadratic programs: ",
                         gurobi_msg, "; ",
                         cplexapi_msg, "; ",
                         rcplex_msg, "; ",
                         limsolve_msg, ";",
                         lpsolveapi_msg, ".")),
             call. = FALSE)
      } else if (norm == 2){
        stop(gsub("\\s+", " ",
                  paste0("This function with L2-norm in the estimation 
                         procedure is only incompatible with 'gurobi'. ", 
                         "Please install ", gurobi_msg, " to obtain the
                       bounds of the problem subject to shape restriction.")), 
             call. = FALSE)
      }
    }
    if (progress == TRUE){
      cat(paste("No solver solver is suggested by the user. The solver", 
                x, "is chosen.\n", sep = ""))
    }
  } else if (x == "gurobi"){
    # = = = = = = 
    # Step 2b: If a solver name is provided by the user
    # = = = = = = 
    ## Case 1: If the user specified the solver as 'gurobi'
    solver = gurobi_optim
  } else if (x == "limsolve"){
    ## Case 2: If the user specified the solver as 'limSolve'
    solver = limsolve_optim
  } else if (x == "rcplex"){
    ## Case 3: If the user specified the solver as 'rcplex'
    solver = rcplex_optim
  } else if (x == "cplexapi"){
    ## Case 4: If the user specified the solver as 'cplexapi'
    solver = cplexapi_optim
  } else if (x == "lpsolveapi" & norm == 1){
    ## Case 5: If the user specified the solver as 'lpsolveapi'
    solver = lpsolveapi_optim
  } else {
    ## Case 6: If the user specified a solver that is not compatible
    if (norm == 1){
      stop(gsub("\\s+", " ",
                paste0("This function is incompatible with '", x, 
                       "' when L1-norm is chosen in the estimation procedure. 
                       Please install one of the following packages to solve 
                       the linear and quadratic programs: ",
                       gurobi_msg, "; ",
                       cplexapi_msg, "; ",
                       rcplex_msg, "; ",
                       limsolve_msg, ".")),
           call. = FALSE)
    } else if (norm == 2){
      stop(gsub("\\s+", " ",
                paste0("This function with L2-norm in the estimation procedure
                       is only incompatible with 'gurobi'. ", 
                       "Please install ", gurobi_msg, " to obtain the
                       bounds of the problem subject to shape restriction.")), 
           call. = FALSE)
    }
  }
  
  # = = = = = = 
  # Step 3: Returns the function name that corresponds to the solver
  # = = = = = = 
  return(list(solver = solver,
              solver_name = x))
}