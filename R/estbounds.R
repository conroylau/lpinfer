#' Estimate bounds with shape restrictions
#' 
#' @description This function computes the bound of the linear program 
#'    subject to shape constraints. This function also offers an option
#'    to estimate the shape constraints using a two-step procedure and 
#'    some tolerance level.
#' 
#' @import Matrix gurobi
#' 
#' @param beta_obs Observed value of parameter of interest.
#' @param A_shp_eq Matrix representing equality shape constraints.
#' @param A_shp_ineq Matrix representing inequality shape constraints.
#' @param beta_shp_eq RHS vector in equality shape constraints.
#' @param beta_shp_ineq RHS vector in inequality shape constraints.
#' @param lnorm Norm used in the optimization problem.
#' @param kappa Parameter used in the second step of the two-step procedure 
#'    for obtaining the solution subject to the shape constraints.
#' @param estimate Boolean variable to indicate whether the estimated 
#'    problem should be considered.
#' @inheritParams dkqs
#' @inheritParams invertci
#' 
#' @return Returns the bounds subject to the shape constraints.
#'   \item{ub}{Upper bound with shape constraints}
#'   \item{lb}{Lower bound with shape constraints}
#'   \item{est}{Indicator of whether estimation is involved in the 
#'   estimation}
#'   \item{call}{The function that has been called.}
#'   \item{lnorm}{Norm used in the optimization problem.}
#'   
#'    
#' @export
#' 

estbounds <- function(df, func_obs, A_obs, A_tgt, beta_tgt,
                      A_shp_eq, A_shp_ineq, beta_shp_eq, beta_shp_ineq,
                      kappa = 1e-5, lnorm = 2, solver = "gurobi", 
                      estimate = TRUE, progress = TRUE){
  
  #### Step 1: Obtain call, check and update input
  # Obtain call information
  call = match.call()
  # Check and update
  estbounds_return = estbounds_check(df, func_obs, A_obs, A_tgt, beta_tgt,
                                     A_shp_eq, A_shp_ineq, beta_shp_eq, 
                                     beta_shp_ineq, kappa, lnorm, solver, 
                                     estimate, progress)
  # Update the input
  df = estbounds_return$df
  A_obs = estbounds_return$A_obs
  A_tgt = estbounds_return$A_tgt
  beta_obs = estbounds_return$beta_obs
  beta_tgt = estbounds_return$beta_tgt
  A_shp_eq = estbounds_return$A_shp_eq
  A_shp_ineq = estbounds_return$A_shp_ineq
  beta_shp_eq = estbounds_return$beta_shp_eq
  beta_shp_ineq = estbounds_return$beta_shp_ineq
  solver = estbounds_return$solver
  
  #### Step 2: Construct the bounds
  # Default - Boolean variable of whether the answer to the scenario 1 is 
  # feasible or not
  bound0infe = FALSE
  
  ### Scenario 1: Estimate = FASLE, i.e. solve the exact problem
  if (estimate == FALSE){
    ub_shp0 = estbounds_original(A_obs, A_tgt, beta_tgt, beta_obs, A_shp_eq, 
                                 A_shp_ineq, beta_shp_eq, beta_shp_ineq, 
                                 "max", solver)
    lb_shp0 = estbounds_original(A_obs, A_tgt, beta_tgt, beta_obs, A_shp_eq, 
                                 A_shp_ineq, beta_shp_eq, beta_shp_ineq, 
                                 "min", solver)
    ub = ub_shp0$objval
    lb = lb_shp0$objval
    if (progress == TRUE){
      cat(sprintf("True bounds subject to shape constraints: [%s, %s]\n", 
                  lb_shp0$objval, ub_shp0$objval))      
    }
    
    # Store indicator of whether the estimation procedure should be conducted
    if (is.numeric(ub) == FALSE | is.numeric(lb) == FALSE){
      bound0infe = TRUE
      if (progress == TRUE){
        cat(paste("The original problem is infeasible. ",
                  "The estimated bounds will be displayed.", sep = ""))
      }
    } else {
      est = FALSE
    }
  }
  
  ### Scenario 2: Estimate = TRUE or scenario 1 is infeasible
  if (estimate == TRUE | bound0infe == TRUE){
    ## Update constraints
    # Zero lower bound
    lb_zero = rep(0,ncol(A_tgt))
    # Combine linear constraints 
    A1 = rbind(A_shp_eq, A_shp_ineq)
    rhs1 = rbind(beta_shp_eq, beta_shp_ineq)
    # Combine sense of models
    sense1 = rbind(rep("=", nrow(A_shp_eq)))
    if (is.null(A_shp_ineq) == FALSE){
      sense1 = rbind(sense1, rep("<=", nrow(A_shp_ineq)))
    }
    
    ## Solve model
    if (lnorm == 1){
      ## L1-norm
      stop("Work-in-progress")
    } else if (lnorm == 2){
      ## L2-norm
      # Stage one of the problem
      estbounds21 = estbounds1_L2(A_obs, beta_obs, A1, rhs1, sense1, 
                                  lb_zero, solver)
      # Stage two of the problem
      estbounds_ub = estbounds2_L2(estbounds21, A_obs, beta_obs, "max", 
                                   kappa, solver)
      estbounds_lb = estbounds2_L2(estbounds21, A_obs, beta_obs, "min", 
                                   kappa, solver)
    } else if (lnorm == "sup"){
      ## sup-norm
      stop("Work-in-progress")      
    }
    
    # Store results
    ub = estbounds_ub$objval
    lb = estbounds_lb$objval
    
    ## Print results
    if (progress == TRUE){
      cat(sprintf("Estimated bounds subject to shape constraints: [%s, %s]\n", 
                  round(lb, digits = 5), 
                  round(ub, digits = 5))) 
    }
    ## Store indicator variable that the result is estimated 
    est = TRUE
  }
  
  #### Step 3: Assign the return list and define class of output
  output = list(ub = ub,
                lb = lb,
                est = est,
                call = call,
                lnorm = lnorm)
  
  attr(output, "class") = "estbounds"
  
  invisible(output)
}

#' Computes the true bounds with shape contraints
#' 
#' @description The function computes the true bound subject to the shape 
#'    constraints without approximation.
#' 
#' @param original_sense Sense of the contraints to compute the true bound.
#' @inheritParams dkqs
#' @inheritParams invertci
#' @inheritParams estbounds
#' 
#' @return Returns the solution to the linear program.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#' 
#' @export
#' 
estbounds_original <- function(A_obs, A_tgt, beta_tgt, beta_obs, A_shp_eq, 
                               A_shp_ineq, beta_shp_eq, beta_shp_ineq,
                               original_sense, solver){
  
  #### Step 1: Problem set-up
  A_original = rbind(A_obs, A_shp_eq, A_shp_ineq)
  beta_original = rbind(beta_obs, beta_shp_eq, beta_shp_ineq)
  
  # Sense contraints
  sense_original = c(rep("=", nrow(A_obs)), rep("=", nrow(A_shp_eq)))
  if (is.null(A_shp_ineq) == FALSE){
    sense_original = c(sense_original, rep("<=", nrow(A_shp_ineq)))
  }
  # Zero lower bound
  lb_zero = rep(0, ncol(A_tgt))
  
  #### Step 2: Formulate the argument for optimization
  oarg = list(Af = NULL,
              bf = A_tgt,
              nf = NULL,
              A = A_original,
              rhs = beta_original,
              sense = sense_original,
              modelsense = original_sense,
              lb = lb_zero,
              qc = NULL)
  

  #### Step 3: Solve the model
  ans = do.call(solver, oarg)
  
  #### Step 4: Return result
  invisible(list(objval = ans$objval,
                 x = ans$x))
}

#' Estimates the bounds with shape contraints (Stage 1 with \eqn{\ell^2}-norm)
#' 
#' @description This function evaluates the solution to stage 1 of the 
#'    two-step procedure obtaining the estimated bound. \eqn{\ell^2}-norm 
#'    is used in the objective function.
#' 
#' @param A1 Constraint matrix for step 1 in the linear program of bound
#'    estimation.
#' @param rhs1 RHS vector for step 1 in the linear program of bound
#'    estimation.
#' @param sense1 Sense vector for step 1 in the linear program of bound
#'    estimation.
#' @inheritParams dkqs
#' 
#' @return Returns the solution to the first step of the two-step procedure 
#'    and argument for the linear/quadratic program.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{larg}{Arguments for the linear/quadratic program.}
#' 
#' @export
#' 
estbounds1_L2 <- function(A_obs, beta_obs, A1, rhs1, sense1, lb, solver){
  #### Step 1: Set up argument for optimizer
  l2_arg = list(Af = A_obs,
                bf = beta_obs,
                nf = 1,
                A = A1,
                rhs = rhs1,
                sense = sense1,
                modelsense = "min",
                lb = lb,
                qc = NULL)
  
  #### Step 2: Solve the model
  ans = do.call(solver, l2_arg) 
  
  #### Step 3: Return results
  invisible(list(objval = ans$objval, 
                 x = ans$x,
                 larg = l2_arg))
}

#' Estimates the bounds with shape contraints (Stage 2 with \eqn{\ell^2}-norm)
#' 
#' @description This function evaluates the solution to stage 2 of the 
#'    two-step procedure obtaining the estimated bound. \eqn{\ell^2}-norm 
#'    is used in the constraint
#' 
#' @param firststepsoln List of solutions to the first step problem.
#' @inheritParams gurobi_optim
#' @inheritParams estbounds
#' @inheritParams dkqs
#' 
#' @return Returns the solution to the second step of the two-step procedure.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#' 
#' @export
#' 
estbounds2_L2 <- function(firststepsoln, A_obs, beta_obs, modelsense, 
                          kappa, solver){
  #### Step 1: Extract information from the first-stage solution
  Qhat = firststepsoln$objval
  larg = firststepsoln$larg
  
  #### Step 2: Construct the quadratic inequality constraint
  step2_qc = list()
  if (is.null(A_obs) == FALSE){
    step2_qc$Qc = t(A_obs) %*% A_obs  
    step2_qc$q = -2*t(A_obs) %*% beta_obs
    step2_qc$rhs = Qhat * (1+kappa) - t(beta_obs) %*% beta_obs
    step2_qc$sense = "<="
  } else {
    step2_qc = NULL
  }
  qc_all = list(step2_qc)
  
  # Update quadratic constraint
  larg$qc = qc_all
  
  #### Step 3: Update objective function
  larg$Af = 0
  larg$bf = A_tgt
  larg$nf = 0
  
  #### Step 4: Update model sense based on max or min in step 2
  larg$modelsense = modelsense
  
  #### Step 5: Solve the model
  step2_ans = do.call(solver, larg)
  
  #### Step 6: Return results
  return(list(objval = step2_ans$objval,
              x = step2_ans$x))
}

#' Checks and updates the input of the function \code{estbounds}
#' 
#' @description This function checks and updates the input from the user for 
#'    the function \code{estbounds}. If there is any invalid input, this  
#'    function will terminate the procedure and generate appropriate error 
#'    messages.
#'    
#' @inheritParams estbounds
#' @inheritParams dkqs
#' @inheritParams invertci
#' 
#' @return Returns the updated value of the parameters back to the function 
#'    \code{estbounds}. The following information are updated:
#'    \itemize{
#'       \item{\code{df}}
#'       \item{\code{A_obs}}
#'       \item{\code{A_tgt}}
#'       \item{\code{beta_obs}}
#'       \item{\code{beta_tgt}}
#'       \item{\code{A_shp_eq}}
#'       \item{\code{beta_shp_eq}}
#'       \item{\code{A_shp_ineq}}
#'       \item{\code{beta_shp_ineq}}
#'       \item{\code{solver}}
#'    }
#' 
#' @export
#' 
estbounds_check <- function(df, func_obs, A_obs, A_tgt, beta_tgt,
                            A_shp_eq, A_shp_ineq, beta_shp_eq, beta_shp_ineq,
                            kappa, lnorm, solver, estimate, progress){

  ### Part 1. Check the data frame
  if (class(df) %in% c("data.frame", "matrix") == TRUE){
    df = as.data.frame(df)  
  } else {
    stop(gsub("\\s+", " ",
              "The data povided 'df' must either be a data.frame, a data.table, 
               or a matrix."), call. = FALSE)    
  }
  
  ### Part 2. Check the function
  if (class(func_obs) != "function"){
    stop("The input of 'func_obs' has to be a function.", call. = FALSE)
  } else{
    beta_obs = func_obs(df)
    beta_obs = as.matrix(beta_obs)
    # Check if the output is numeric
    if (is.numeric(beta_obs[,1]) == FALSE){
      stop("The output of 'func_obs' has to be numeric.", call. = FALSE)
    } else{
      if (dim(beta_obs)[2] != 1){
        stop("The output of 'func_obs' has to be a column vector", 
             call. = FALSE)
      } else if (dim(beta_obs)[1] != dim(A_obs)[1]){
        stop("The number of rows in the output of 'func_obs' has to be the 
             same as the number of rows in the argument 'A_obs'.", 
             call. = FALSE)
      }
    }
  }
  
  #### Part 3. Check beta_tgt
  if (!(is.numeric(beta_tgt) == TRUE & length(beta_tgt) == 1)) {
    stop("The argument 'beta_tgt' must be a scalar.", call. = FALSE)
  }  
  
  #### Part 4. Check matrices and vectors
  # Check A_shp_eq and beta_shp_eq
  eq_return = estbounds_check_Ab(A_shp_eq, beta_shp_eq, "A_shp_eq", 
                                 "beta_shp_eq")
  A_shp_eq = eq_return$A_updated
  b_shp_eq = eq_return$b_updated
  # Check A_shp_ineq and beta_shp_ineq
  ineq_return = estbounds_check_Ab(A_shp_ineq, beta_shp_ineq, "A_shp_ineq", 
                                   "beta_shp_ineq")
  A_shp_ineq = ineq_return$A_updated
  b_shp_ineq = ineq_return$b_updated
  # Check A_tgt and beta_tgt
  tgt_return = estbounds_check_Ab(A_tgt, beta_tgt, "A_tgt", "beta_tgt")
  A_tgt = tgt_return$A_updated
  beta_tgt = tgt_return$b_updated
  # Check A_obs and beta_obs
  obs_return = estbounds_check_Ab(A_obs, beta_obs, "A_obs", "beta_obs")
  A_obs = obs_return$A_updated
  beta_obs = obs_return$b_updated

  #### Step 5: Check 'kappa'
  if (!(is.numeric(kappa) == TRUE & length(kappa) == 1 & kappa >= 0)) {
    stop("The argument 'kappa' must be a nonnegative scalar.", call. = FALSE)
  } 
  
  #### Step 6: Check 'lnorm'
  if (lnorm != 1 & lnorm != 2 & lnorm != "sup"){
    stop("Only 1-norm, 2-norm and sup-norm is supported in this function.", 
         call. = FALSE)
  }
  
  #### Part 7: Check solver - only 'gurobi' can be used to obtain the bounds 
  #### of the shape restriction
  # Turn the name of solver to lower case 
  solver = tolower(solver)
  
  # Case 1: If no solver name is provided by the user
  if (is.null(solver) == TRUE){
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
      solver = "gurobi"
    } else {
      stop(gsub("\\s+", " ",
                paste0("This function is only incompatible with 'gurobi'. 
                       Please install 'gurobi' to obtain the bounds of the 
                       problem subject to shape restriction: gurobi 
                       (version 8.1-1 or later).")), call. = FALSE)
    }
    if (progress == TRUE){
      cat(paste("No solver solver is suggested by the user. The solver", 
                "'gurobi' is automatically selected.\n", sep = ""))
    }
  } else if (solver == "gurobi"){
    # Case 2: If the user specified the solver as 'gurobi'
    solver = gurobi_optim
  } else {
    # Case 3: If the user specified a solver that is not compatible
    stop(gsub("\\s+", " ",
              paste0("This function is only incompatible with ", solver, 
                     ". Please install 'gurobi' to obtain the bounds of the 
                     problem subject to shape restriction: gurobi 
                     (version 8.1-1 or later).")), call. = FALSE)
  }
  
  #### Step 8 Check estimate
  if (!(estimate == TRUE | estimate == FALSE)){
    stop("The argument 'estimate' has to be boolean.", call. = FALSE)
  }
  
  #### Step 9. Check progress
  if (!(progress == TRUE | progress == FALSE)){
    stop("The argument 'progress' has to be boolean.", call. = FALSE)
  }
  
  #### Step 10. Return the updated information
  invisible(list(df = df,
                 A_obs = A_obs,
                 A_tgt = A_tgt,
                 beta_obs = beta_obs,
                 beta_tgt = beta_tgt,
                 A_shp_eq = A_shp_eq,
                 A_shp_ineq = A_shp_ineq,
                 beta_shp_eq = beta_shp_eq,
                 beta_shp_ineq = beta_shp_ineq,
                 solver = solver))
}

#' Check the constraint matrix and the corresponding vector
#' 
#' @description This function checks the constraint matrix \eqn{\bm{A}} and 
#'    the corresponding vector \eqn{\bm{\beta}}
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
estbounds_check_Ab <- function(A, b, Aname, bname){
  if (is.null(A) + is.null(b) == 1){
    #### Step 1: Check if A and b must be both NULL or both non-NULL
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
    #### Step 2: Both A and b are NULL. Nothing else to do
  } else {
    #### Step 3:  Checks for the case where both A and b are non-NULL
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
      stop(sprintf("The argument '%s' cannot have more than 1 column", bname), 
           call. = FALSE)
    }
    
    ## Part 4: Ensure that the number of rows of A and b are identical
    if (nrow(A) != length(b)){
      stop(sprint("The number of rows of '%s' has to be equal to the number
                  of rows of '%s.", Aname, bname), call. = FALSE)
    }
    
    #### Step 4: Update class
    # Ensure that both A and b to ensure that they are both matrices
    A = matrix_list[[1]]
    b = matrix_list[[2]]
  }
  
  #### Step 5: Return results
  return(list(A_updated = A,
              b_updated = b))
}

#' Print results from \code{estbounds}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{estbounds}.
#'    
#' @param x Object returned from \code{estbounds}.
#' @param ... Additional arguments.
#' 
#' @return Nothing is returned. This function prints results from 
#'    \code{estbounds}.
#' 
#' @export
#' 
print.estbounds <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  
  if (x$est == TRUE){
    #### Case 1: Report the estimated bounds
    if (is.numeric(x$lnorm) == TRUE){
      cat(sprintf("Norm used in optimization problem: L%s-norm \n", x$lnorm))
    } else {
      cat(sprintf("Norm used in optimization problem: %s-norm \n", x$lnorm)) 
    }
    cat(sprintf("Estimated bounds subject to shape constraints: [%s, %s] \n", 
                round(x$lb, digits = 5), round(x$ub, digits = 5)))
  } else {
    #### Case 2: Report the true bounds
    cat(sprintf("True bounds subject to shape constraints: [%s, %s] \n", 
                x$lb, x$ub)) 
  }
}

#' Summary of results from \code{estbounds}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{estbounds}. This is a wrapper of the function
#'    \code{print.estbounds}.
#'    
#' @param x Object returned from \code{estbounds}.
#' @param ... Additional arguments.
#' 
#' @return Nothing is returned. This function prints results from 
#'    \code{estbounds}.
#' 
#' @export
#' 
summary.estbounds <- function(x, ...){
  #### Call theprint function
  print(x)
}


       