#' Computes the \eqn{p}-value of the subsampling procedure
#' 
#' @description This function conducts inference and returns the 
#'   \eqn{p}-value using the subsampling procedure.
#' 
#' @import gurobi cplexAPI Rcplex Momocs limSolve foreach doMC parallel Momocs
#' 
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @param func_var Function that generates the asymptotic variance 
#'   matrix of the estimator \eqn{\hat{\beta}_{\mathrm{obs}}}.
#' @param phi Power for the sample. \eqn{n^\phi} represents the size
#'   of each subsample.
#' @param alpha Significance level.
#' 
#' @return Returns a list of output calculated from the function:
#'   \item{p_val}{\eqn{p}-value.}
#'   \item{decision}{Decision of the test.}
#'   \item{alpha}{Significance level.}
#'   \item{T_n}{Test statistic \eqn{T_n}.}
#'   \item{T_sub}{The list of test statistics from the subsampling procedure.}
#'   \item{solver}{Solver used in solving the linear and quadratic programs.}
#'   \item{cores}{Number of cores used.}
#'   \item{call}{The function that has been called.}
#'   \item{norm}{Norm used.} 
#' 
#' @export
#' 
subsample <- function(data, A_obs, func_obs, func_var, 
                      A_shp, beta_shp, A_tgt, beta_tgt, 
                      bs_seed = 1, R = 100, p_sig = 2, solver = NULL, 
                      cores = 8, lnorm = 2, phi = 2/3, alpha = .05,
                      progress = FALSE){
  
  # = = = = = = 
  # Step 1: Obtain call, check and update the dependencies
  # = = = = = = 
  ## Obtain the call information
  call = match.call()
  
  ## Check and update 
  checkupdate = subsample_check(data, A_obs, func_obs, func_var, 
                                A_shp, beta_shp, A_tgt, beta_tgt, 
                                bs_seed, R, p_sig, solver, cores, lnorm, 
                                phi, progress)
  
  ## Update information obtained from check
  # Data frame
  data = checkupdate$data
  # Matrices and vectors
  A_obs = checkupdate$A_obs
  A_shp = checkupdate$A_shp
  A_tgt = checkupdate$A_tgt
  beta_shp = checkupdate$beta_shp
  beta_tgt = checkupdate$beta_tgt
  # Solver
  solver = checkupdate$solver

  # = = = = = =
  # Step 2: Solve for T_n
  # = = = = = =
  ## Solve the main problem with the full sample
  Treturn = subsample_prob(data, func_obs, func_var, A_obs, A_shp, A_tgt, 
                           beta_shp, beta_tgt, lnorm, solver)
  
  # = = = = = = 
  # Step 3: Subsampling procedure
  # = = = = = =  
  n = nrow(data)
  m = floor(n^(phi))
  if (cores == 1){
    # One core
    T_subsample = subsample_onecore(data, bs_seed, R, func_obs, func_var, 
                                    A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                                    lnorm, solver, progress, m)
    
  } else {
    # Many cores
    T_subsample = subsample_manycores(data, bs_seed, R, func_obs, func_var, 
                                      A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                                      lnorm, solver, progress, m)
  }
  
  # = = = = = = 
  # Step 4: Compute the p-value (using the p_eval function in dkqs)
  # = = = = = = 
  pval_return = p_eval(T_subsample$T_sub, Treturn$objval, p_sig, alpha)
  p_val = pval_return$p
  decision = pval_return$decision
  
  # = = = = = = 
  # Step 5: Close the progress bar that is used in the subsampling procedure
  # = = = = = = 
  if (progress == TRUE){
    close(T_subsample$pb) 
    cat("                                            ")
  }
  
  # = = = = = = 
  # Step 6: Assign the return list
  # = = = = = = 
  output = list(p_val = as.numeric(p_val),
                decision = decision,
                alpha = alpha,
                T_n = as.numeric(Treturn$objval),
                T_sub = T_subsample$T_sub,
                solver = checkupdate$solver_name,
                cores = cores,
                call = call,
                norm = lnorm)
  
  attr(output, "class") = "subsample"
  
  return(output)
}

#' Formulate and solve the subsampling problem
#' 
#' @description Based on the sample data given by the data frame \code{data}, 
#'   this function formulates and solves linear/quadratic program in the 
#'   subsampling procedure.
#'   
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' 
#' @return Returns a list of output that are obtained from the optimizer:
#'   \item{x}{Optimal point calculated from the optimizer.}
#'   \item{objval}{Optimal value calculated from the optimizer.}
#'   \item{larg}{List of arguments passed to the optimizer.}
#'   \item{beta}{The beta vector \eqn{\widehat{\bm{\beta}}_{\mathrm{obs}}}
#'     used in the optimization problem that is calculated from the 
#'     \code{func_obs} function with data \code{data}.}
#'   \item{omega}{The Omega matrix \eqn{\widehat{\bm{\Omega}}_n} used in the 
#'     optimization problem that is calculated from the \code{func_obs} 
#'     function with data \code{data}.}
#'     
#' @export 
#' 
subsample_prob <- function(data, func_obs, func_var, A_obs, A_shp, A_tgt, 
                           beta_shp, beta_tgt, lnorm, solver){
  # = = = = = = 
  # Step 1: Obtain parameters from the data frame
  # = = = = = = 
  n = nrow(data)
  beta_obs_hat = func_obs(data)
  omega_hat = func_var(data)
  
  # = = = = = = 
  # Step 2: Define the inverse omega matrix
  # = = = = = = 
  # Obtain the inverse of the diagonal etnreis
  diag_omega = diag(omega_hat)
  g = 1/diag_omega
  # Replace the entries by 0 for those that are equal to zero in 1/Omega
  g[diag_omega == 0] = 0
  G = diag(g)
  # Create the new A and b matrices
  GA = G %*% A_obs
  Gb = G %*% beta_obs_hat
  
  # = = = = = = 
  # Step 3: Form the objective function and constraints
  # = = = = = = 
  # Model sense
  modelsense_new = "min"
  # Set the objective function and constraints
  if (lnorm == 1){
    ### L1-norm
    # Objective function - cost matrix
    c_new = c(rep(0, ncol(A_obs)), rep(1, k), rep(-1, k))
    # Constraints
    A_zero = matrix(rep(0, k^2), nrow = k)
    A1_shp = cbind(A_shp, A_zero, A_zero)
    A1_tgt = cbind(A_tgt, A_zero, A_zero)
    A1_obs = cbind(GA, -diag(k), -diag(k))
    A_new = rbind(A1_shp, A1_tgt, A1_obs)
    # RHS vector
    rhs_new = c(beta_shp, beta_tgt, Gb)
    # Lower bounds
    lb_new = rep(0, length(c_new))
    # Sense
    sense_new = rep("=", nrow(A_new))
    
    # Set the list to pass to the solver
    l_arg = list(Af = NULL,
                 bf = c_new,
                 nf = sqrt(n),
                 A = A_new,
                 rhs = rhs_new,
                 sense = sense_new,
                 modelsense = modelsense_sense,
                 lb = lb_new)
    
  } else if (lnorm == 2){
    ### L2-norm
    # Constraints
    A_new = rbind(A_shp, A_tgt)
    # RHS vector
    rhs_new = c(beta_shp, beta_tgt)
    # Lower bounds
    lb_new = rep(0, ncol(A_shp))
    # Sense
    sense_new = rep("=", nrow(A_new))
    
    # Set the list to pass to the solver
    l_arg = list(Af = GA,
                 bf = Gb,
                 nf = sqrt(n),
                 A = A_new,
                 rhs = rhs_new,
                 sense = sense_new,
                 modelsense = modelsense_new,
                 lb = lb_new)
  }
  
  # = = = = = = 
  # Step 4: Solve the model and return the results
  # = = = = = = 
  # Solve the model
  ans = do.call(solver, l_arg)
  
  # Return the results
  invisible(list(x = ans$x,
                 objval = ans$objval,
                 larg = l_arg,
                 beta = beta_obs_hat,
                 omega = omega_hat))
}

#' Subsampling procedure without parallel programming
#' 
#' @description This function computes the list of test statistics that are 
#'   obtained from the subsampling procedure without using parallel 
#'   programming.
#'   
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample 
#' @param m Size of each subsample.
#' 
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T_sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#' 
#' @export
#' 
subsample_onecore <- function(data, bs_seed, R, func_obs, func_var, 
                              A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                              lnorm, solver, progress, m){
  # = = = = = = 
  # Step 1: Initialize the vectors and the progress bar
  # = = = = = = 
  # Initialize the vectors
  T_sub = NULL
  beta_sub = NULL
  # Initialize the progress bar
  if (progress == TRUE){
    pb = utils::txtProgressBar(min = 0, max = R, style = 3, width = 20)
    cat("\r") 
  } else {
    pb = NULL
  }
  
  # = = = = = = 
  # Step 2: Conduct the subsampling procedure
  # = = = = = = 
  for (i in 1:R){
    ## (2.1) Set the seed
    set.seed(bs_seed + i)
    ## (2.2) Draw the subsample
    df_sub = as.data.frame(Momocs::sample_n(data, m, replace = FALSE))
    # Re-index the rows
    rownames(df_sub) = 1:nrow(df_sub)
    ## (2.3) Compute the bootstrap estimates
    # Compute the value of beta_bs_star using the function func_obs
    sub_return = subsample_prob(df_sub, func_obs, func_var, 
                                A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                                lnorm, solver)
    T_sub = c(T_sub, sub_return$objval)
    beta_sub = cbind(beta_sub, sub_return$beta)
    ## (2.4) Update progress bar
    if (progress == TRUE){
      if (i != R){
        utils::setTxtProgressBar(pb, i)
        cat("\r\r")  
      } else {
        utils::setTxtProgressBar(pb, i)
        cat("\r\b")          
      }
    }
  }
  
  # = = = = = = 
  # Step 3: Return the results
  # = = = = = = 
  return(list(T_sub = T_sub,
              pb = pb))
}

#' Subsampling procedure with parallel programming
#' 
#' @description This function computes the list of test statistics that are 
#'   obtained from the subsampling procedure using parallel programming.
#'   
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample 
#' @param m subsample_onecore
#' 
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T_sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#' 
#' @export
#' 
subsample_manycores <- function(data, bs_seed, R, func_obs, func_var, 
                                A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                                lnorm, solver, progress, m){
  # = = = = = = 
  # Step 1: Initialize the parallel programming package
  # = = = = = = 
  options(warn=-1)
  # Assign dopar
  `%dopar%` = foreach::`%dopar%`
  # Register core
  doMC::registerDoMC(cores)
  
  # = = = = = = 
  # Step 2: Initialize the vectors and the progress bar
  # = = = = = = 
  # Initialize the vectors
  T_sub = NULL
  beta_sub = NULL
  # Initialize the progress bar
  if (progress == TRUE){
    # Initialize the counter
    cl = PtProcess::makeSOCKcluster(8)
    doSNOW::registerDoSNOW(cl)
    # Set the counter and progress bar
    pb = utils::txtProgressBar(max=R, style=3, width = 20) 
    
    cat("\r")
    progress <- function(n){
      utils::setTxtProgressBar(pb, n) 
      if (n < R){
        cat("\r\r") 
      } else {
        cat("\r\b")     
      }
    }
    opts = list(progress=progress) 
  } else {
    pb = NULL
    opts = NULL
  }
  
  # Comb function for using parallel programming
  comb <- function(x, ...) {
    lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), 
                                                      function(y) y[[i]])))
  }
  
  # = = = = = = 
  # Step 3: Conduct the subsampling procedure
  # = = = = = = 
  listans = foreach::foreach(i=1:R, .multicombine = TRUE, 
                             .combine="comb", .options.snow=opts,
                             .packages='linearprog') %dopar% {
    ## (3.1) Set the seed
    set.seed(bs_seed + i)
    ## (3.2) Draw the subsample
    df_sub = as.data.frame(Momocs::sample_n(data, m, replace = FALSE))
    # Re-index the rows
    rownames(df_sub) = 1:nrow(df_sub)
    ## (3.3) Compute the bootstrap estimates
    # Compute the value of beta_bs_star using the function func_obs
    sub_return = subsample_prob(df_sub, func_obs, func_var, 
                                A_obs, A_shp, A_tgt, beta_shp, beta_tgt, 
                                lnorm, solver)
    T_sub = data.frame(sub_return$objval)
    beta_sub = data.frame(c(sub_return$beta))
    ## (3.4) Combine the results
    list(T_sub, beta_sub)
  }
  
  # = = = = = = 
  # Step 4: Retrieve results from output
  # = = = = = = 
  T_sub = as.vector(unlist(listans[[1]]))
  beta_sub = data.frame(matrix(as.matrix(listans[[2]]), 
                               nrow = length(func_obs(data)), 
                               byrow = FALSE))
  
  # = = = = = = 
  # Step 5: Return the results
  # = = = = = = 
  return(list(T_sub = T_sub,
              pb = pb))
}

#' Print results from \code{subsample}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{subsample}.
#'    
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#' 
#' @return Print the basic set of results from \code{subsample}.
#' 
#' @export
#' 
print.subsample <- function(x, ...){
  cat("\r\r")
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  if (x$decision == 1){
    cat(sprintf("The null hypothesis is rejected at the %s level.\n\n", 
                paste0(x$alpha*100, "%")))
  } else {
    cat(sprintf("The null hypothesis cannot be rejected at the %s level.\n\n", 
                paste0(x$alpha*100, "%")))
  }
  cat(sprintf("Test statistic: %s.             \n", round(x$T_n, digits = 5)))  
  cat(sprintf("p-value: %s.\n", round(x$p_val, digits = 5)))
  cat(sprintf("Linear and quadratic programming solver used: %s.\n", x$solver))
  if (x$norm == 1){
    cat(sprintf("Norm used in the optimization problem: L1-norm.\n")) 
  } else if (x$norm == 2){
    cat(sprintf("Norm used in the optimization problem: L2-norm.\n")) 
  }
  cat(sprintf("Number of cores used: %s.\n", x$cores))
}

#' Summary of results from \code{subsample}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{subsample}. This is a wrapper of the \code{print} command.
#'    
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{subsample}.
#' 
#' @export
#' 
summary.subsample <- function(x, ...){
  print(x)
}

#' Checks and updates the input from \code{subsample}
#'
#' @description This function checks and updates the input of the user. If 
#'    there is any invalid input, this function will be terminated and 
#'    generates appropriate error messages. This function is mainly calling
#'    other functions in the \code{checks} files to conduct the checks and
#'    updates.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' 
#' @return Returns the updated value of the parameters back to the function 
#'    \code{subsample}. 
#'   \item{data}{Upated data in class \code{data.frame}}
#'   \item{A_obs}{Updated "observed" matrix in the \code{matrix} class.}
#'   \item{beta_obs}{Obtain \eqn{\widehat{\bm{\beta}}_{\mathrm{obs}}} 
#'      that is obtained from the function \code{func_obs}.}
#'   \item{A_shp}{Updated matrix for shape constraints in the \code{matrix} 
#'     class.}
#'   \item{beta_shp}{Updated RHS vector for shape constraints in the 
#'     \code{matrix} class.}
#'   \item{A_tgt}{Updated "target" matrix in the \code{matrix} class.}
#'   \item{beta_tgt}{Updated "target" vector in the \code{matrix} class}
#'   \item{solver}{Name of the function that corresponds to the optimizer.}
#'   \item{solver_name}{Updated name of solver in lower case.}
#'   \item{cores}{Updated number of cores to be used in the parallelization
#'      of the for-loop in the bootstrap procedure.}
#' 
#' @export
#' 
subsample_check <- function(data, A_obs, func_obs, func_var, 
                            A_shp, beta_shp, A_tgt, beta_tgt, bs_seed, R, 
                            p_sig, solver, cores, lnorm, phi, progress){
  
  # = = = = = = 
  # Step 1: Conduct the checks
  # = = = = = = 
  # Check the data frame
  data = check_dataframe(data, "data")
  
  # Check function and export output
  beta_obs = check_func(func_obs, A_obs, data, "func_obs", "A_obs", "col")
  check_func(func_var, A_obs, data, "func_var", "A_obs", "square")
  
  # Check the matrices A vs b
  shp_return = check_Ab(A_shp, beta_shp, "A_shp", "beta_shp")
  tgt_return = check_Ab(A_tgt, beta_tgt, "A_tgt", "beta_tgt")
  obs_return = check_Ab(A_obs, beta_obs, "A_obs", "beta_obs")
  
  # Check if the variable is a positive integer
  check_positiveinteger(bs_seed, "bs_seed")
  check_positiveinteger(R, "R")
  check_positiveinteger(p_sig, "p_sig")
  check_positiveinteger(cores, "cores")
  
  # Check if lnorm is either 1 or 2
  check_norm(lnorm, "lnorm")
  
  # Check if phi is in the range (0,1)
  check_numrange(phi, "phi", "open", 0, "open", 1)
  
  # Check if the variable is boolean
  check_boolean(progress, "progress")
  
  # Check solveresssssss
  solver_return = check_solver(solver, "solver")
  
  # = = = = = = 
  # Step 2: Return the updated objects
  # = = = = = = 
  return(list(data = data,
              A_obs = obs_return$A,
              beta_obs = obs_return$b,
              A_shp = shp_return$A,
              beta_shp = shp_return$b,
              A_tgt = tgt_return$A,
              beta_tgt = tgt_return$b,
              solver = solver_return$solver,
              solver_name = solver_return$solver_name))
}