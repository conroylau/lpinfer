#' Compute p-value of a quadratic program
#'
#' @description This function runs the procedure suggested by Torgovitsky (2019)
#'    using the cone-tightening procedure proposed by Deb, Kitamura, Quah and
#'    Stoye (2018).
#'
#' @import slam gurobi car modelr
#'
#' @param df The dataframe that contains the sample data.
#' @param A_obs The observed matrix in the linear program.
#' @param A_tgt The "target matrix" in the linear program.
#' @param func_obs The function that generates the required 
#'    \eqn{\beta_{\mathrm{obs}}}.
#' @param beta_tgt The value of \eqn{t} in the null hypothesis.
#' @param bs_seed The starting value of the seed in bootstrap.
#' @param bs_num The total number of bootstraps.
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#' @param tau_input The value of tau chosen by the user.
#' @param lpsolver The name of the linear/quadratic programming package to 
#'    be used to obtain the solutions to the linear/quadratic programs. This
#'    function supports \code{'gurobi'}, \code{'cplexapi'} and 
#'    \code{'quadprog'}.
#'
#' @returns Returns the \eqn{p}-value, the value of tau used \eqn{\tau^\ast}, 
#'   test statistic \eqn{T_n(\tau_n)} and the list of bootstrap test statistics 
#'   \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}. 
#'
#' @export
# library(slam)
# library(gurobi)
# library(car)
dkqs_cone <- function(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed = 1,
                      bs_num = 100, p_sig = 2, tau_input = .5, lpsolver = NULL){
  
  #### Step 1: Check the dependencies
  checkupdate = dkqs_cone_check(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                                bs_num, p_sig, tau_input, lpsolver)
  # Update and return the information returned from the function dkqs_cone_check
  df = checkupdate$df
  A_obs = checkupdate$A_obs
  A_tgt = checkupdate$A_tgt
  lpsolver = checkupdate$lpsolver
  # Display the lpsolver being used
  cat(paste("LP solver used: ", lpsolver, ".\n", sep = ""))
  
  #### Step 2: Initialization
  # Initialization
  N = nrow(df)
  J = length(unique(df[,"Y"])) - 1
  # Compute beta_obs_hat using the function defined by user
  beta_obs_hat = checkupdate$beta_obs_hat
  # Assign the solver used
  if (lpsolver == "gurobi"){
    f_solver = gurobi_optim
  } else if (lpsolver == "cplexapi"){
    f_solver = cplexapi_optim
  } else if(lpsolver == "quadprog"){
    f_solver = quadprog_optim
  }

  #### Step 3: Choose the value of tau
  tau_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "tau", N, 
                         f_solver)$objval
  if (tau_input > tau_return){
    tau = tau_return
  } else if (tau_input <= tau_return){
    tau = tau_input
  } else {
    # Error message when the problem is infeasible.
    stop("The problem is infeasible. Choose other values of tau.")
  }

  #### Step 4: Solve QP (4) in Torgovitsky (2019)
  full_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "test", N, 
                          f_solver)
  x_star = full_return$x
  # Return and stop program if it is infeasible
  if (is.null(x_star) == TRUE){
    stop("The problem is infeasible. Choose other values of tau.")
  }
  s_star = A_obs %*% x_star
  # T_star is the test statistic used
  T_star = full_return$objval
  
  #### Step 5: Compute the bootstrap estimates
  # T_bs is the list of bootstrap test statistics used
  T_bs = beta_bs(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt, func_obs, 
                 beta_obs_hat, beta_tgt, tau, N, f_solver)
  #### Step 6: Compute the p-value
  # decision = 1 refers to rejected, decision = 0 refers to not rejected
  p_val = p_eval(T_bs, T_star, p_sig)
  tau = round(tau, digits = 5)
  
  cat(paste("-----------------------------------", "\n"))
  cat(paste("p-value: ", p_val, ".\n", sep = ""))
  cat(paste("Value of tau used in LP: ", tau, ".\n", sep = ""))
  invisible(list(p = p_val, tau = tau, T_star = T_star, T_bs = T_bs))
}

#' Formulating and solving quadratic programs
#'
#' @description This function formulates and solves the quadratic programs (4)
#'    and (6) in Torgovitsky (2019).
#'
#' @param tau The RHS vector fo the linear constraints.
#' @param problem The sense of the linear constraints.
#' @inheritParams dkqs_cone
#'
#' @returns Returns the solution to the quadratic program.
#'
#' @export
prog_cone <- function(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, problem, n,
                      f_solver){
  #### Step 1: Formulation of the objective function for (5)
  rn = dim(A_tgt)[1]
  cn = dim(A_tgt)[2]
  obj2 = t(A_obs) %*% A_obs * n
  obj1 = -2 * t(A_obs) %*% beta_obs_hat * n
  obj0 = t(beta_obs_hat) %*% beta_obs_hat * n

  #### Step 2: Formulation of constraints
  ones = matrix(rep(1, cn), nrow = 1)
  lb = matrix(rep(0, cn), nrow = 1)
  # Obtain parameters
  theta_up = gurobi_optim(NULL, A_tgt, NULL, ones, c(1), "=", "max", lb)
  theta_down = gurobi_optim(NULL, A_tgt, NULL, ones, c(1), "=", "min", lb)
  # Obtain required set of indices
  x_fullind = 1:cn
  ind_up = which(A_tgt %in% theta_up$objval)
  ind_down = which(A_tgt %in% theta_down$objval)
  ind_0 = x_fullind[-c(ind_up, ind_down)]
  # Updated lb for certain indices
  rhs_up = (beta_tgt - theta_down$objval) * tau / length(c(ind_0, ind_up))
  rhs_down = (theta_up$objval - beta_tgt) * tau / length(c(ind_0, ind_down))
  rhs_0 = (1 - rhs_up * length(ind_up) - rhs_down * length(ind_down)) *
            tau / length(ind_0)
  # LP rhs and sense
  lp_rhs = c(beta_tgt, 1)
  lp_sense = c("=", "=")

  #### Step 3: Solve the QP
  # If problem == "test", this function solves linear program (4)
  # If problem == "bootstrap", this function solves linear program (5)
  # If program == "tau", this function solves linear program (6)
  if (problem == "test"){
    ans = f_solver(obj2, obj1, obj0, rbind(A_tgt, ones), lp_rhs, lp_sense,
                       "min", lb)
  } else if (problem == "bootstrap"){
  # Update lb
    lb_new = lb
    lb_new[ind_up] = rhs_up
    lb_new[ind_down] = rhs_down
    lb_new[ind_0] = rhs_0
    ### Use Gurobi to find the optimum
    ans = f_solver(obj2, obj1, obj0, rbind(A_tgt, ones), lp_rhs, lp_sense,
                       "min", lb_new)
  } else if (problem == "tau"){
    # Update matrix A
    lp_lhs_tau = rbind(A_tgt, ones)
    lp_lhs_tau = cbind(matrix(c(0,0), nrow = 2), lp_lhs_tau)
    # Add one unit because of the additional position for tau
    len_tau = cn + 1
    lb_tau = rep(0, len_tau)
    lp_rhs_tau = lp_rhs
    lp_sense_tau = lp_sense
    # Inequality constraints for ind_up
    for (i in 1:length(ind_up)){
      new_const = tau_constraints(len_tau, rhs_up, -1, ind_up[i]+1, 0, "<=",
                                  lp_lhs_tau, lp_rhs_tau, lp_sense_tau)
      lp_lhs_tau = new_const$lp_lhs_tau
      lp_rhs_tau = new_const$lp_rhs_tau
      lp_sense_tau = new_const$lp_sense_tau
    }
    # Inequality constraints for ind_down
    for (i in 1:length(ind_down)){
      new_const = tau_constraints(len_tau, rhs_down, -1, ind_down[i]+1, 0, "<=",
                                  lp_lhs_tau, lp_rhs_tau, lp_sense_tau)
      lp_lhs_tau = new_const$lp_lhs_tau
      lp_rhs_tau = new_const$lp_rhs_tau
      lp_sense_tau = new_const$lp_sense_tau
    }
    # Inequality constraints for ind_0
    for (i in 1:length(ind_0)){
      new_const = tau_constraints(len_tau, rhs_0, -1, ind_0[i]+1, 0, "<=",
                                  lp_lhs_tau, lp_rhs_tau, lp_sense_tau)
      lp_lhs_tau = new_const$lp_lhs_tau
      lp_rhs_tau = new_const$lp_rhs_tau
      lp_sense_tau = new_const$lp_sense_tau
    }
    ans = gurobi_optim(NULL, c(1,rep(0, cn)), NULL, lp_lhs_tau, lp_rhs_tau,
                  lp_sense_tau, "max", lb_tau)
  }
  return(ans)
}

#' Gurobi solver for quadratic and linear programs
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using Gurobi.
#'
#' @param obj2 The matrix that corresponds to the second-order term in the
#'    quadratic program.
#' @param obj1 The vector that corresponds to the first-order term in the
#'    quadratic or linear program.
#' @param obj0 The constant term in the quadratic or linear program.
#' @param A The constraint matrix.
#' @param rhs The RHS vector fo the linear constraints.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#'
#' @returns Returns the solution to the quadratic or linear program.
#'
#' @export
gurobi_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  # Gurobi set-up
  model = list()
  model$Q = obj2
  model$obj = obj1
  model$objcon = obj0
  model$A = A
  model$rhs = rhs
  model$sense = sense 
  model$modelsense = modelsense
  model$lb = lb
  # Result of the LP or QP, and return result
  params = list(OutputFlag=0)
  gurobi_result = gurobi(model, params)
  return(gurobi_result)
}

#' cplexapi solver for quadratic and linear programs
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the \code{cplexAPI} package.
#'
#' @param obj2 The matrix that corresponds to the second-order term in the
#'    quadratic program.
#' @param obj1 The vector that corresponds to the first-order term in the
#'    quadratic or linear program.
#' @param obj0 The constant term in the quadratic or linear program.
#' @param A The constraint matrix.
#' @param rhs The RHS vector fo the linear constraints.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#'
#' @returns Returns a list of output from \code{cplexAPI}.
#'
#' @export
cplexapi_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  ### Step 1: Translate the notations that are used for gurobi into those to be 
  ### used in cplexAPI
  # Model sense
  modelsense[modelsense == "min"] == cplexAPI::CPX_MIN
  modelsense[modelsense == "max"] == cplexAPI::CPX_MAX
  # Inequality/equality signs
  sense[sense == "<="] == "L"
  sense[sense == ">="] == "G"
  sense[sense == "=="] == "E"
  # Bounds
  sense[sense == Inf] == cplexAPI::CPX_INFBOUND
  sense[sense == -Inf] == - cplexAPI::CPX_INFBOUND
  
  ### Step 2: cplexAPI set-up
  env = cplexAPI::openEnvCPLEX()
  cplexAPI::setDblParmCPLEX(env, 1016, 1e-06)
  prob = cplexAPI::initProbCPLEX(env)
  cplexAPI::chgProbNameCPLEX(env, prob, "sample")
  
  ### Step 3: Model set-up
  model = list()
  model$Q = obj2
  model$obj = obj1
  model$objcon = obj0
  model$A = A
  model$rhs = rhs
  model$sense = sense 
  model$modelsense = modelsense
  model$lb = lb
  
  ## ~~~~ Note:
  ## Check the notations of the followings
  ##
  ##
  
  ### Step 4: Solve LP/QP
  cplexAPI::copyLpwNamesCPLEX(env = env,
                              lp = prob,
                              nCols = ncol(lpobj$A),
                              nRows = nrow(lpobj$A),
                              lpdir = lpdir,
                              objf = lpobj$obj,
                              rhs = lpobj$rhs,
                              sense = sense,
                              matbeg = beg,
                              matcnt = cnt,
                              matind = ind,
                              matval = val,
                              lb = lb,
                              ub = ub)
  cplexAPI::lpoptCPLEX(env, prob)
  solution = cplexAPI::solutionCPLEX(env, prob)
  cplexAPI::delProbCPLEX(env, prob)
  cplexAPI::closeEnvCPLEX(env)
  if (typeof(solution) == "S4") {
    if (attr(solution, "class") == "cplexError") {
      status = 0
      solution = list()
      solution$objval = NA
      solution$x = NA
    }
  }  else {
    if (solution$lpstat == 1){
      status = 1
    }
    if (solution$lpstat != 1){
      status = 0
    }
  }
  return(list(objval = solution$objval,
              x = solution$x))
}

#' quadprog solver for quadratic and linear programs
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the \code{quadprog} package.
#'    
#' @require
#'
#' @param obj2 The matrix that corresponds to the second-order term in the
#'    quadratic program.
#' @param obj1 The vector that corresponds to the first-order term in the
#'    quadratic or linear program.
#' @param obj0 The constant term in the quadratic or linear program.
#' @param A The constraint matrix.
#' @param rhs The RHS vector fo the linear constraints.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#'
#' @returns Returns a list of output from \code{quadprog}.
#'
#' @export
quadprog_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  ### Step 1: Set the objective function
  # Set objective function
  Dmat = 2 * obj2
  dvec = obj1
  # Negate the objective function if modelsense == max because the default of
  # solve.QP is to minimize the objective function
  if (modelsense == "max"){
    Dmat = - Dmat
    dvec = - dvec
  }
  
  ### Step 2: Set the constraints
  # Constraints (Set with negative sign because they are in <= form but the 
  # default is in >= form).
  Amat = A
  bvec = rhs
  # Change the lower bounds to inequality constriants
  lb_Amat = diag(length(lb))
  lb_bvec = lb
  # Update constraint matrices
  Amat = rbind(Amat, lb_Amat)
  bvec = c(bvec, lb_bvec)
  # Number of equality constraints
  #meq = min(nrow(A), 2)
  # upper bound
  uvec = c(rhs, rep(Inf, length(lb_bvec)))
  # solution = solve.QP(Dmat,dvec,t(Amat),bvec=bvec)
  settings = osqpSettings(verbose = FALSE, 
                          max_iter = 1e5, 
                          eps_abs=1e-8, 
                          eps_rel = 1e-8)
  # Solve the quadratic program
  osqp_solution = solve_osqp(Dmat, dvec, Amat, l=bvec, u=uvec, pars=settings)
  # Compute the real solution because obj0 was not included in the objective
  # function
  osqp_solution_real = osqp_solution$info$obj_val + obj0
  return(list(objval = osqp_solution_real,
              x = osqp_solution$x))
}

#' Computes the test statistics via bootstrapping
#'
#' @description This function computes the test statistics via bootstrapping.
#'
#' @require modelr
#'
#' @param J The number of distinct nonzero values in vector \eqn{\bm{y}}.
#' @param s_star The value of \eqn{s^\ast} in the cone-tightening procedure.
#' @param beat_obs_hat The value of \eqn{\hat{\beta}_{\mathrm{obs}}} based on
#'    the function supplied by the user.
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @returns Returns the list of test statistics obtained from bootstrap, i.e.
#'    \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.
#'
#' @export
beta_bs <- function(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt, func_obs, 
                    beta_obs_hat, beta_tgt, tau, N, f_solver){
  T_bs = NULL
  # Loop through all indices in the bootstrap
  for (i in 1:bs_num){
    #### Step 1: Set the seed
    set.seed(bs_seed + i)
    ####  Step 2: Draw the subsample
    df_bs = as.data.frame(resample_bootstrap(as.data.frame(df)))
    # Re-index the rows
    rownames(df_bs) = 1:nrow(df_bs)
    ####  Step 3: Compute the bootstrap estimates
    # Compute the value of beta_bs_star using the function func_obs
    beta_bs_star = func_obs(df_bs)
    ####  Step 4: Compute the bootstrap test statistic
    beta_bs_bar = beta_bs_star - beta_obs_hat + s_star
    T_bs_i = prog_cone(A_obs, A_tgt, beta_bs_bar, beta_tgt, tau, "bootstrap", 
                       N, f_solver)$objval
    T_bs = c(T_bs, T_bs_i)
  }
  # Return the general solution
  return(T_bs)
}

#' Auxiliary function to calculate the p-value
#'
#' @description This function computes the \eqn{p}-value of the test based on the
#'    bootstrap estimates.
#'
#' @param T_bs The test statistics obtained from bootstrap.
#' @param T_star The test statistics obtained from quadratic program (5).
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#'
#' @returns Returns the \eqn{p}-value.
#'
#' @export
p_eval <- function(T_bs, T_star, p_sig){
  # Initialization
  p_val = NULL
  decision = 1 # decision = 1: rejected, decision = 0: not rejected
  alpha = 0
  while (decision != 0 & alpha <= 1){
    T_quan = as.numeric(quantile(T_bs, probs=c(1 - alpha)))
    if (T_star >= T_quan){
      decision = 0
    }
    alpha = alpha + 10^(-p_sig)
  }
  p_val = round(alpha - 10^(-p_sig), digits = p_sig)
  return(p_val)
}

#' Auxiliary function to create the constraints for the linear program for tau
#'
#' @description This function computes the value of tau based on Kamat (2018),
#'   i.e. linear program (6) of Torgovitsky (2019).
#'
#' @param length_tau The number of variables in the constraint.
#' @param coeff_tau The coefficient in front of tau in the constraint.
#' @param coeff_x The coefficient in front of x_i in the constraint.
#' @param ind_x The index of x_i, i.e. the value of i.
#' @param rhs The RHS coefficient of the new constraint.
#' @param sense The equality or inequality symbol to be used in the new
#'    constraint.
#' @param lp_lhs_tau The constraint matrix to be updated.
#' @param lp_rhs_tau The RHS vector of the linear constraints to be updated.
#' @param lp_sense_tau The sense vector fo the linear constraints to be updated
#'
#' @returns Returns the updated matrix of the constraint matrix, RHS vector of
#'    of the linear constraints and the sense vector.
#'
#' @export
tau_constraints <- function(length_tau, coeff_tau, coeff_x, ind_x, rhs, sense,
                            lp_lhs_tau, lp_rhs_tau, lp_sense_tau){
  temp = rep(0, length_tau)
  temp[1] = coeff_tau
  temp[ind_x] = coeff_x
  # Update the lhs, rhs and sense of the constraints
  lp_lhs_tau = rbind(lp_lhs_tau, c(temp))
  lp_rhs_tau = c(lp_rhs_tau, 0)
  lp_sense_tau = c(lp_sense_tau, sense)
  return(list(lp_lhs_tau = lp_lhs_tau,
              lp_rhs_tau = lp_rhs_tau,
              lp_sense_tau = lp_sense_tau))
}

#' Check the input of the user
#'
#' @description This function checks the input of the user. If there is any 
#'    invalid input, the program will sotp the program and generates 
#'    appropriate error messages.
#'
#' @inheritParams dkqs_cone
#'
#' @export
dkqs_cone_check <- function(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                            bs_num, p_sig, tau_input, lpsolver){
  # Check the dataframe provided by the user
  if (class(df) %in% c("data.frame", "matrix") == TRUE){
    df = as.data.frame(df)  
  } else {
    stop(gsub("\\s+", " ",
              "The data povided 'df' must either be a data.frame, data.table, 
               matrix."), call. = FALSE)    
  }
  
  # Check the matrices A_obs and A_tgt provided by the user
  matrix_names = c("A_obs", "A_tgt")
  matrix_list = list(A_obs, A_tgt)
  for (i in 1:2){
    # Check the format of the matrices
    if (class(matrix_list[[i]]) %in% c("data.frame", "matrix") == FALSE){
      stop(gsub("\\s+", " ",
                paste0("The argument '", matrix_names[i], "' must either be a 
                data.frame, data.table, or matrix.")), call. = FALSE)   
    } else{
      # Ensure the variable is in matrix form
      matrix_list[[i]] = as.matrix(matrix_list[[i]])
      # Check whether the matrices are numeric
      if (is.numeric(matrix_list[[i]]) == FALSE){
        stop(paste0("The argument '", matrix_names[i], "' has to be numeric."),
             call. = FALSE)
      } 
    }
  }
  if (dim(A_tgt)[1] != 1){
    stop("The argument 'A_tgt' has to be a column vector", call. = FALSE)
  }
  
  # Update A_obs and A_tgt to ensure that they are both matrices
  A_obs = matrix_list[[1]]
  A_tgt = matrix_list[[2]]

  # Check the function provided by the user
  if (class(func_obs) != "function"){
    stop("The input of 'func_obs' has to be a function.", call. = FALSE)
  } else{
    beta_obs_hat = func_obs(df)
    beta_obs_hat = as.matrix(beta_obs_hat)
    # Check if the output is numeric
    if (is.numeric(beta_obs_hat[,1]) == FALSE){
      stop("The output of 'func_obs' has to be numeric.", call. = FALSE)
    } else{
      if (dim(beta_obs_hat)[2] != 1){
        stop("The output of 'func_obs' has to be a column vector", 
        call. = FALSE)
      } else if (dim(beta_obs_hat)[1] != dim(A_obs)[1]){
        stop("The number of rows in the output of 'func_obs' has to be the 
             same as the number of rows in the argument 'A_obs'.", 
             call. = FALSE)
      }
    }
  }
  
  # Check the seed for bootstrap
  if (!(is.numeric(beta_tgt) == TRUE & length(beta_tgt) == 1)) {
    stop("The argument 'beta_tgt' must be a scalar.", call. = FALSE)
  }  

  # Check the seed for bootstrap
  if (!(is.numeric(bs_seed) == TRUE & length(bs_seed) == 1)) {
    stop("The seed to be used in the bootstrap ('bs_seed') must be a scalar.", 
         call. = FALSE)
  }
  
  # Check the number of bootstraps
  if ((is.numeric(bs_num) == TRUE & length(bs_num) == 1 & bs_num >= 0) 
      == FALSE){
    stop("The number of bootstrap ('bs_num') must be a positive integer.",
         call. = FALSE)
  }
  
  # Check the number of decimal places
  if ((is.numeric(p_sig) == TRUE & length(p_sig) == 1 & p_sig >= 0 & 
       p_sig %%1 == 0) == FALSE){
    stop("The number of decimal places in the p-value ('p_sig') has to be a
         positive integer.", call. = FALSE)
  }
  
  # Check the value of tau_input
  if ((is.numeric(tau_input) == TRUE & length(tau_input) == 1 & 
       tau_input >= 0 & tau_input <= 1) == FALSE){
    stop("The value of tau ('tau_input') has to be in the interval [0,1].", 
         call. = FALSE)
  }
  
  # Check if the user is going to use a linear programming solver that is 
  # supported. If the user does not specify any linear programming solver, the
  # function will assign a linear programming solver if one of the supported
  # solvers is installed.
  if (hasArg(lpsolver) == TRUE){
    lpsolver = tolower(lpsolver)
  }
  # Case 1: If no lpsolver is specified by the user
  if (is.null(lpsolver) == TRUE){
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
      lpsolver = "gurobi"
    } else if (requireNamespace("quadprog", quietly = TRUE) == TRUE){
      lpsolver = "quadprog"
    } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE){
      lpsolver = "cplexAPI"
    } else {
      stop(gsub("\\s+", " ",
                "Please install one of the following packages required for
                estimation:
                gurobi (version 7.5-1 or later);
                cplexAPI (version 1.3.3 or later);
                quadprog (version 1.5.8 or later)."),
           call. = FALSE)
    }
  } else{
  # Case 2: If user specifies a package that is not supported by the function
    if ((lpsolver %in% c("gurobi", "cplexapi", "quadprog")) == FALSE){
      stop(gsub("\\s+", " ",
           paste0("This function is incompatible with linear programming
                  package '", lpsolver, "'. Please install one of the
                  following linear programming packages instead:
                  gurobi (version 7.5-1 or later);
                  cplexAPI (version 1.3.3 or later);
                  quadprog (version 1.5.8 or later).")),
           call. = FALSE)
    }
  }
  
  return(list(df = df, 
              A_tgt = A_tgt,
              A_obs = A_obs,
              beta_obs_hat = beta_obs_hat,
              lpsolver = lpsolver))
}
