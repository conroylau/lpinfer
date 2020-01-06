#' Computes the p-value of a quadratic program
#'
#' @description This module conducts inference in quadratic programs using the 
#'    procedure suggested by Torgovitsky (2019) that incorporates the 
#'    cone-tightening procedure proposed by Deb, Kitamura, Quah and
#'    Stoye (2018).
#'
#' @import slam gurobi car modelr
#'
#' @param df The dataframe that contains the sample data.
#' @param A_obs The "observed matrix" in the quadratic program 
#'    \eqn{A_{\mathrm{obs}}.
#' @param A_tgt The "target matrix" in the quadratic program 
#'    \eqn{A_{\mathrm{tgt}}.
#' @param func_obs The function that generates the required 
#'    \eqn{\hat{\beta}_{\mathrm{obs}}}.
#' @param beta_tgt The value of \eqn{\hat{\beta}_{\mathrm{tgt}}} (i.e. the value 
#'    of \eqn{t} in the missing data problem) in the null hypothesis.
#' @param bs_seed The starting value of the seed in bootstrap.
#' @param bs_num The total number of bootstraps \eqn{B}.
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#' @param tau_input The value of tau chosen by the user.
#' @param lpsolver The name of the linear programming solver package to 
#'    be used to obtain the solutions to the linear programs. The linear 
#'    programming solvers supported by this module are `\code{cplexAPI}', 
#'    `\code{gurobi}' and `\code{lpsolveAPI}'.
#' @param qpsolver The name of the quadratic programming solver package to 
#'    be used to obtain the solutions to the quadratic programs. The quadratic 
#'    programming solvers supported by this module are `\code{cplexAPI}', 
#'    `\code{gurobi}' and `\code{osqp}'.
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
                      bs_num = 100, p_sig = 2, tau_input = .5, lpsolver = NULL,
                      qpsolver = NULL){
  
  #### Step 1: Check the dependencies
  checkupdate = dkqs_cone_check(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                                bs_num, p_sig, tau_input, lpsolver, qpsolver)
  # Update and return the information returned from the function dkqs_cone_check
  # (a) Dataframe
  df = checkupdate$df
  # (b) Matrices for linear or quadratic programs
  A_obs = checkupdate$A_obs
  A_tgt = checkupdate$A_tgt
  # (c) Solver for linear or quadratic programs
  lpsolver = checkupdate$lpsolver
  qpsolver = checkupdate$qpsolver
  # Display the lpsolver being used
  cat(paste("Linear programming solver used: ", lpsolver, ".\n", sep = ""))
  cat(paste("Quadratic programming solver used: ", qpsolver, ".\n", sep = ""))
  
  #### Step 2: Initialization
  # Initialization
  N = nrow(df)
  J = length(unique(df[,"Y"])) - 1
  # Compute beta_obs_hat using the function defined by user
  beta_obs_hat = checkupdate$beta_obs_hat
  ### Assign the solver used
  # Linear programming
  if (lpsolver == "gurobi"){
    lp_solver = gurobi_optim
  } else if (lpsolver == "cplexapi"){
    lp_solver = cplexapi_optim
  } else if(lpsolver == "lpsolveapi"){
    lp_solver = lpprog_optim
  }
  # Quadratic programming
  if (qpsolver == "gurobi"){
    qp_solver = gurobi_optim
  } else if (qpsolver == "cplexapi"){
    qp_solver = cplexapi_optim
  } else if(qpsolver == "osqp"){
    qp_solver = osqp_optim
  }

  #### Step 3: Choose the value of tau
  tau_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "tau", N, 
                         lp_solver, qp_solver)$objval
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
                          lp_solver, qp_solver)
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
                 beta_obs_hat, beta_tgt, tau, N, lp_solver, qp_solver)
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
#' @description This function formulates the matrices and vectors, and solves 
#'    the quadratic programs (4) or linear programs (5) and (6) in Torgovitsky 
#'    (2019).
#'
#' @param tau The value of tau to be used in the linear program.
#' @param problem The problem that the function will be solved.
#' @param lp_solver Name of the solver that solves the linear programs.
#' @param qp_solver Name of the solver that solves the quadratic programs.
#' @param beta_obs_hat The value of \eqn{\hat{\beta}_{\mathrm{obs}}} based on
#'    the function supplied by the user.
#' @param n The number of observations in the dataframe.
#' @inheritParams dkqs_cone
#'
#' @returns Returns the solution to the quadratic program.
#' 
#' @details The argument \code{problem} must be one of the followings:
#' \itemize{
#'   \item{\code{test} --- this computes the quadratic program for the test 
#'   statistic, i.e. quadratic program (4) in Torgovitsky (2019)}
#'   \item{\code{bootstrap} --- this computes the quadratic program for the 
#'   bootstrap test statistics, i.e. quadratic program (5) in 
#'   Torgovitsky (2019)}
#'   \item{\code{tau} --- this computes the value of tau based on the procedure
#'   suggested by Kamat (2018), i.e. linear program (6) in Torgovitsky (2019)}
#' }
#'
#' @export
prog_cone <- function(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, problem, n,
                      lp_solver, qp_solver){
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
  theta_up = lp_solver(NULL, A_tgt, NULL, ones, c(1), "=", "max", lb)
  theta_down = lp_solver(NULL, A_tgt, NULL, ones, c(1), "=", "min", lb)
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
  # RHS vector and sense for the linear or quadratic program
  lp_rhs = c(beta_tgt, 1)
  lp_sense = c("=", "=")

  #### Step 3: Solve the QP
  # If problem == "test", this function solves linear program (4)
  # If problem == "bootstrap", this function solves linear program (5)
  # If program == "tau", this function solves linear program (6)
  if (problem == "test"){
    ans = qp_solver(obj2, obj1, obj0, rbind(A_tgt, ones), lp_rhs, lp_sense,
                       "min", lb)
  } else if (problem == "bootstrap"){
  # Update lb
    lb_new = lb
    lb_new[ind_up] = rhs_up
    lb_new[ind_down] = rhs_down
    lb_new[ind_0] = rhs_0
    ### Use QP solver to find the optimum
    ans = qp_solver(obj2, obj1, obj0, rbind(A_tgt, ones), lp_rhs, lp_sense,
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
    ans = lp_solver(NULL, c(1,rep(0, cn)), NULL, lp_lhs_tau, lp_rhs_tau,
                  lp_sense_tau, "max", lb_tau)
  }
  return(ans)
}

#' Gurobi solver for quadratic and linear programs
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the `\code{Gurobi}' package.
#'    
#' @require gurobi
#'
#' @param obj2 The matrix that corresponds to the second-order term in the
#'    quadratic program. This input will be set to \code{NULL} if the problem
#'    considered is a linear program.
#' @param obj1 The vector that corresponds to the first-order term in the
#'    quadratic or linear program.
#' @param obj0 The constant term in the quadratic or linear program.
#' @param A The constraint matrix.
#' @param rhs The RHS vector for the linear constraints.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#'
#' @returns Returns the optimal objective value and the corresponding argument
#'   to the quadratic or linear program.
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
  # Result of the linear or quadratic program, and return result
  params = list(OutputFlag=0)
  gurobi_result = gurobi(model, params)
  return(gurobi_result)
}

#' cplexAPI solver for quadratic and linear programs
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the `\code{cplexAPI}' package.
#'    
#' @inheritParams gurobi_optim
#'
#' @returns Returns the optimal objective value and the corresponding argument
#'   to the quadratic or linear program.
#'
#' @export
cplexapi_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  ### Step 1: Translate the notations that are used for gurobi into those to
  ### be used in cplexAPI
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

#' lpSolveAPI solver for linear programs
#'
#' @description This function computes the solution to the linear program
#'    using the `\code{lpsolveAPI}' package.
#'    
#' @require lpSolveAPI
#'
#' @inheritParams gurobi_optim
#'
#' @returns Returns the optimal objective value and the corresponding argument
#'   to the linear program.
#'   
#' @details The package `\code{lpSolveAPI}' cannot be used to solve quadratic 
#'   programs. Hence, the argument \code{obj2} is not used in the function.
#'
#' @export
lpprog_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  ### Step 1: Update the constraint matrices
  # Change the lower bounds to inequality constriants
  lb_Amat = diag(length(lb))
  lb_bvec = lb
  # Update constraint matrices
  A = rbind(A, lb_Amat)
  rhs = c(rhs, lb_bvec)
  sense = c(sense, rep(">=", length(lb_bvec)))
  
  ### Step 2: Basic set-ups of LP
  # Constraints numbers
  lprec = make.lp(nrow = nrow(A), ncol = ncol(A))
  # Model sense
  lp.control(lprec, sense=modelsense)
  # Types of decision variables
  set.type(lprec, 1:ncol(A), type=c("real"))
  set.objfn(lprec, obj1)
  
  ### Step 3: Define the constraints
  for (i in 1:nrow(A)){
    add.constraint(lprec, A[i, ], sense[i], rhs[i])
    
  }
  
  ### Step 4: Solve and obtain solution of LP
  solve(lprec)
  x = get.variables(lprec)
  objval = get.objective(lprec)
  return(list(x = x,
              objval = objval))
}

#' osqp solver for quadratic programs
#'
#' @description This function computes the solution to the quadratic program
#'    using the `\code{osqp}' package.
#'    
#' @require osqp
#'
#' @inheritParams gurobi_optim
#'
#' @returns Returns the optimal objective value and the corresponding argument
#'   to the linear program.
#'   
#' @details The package `\code{osqp}' cannot be used to solve linear programs. 
#'   Hence, the argument \code{obj2} cannot be \code{NULL}.
#'
#' @export
osqp_optim <- function(obj2, obj1, obj0, A, rhs, sense, modelsense, lb){
  ### Step 1: Set the objective function
  # Set objective function
  # It is multiplied by 2 because of the term 1/2 in the quadratic term in the 
  # osqp package.
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

#' Computes the bootstrap test statistics
#'
#' @description This function computes the bootstrap test statistics.
#'
#' @require modelr
#'
#' @param J The number of distinct nonzero values in vector \eqn{\bm{y}}.
#' @param s_star The value of 
#'    \eqn{\hat{s}^\ast \equiv A_{\mathrm{obs}}\hat{\bm{x}}_n^\ast} 
#'    in the cone-tightening procedure.
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @returns Returns the list of bootstrap test statistics, i.e.
#'    \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.
#'
#' @export
beta_bs <- function(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt, func_obs, 
                    beta_obs_hat, beta_tgt, tau, n, lp_solver, qp_solver){
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
                       n, lp_solver, qp_solver)$objval
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
#' @description This function generates the constraints in the linear program
#'   for computing the value of tau based on the procedure suggested by Kamat 
#'   (2018), i.e. linear program (6) of Torgovitsky (2019).
#'
#' @param length_tau The number of variables in the constraint.
#' @param coeff_tau The coefficient in front of tau in the constraint.
#' @param coeff_x The coefficient in front of \eqn{x_i} in the constraint.
#' @param ind_x The index of \eqn{x_i}, i.e. the value of \eqn{i}.
#' @param rhs The RHS vector of the new constraint.
#' @param sense The equality or inequality symbol to be used in the new
#'    constraint.
#' @param lp_lhs_tau The constraint matrix to be updated.
#' @param lp_rhs_tau The RHS vector of the linear constraints to be updated.
#' @param lp_sense_tau The sense vector fo the linear constraints to be updated
#'
#' @returns Returns the list of updated constraint matrix, RHS vector of the 
#'   linear constraints and the sense vector.
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

#' Checks and update the input
#'
#' @description This function checks and updates the input of the user. If 
#'    there is any invalid input, this function will be terminated and generates 
#'    appropriate error messages.
#'
#' @inheritParams dkqs_cone
#'
#' @export
dkqs_cone_check <- function(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                            bs_num, p_sig, tau_input, lpsolver, qpsolver){
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
  
  # Check if the user supplied a name of linear or quadratic programming solver
  # that is supported by the function. If the user does not specify any linear
  # programming solver, the function will assign a linear or quadratic
  # programming solver that is supported. 
  solvers = list(lpsolver, qpsolver)
  solver_incomp = c(FALSE, FALSE)
  # Package installation message generator
  gurobi_msg = "gurobi (version 8.1-1 or later)"
  cplexapi_msg = "cplexAPI (version 1.3.3 or later)"
  osqp_msg = "osqp (version 0.6.0.3 or later)"
  lpsolveapi_msg = "lpSolveAPI (version 5.5.2.0 or later)"
  # In the for-loop, i = 1 refers to the linear program and i = 2 refers to the 
  # quadratic program
  for (i in 1:2){
    if (hasArg(solvers[[i]]) == TRUE){
      solvers[[i]] = tolower(hasArg(solvers[[i]]))
    }
    # Case 1: If no lpsolver is specified by the user
    if (is.null(solvers[[i]]) == TRUE){
      if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
        solvers[[i]] = "gurobi"
      } else if (requireNamespace("lpsolveapi", quietly = TRUE) == TRUE & 
                 i == 1){
        solvers[[i]] = "lpsolveapi"
      } else if (requireNamespace("osqp", quietly = TRUE) == TRUE & i == 2){
        solvers[[i]] = "osqp"
      } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE){
        solvers[[i]] = "cplexapi"
      } else {
        stop(gsub("\\s+", " ",
                  paste0("Please install one of the following packages required
                         for estimation: ",
                         gurobi_msg, "; ",
                         cplexapi_msg, "; ",
                         osqp_msg, " AND ", lpsolveapi_msg, ".")),
             call. = FALSE)
      }
    } else{
    # Case 2: If user specifies a package that is not supported by the function
      if ((lpsolver %in% c("gurobi", "cplexapi", "osqp", "lpsolveapi")) 
          == FALSE){
      solver_incomp[i] = TRUE
      }
    }
  }
  # Display error message if user entered a linear or quadratic program solver 
  # that is not supported by the function
  lp_solver_incomp = NULL
  qp_solver_incomp = NULL
  connector_incomp = NULL
  if (sum(solver_incomp) != 0){
    if (solver_incomp[1] == TRUE){
      lp_solver_incomp = paste0("linear programming solver package '", 
                                solvers[1], "' ") 
    }
    if (solver_incomp[2] == TRUE){
      qp_solver_incomp = paste0("quadratic programming solver package '", 
                                solvers[2], "' ") 
    }
    if (sum(solver_incomp) == 2){
      connector_incomp = " and "
    }
    stop(gsub("\\s+", " ",
              paste0("This function is incompatible with ", lp_solver_incomp, 
                    connector_incomp, qp_solver_incomp,
                    ". Please install one of the following packages instead: ",
                    gurobi_msg, "; ",
                    cplexapi_msg, "; ",
                    osqp_msg[solver_incomp[1]], connector_incomp, 
                    lpsolveapi_msg[solver_incomp[2]], ".")),
         call. = FALSE)
  }
  
  # Case 3 - If user specified a linear (resp. quadratic) programming solver for 
  # a quadratic (resp linear) programming solver.
  if (solvers[[1]] == "osqp"){
    stop(gsub("\\s+", " ",
              paste0("The package 'osqp' cannot be used to solve linear programs
                     Please use one of the following packges instead: ",
                     gurobi_msg, "; ",
                     cplexapi_msg, "; ",
                     lpsolveapi_msg, ".")),
         call. = FALSE)
  }
  if (solvers[[2]] == "lpsolveapi"){
    stop(gsub("\\s+", " ",
              paste0("The package 'lpsolveAPI' cannot be used to solve 
                     quadratic programs. Please use one of the following packges 
                     instead: ",
                     gurobi_msg, "; ",
                     cplexapi_msg, "; ",
                     osqp_msg, ".")),
         call. = FALSE)
  }

  # Step 4 - Update the solver names
  lpsolver = solvers[[1]]
  qpsolver = solvers[[2]]
  
  return(list(df = df, 
              A_tgt = A_tgt,
              A_obs = A_obs,
              beta_obs_hat = beta_obs_hat,
              lpsolver = lpsolver,
              qpsolver = qpsolver))
}
