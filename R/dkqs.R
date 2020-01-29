#' Computes the p-value of a quadratic program
#'
#' @description This module conducts inference in quadratic programs using the 
#'    procedure suggested by Torgovitsky (2019) that incorporates the 
#'    cone-tightening procedure proposed by Deb, Kitamura, Quah and
#'    Stoye (2018).
#' @import gurobi cplexAPI Rcplex Momocs limSolve
#'
#' @param df The data being used in the inference.
#' @param A_obs The "observed matrix" in the inference \eqn{A_{\mathrm{obs}}.
#' @param A_tgt The "target matrix" in the inference \eqn{A_{\mathrm{tgt}}.
#' @param func_obs The function that generates the required 
#'    \eqn{\hat{\beta}_{\mathrm{obs}}}.
#' @param beta_tgt The value of \eqn{\hat{\beta}_{\mathrm{tgt}}} (i.e. the 
#'    value of \eqn{t} in the missing data problem) in the null hypothesis.
#' @param bs_seed The starting value of the seed in bootstrap.
#' @param bs_num The total number of bootstraps \eqn{B} to be conducted.
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#' @param tau_input The value of tau chosen by the user.
#' @param solver The name of the linear and quadratic programming solver that 
#'    are used to obtain the solution to linear and quadratic programs. 
#'    The solvers supported by this module are `\code{cplexAPI}', 
#'    `\code{gurobi}', `\code{limSolve}' and `\code{Rcplex}'.
#' @param noisy The boolean variable for whether the result messages should
#'    be displayed in the inference procedure. If it is set as \code{TRUE}, the 
#'    messages are displayed throughout the procedure. Otherwise, the messages
#'    will not be displayed.
#'    
#' @return Returns a list of output calculated from the function:
#'   \item{p_val}{\eqn{p}-value.}
#'   \item{tau}{The value of tau used \eqn{\tau^\ast} in the linear and 
#'      quadratic programs.}
#'   \item{T_n}{Test statistic \eqn{T_n}.}
#'   \item{T_bs}{The list of bootstrap test statistics 
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'   \item{beta_bs_bar}{The list of \eqn{\tau}-tightened re-centered bootstrap 
#'      estimators \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'   \item{lb0}{Logical lower bound of the problem.}
#'   \item{ub0}{Logical upper bound of the problem.}
#'      
#' @export
dkqs_cone <- function(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed = 1,
                      bs_num = 100, p_sig = 2, tau_input = .5, solver = NULL,
                      noisy = TRUE){
  
  #### Step 1: Check and update the dependencies
  checkupdate = dkqs_cone_check(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                                bs_num, p_sig, tau_input, solver, noisy)
  # Update and return the quantities returned from the function dkqs_cone_check
  # (a) Dataframe
  df = checkupdate$df
  # (b) Matrices for linear and quadratic programs
  A_obs = checkupdate$A_obs
  A_tgt = checkupdate$A_tgt
  # (c) Solver for linear and quadratic programss
  solver = checkupdate$solver
  # Display the solver used
  if (noisy == TRUE){
    cat(paste("Linear and quadratic programming solver used: ", solver, ".\n", 
              sep = ""))    
  }
  
  #### Step 2: Initialization
  # Initialization
  n = nrow(df)
  J = length(unique(df[,"Y"])) - 1
  # Compute beta_obs_hat using the function defined by user
  beta_obs_hat = checkupdate$beta_obs_hat
  ### Assign the solver to be used
  if (solver == "gurobi"){
    solver = gurobi_optim
  } else if (solver == "cplexapi"){
    solver = cplexapi_optim
  } else if (solver == "rcplex"){
    solver = rcplex_optim
  } else if (solver == "limsolve"){
    solver = limsolve_optim
  }

  #### Step 3: Choose the value of tau
  tau_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "tau", n,
                         solver)
  if (tau_input > tau_return$objval){
    tau = tau_return$objval
  } else if (tau_input <= tau_return$objval){
    tau = tau_input
  } else {
    # Error message when the problem is infeasible.
    stop("The problem is infeasible. Choose other values of tau.")
  }
  
  #### Step 4: Compute T_n, x_star and s_star
  # Compute T_n
  T_n = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "test", n, 
                          solver)$objval
  # Return and stop program if it is infeasible
  if (is.null(T_n) == TRUE){
    stop("The problem is infeasible. Choose other values of tau.")
  }
  # Compute s_star
  x_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "cone", n, 
                     solver)
  x_star = x_return$x
  s_star = A_obs %*% x_star
  
  #### Step 5: Compute the bootstrap estimates
  # T_bs is the list of bootstrap test statistics used
  T_bs_return = beta_bs(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt, func_obs, 
                 beta_obs_hat, beta_tgt, tau, n, solver)
  T_bs = T_bs_return$T_bs
  beta_bs_bar = T_bs_return$beta_bs_bar_set
  
  #### Step 6: Compute the p-value
  # decision = 1 refers to rejected, decision = 0 refers to not rejected
  p_val = p_eval(T_bs, T_n, p_sig)
  
  #### Step 7: Obtain logical bounds for the function qrci
  lb0 = x_return$lb0
  ub0 = x_return$ub0
  
  #### Step 8: Print results
  if (noisy == TRUE){  
    cat(paste("-----------------------------------", "\n"))
    cat(paste("Test statistic: ", round(T_n, digits = 5), ".\n", sep = ""))
    cat(paste("p-value: ", p_val, ".\n", sep = ""))
    cat(paste("Value of tau used: ", round(tau, digits = 5), ".\n", 
              sep = ""))
  }
  
  invisible(list(p_val = as.numeric(p_val), 
                 tau = as.numeric(tau), 
                 T_n = as.numeric(T_n), 
                 T_bs = T_bs,
                 beta_bs_bar = beta_bs_bar,
                 lb0 = lb0$objval,
                 ub0 = ub0$objval))
}

#' Formulating and solving quadratic programs
#'
#' @description This function formulates the matrices and vectors, and solves 
#'    the quadratic programs (4) or linear programs (5) and (6) in Torgovitsky 
#'    (2019).
#'
#' @param tau The value of tau to be used in the linear program.
#' @param problem The problem that the function will be solved.
#' @param solver Name of the solver that solves the linear and quadratic 
#'    programs.
#' @param beta_obs_hat The value of \eqn{\hat{\beta}_{\mathrm{obs}}} based on
#'    the function supplied by the user.
#' @param n The number of observations in the dataframe.
#' @inheritParams dkqs_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
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
                      solver){
  #### Step 1: Obtain dimension of A_tgt
  rn = dim(A_tgt)[1]
  cn = dim(A_tgt)[2]

  #### Step 2: Formulation of constraints
  ones = matrix(rep(1, cn), nrow = 1)
  lb = matrix(rep(0, cn), nrow = 1)
  # Obtain parameters
  theta_down = solver(NULL, A_tgt, n, ones, c(1), "=", "min", lb)
  theta_up = solver(NULL, A_tgt, n, ones, c(1), "=", "max", lb)
    
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
  # If problem == "cone", this function solves linear program (5)
  # If program == "tau", this function solves linear program (6)
  if (problem == "test"){
    ans = solver(A_obs, beta_obs_hat, n, rbind(A_tgt, ones), lp_rhs, lp_sense,
                 "min", lb)
  } else if (problem == "cone"){
  # Update lb
    lb_new = lb
    lb_new[ind_up] = rhs_up
    lb_new[ind_down] = rhs_down
    lb_new[ind_0] = rhs_0
    ### Find solution using the solver
    ans = solver(A_obs, beta_obs_hat, n, rbind(A_tgt, ones), lp_rhs, lp_sense,
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
    ans = solver(NULL, c(1, rep(0, cn)), n, lp_lhs_tau, lp_rhs_tau,
                 lp_sense_tau, "max", lb_tau)
  }
  
  #### Step 4: Append the logical bounds to the results
  ans$lb0 = theta_down
  ans$ub0 = theta_up
  return(ans)
}

#' LP and QP solver by Gurobi
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the `\code{Gurobi}' package.
#'    
#' @import gurobi
#'
#' @param Af The matrix that is involved in the objective function.
#' @param bf The vector that is involved in the objective function.
#' @param nf The number of observations in the dataframe.
#' @param A The constraint matrix.
#' @param rhs The RHS vector for the linear constraints.
#' @param modelsense The indicator of whether the model is to max or min an
#'    objective function.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
gurobi_optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  ### Step 1: Obtain the coefficients of the objective function
  objective_return = objective_function(Af, bf, nf)
  
  ### Step 2: Gurobi set-up
  model = list()
  model$Q = objective_return$obj2
  model$obj = objective_return$obj1
  model$objcon = objective_return$obj0
  model$A = A
  model$rhs = rhs
  model$sense = sense 
  model$modelsense = modelsense
  model$lb = lb
  
  ### Step 3: Result of the linear or quadratic program, and return result
  params = list(OutputFlag=0)
  solution = gurobi(model, params)
  return(list(objval = as.numeric(solution$objval),
              x = as.numeric(solution$x)))
}

#' LP and QP solver by cplexAPI
#'
#' @description This function computes the solution to the quadratic and linear
#'    programs using the `\code{cplexAPI}' package.
#'    
#' @inheritParams cplexAPI
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
cplexapi_optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  ### Step 1: Obtain the coefficients of the objective function
  objective_return = objective_function(Af, bf, nf)
  
  ### Step 2: Update the notations
  # Model sense
  modelsense[modelsense == "min"] = CPX_MIN
  modelsense[modelsense == "max"] = CPX_MAX
  # Inequality/equality signs
  sense[sense == "<="] = "L"
  sense[sense == ">="] = "G"
  sense[sense == "=="] = "E"
  sense[sense == "="] = "E"
  # Bounds
  lb[lb == Inf] = CPX_INFBOUND
  lb[lb == -Inf] = -CPX_INFBOUND
  ub = rep(CPX_INFBOUND, length(lb))
  
  ### Step 3: cplexAPI environment
  # Model environment
  env = cplexAPI::openEnvCPLEX()
  cplexAPI::setDblParmCPLEX(env, 1016, 1e-06)
  prob = cplexAPI::initProbCPLEX(env)
  cplexAPI::chgProbNameCPLEX(env, prob, "sample")
  
  # Constraint matrices
  cnt = apply(A, MARGIN = 2, function(x) length(which(x != 0)))
  beg = rep(0, ncol(A))
  beg[-1] = cumsum(cnt[-length(cnt)])
  ind = unlist(apply(A, MARGIN = 2, function(x) which(x != 0) - 1))
  val = c(A)
  val = val[val != 0]
  
  ### Step 4: Solve the problem
  # A linear program is identified if obj2 == NULL
  if (is.null(obj2) == TRUE){
    # Solving linear program
    copyLpwNamesCPLEX(env,
                      prob,
                      ncol(A),
                      nrow(A),
                      modelsense,
                      objective_return$obj1,
                      rhs,
                      sense,
                      beg,
                      cnt,
                      ind,
                      val,
                      lb,
                      ub)
    lpoptCPLEX(env, prob)
    solution = solutionCPLEX(env, prob)
  } else {
    # Solving quadratic program
    stop("This version can only solve linear programs by CPLEX at the moment. 
         Please use another solver for quadratic progarms.")
  }
  return(list(objval = as.numeric(solution$objval),
              x = as.numeric(solution$x)))
}

#' LP and QP solver by Rcplex
#'
#' @description This function computes the solution to the linear and quadratic
#'    programs using the `\code{Rcplex}' package.
#'    
#' @import Rcplex
#'
#' @inheritParams gurobi_optim
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
rcplex_optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  ### Step 1: Obtain the coefficients of the objective function
  objective_return = objective_function(Af, bf, nf)
  
  ### Step 2: Update vectors and sense
  # Update sense
  sense[sense == ">="] = "G"
  sense[sense == "<="] = "L"
  sense[sense == "="]  = "E"
  # Define upper bound
  ub = rep(Inf, length(lb))
  # Define Q matrix
  # - Keep Qmat as NULL for linear program
  # - Multiply obj2 by 2 for Q for quadratic program to offset the 1/2 factor
  if (is.null(objective_return$obj2) == TRUE){
    Qmat = objective_return$obj2
  } else {
    Qmat = 2*objective_return$obj2
  }
  
  ### Step 3: Solve model
  solution = Rcplex(cvec = t(objective_return$obj1),
                    Amat = A, 
                    bvec = rhs,
                    Qmat = Qmat,
                    lb = lb,
                    sense = sense,
                    ub = ub,
                    objsense = modelsense,
                    vtype = "C",
                    n = 1)
  
  ### Step 3: Update and return result
  if (is.null(objective_return$obj0) == FALSE){
    objval = solution$obj + objective_return$obj0
  } else {
    objval = solution$obj
  }
  return(list(objval = as.numeric(objval),
              x = as.numeric(solution$xopt)))
}

#' LP and QP solver by limSolve
#' 
#' @description This function computes the solution to linear and quadratic 
#'    programs using the `\code{limSolve}' package.
#' 
#' @import limSolve
#'
#' @inheritParams gurobi_optim
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'   
#' @export
limsolve_optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  ### Step 1: Obtain the coefficients of the objective function
  objective_return = objective_function(Af, bf, nf)
  
  ### Step 2: Update lower bounds
  # Change the lower bounds to inequality constriants
  lb_Amat = diag(length(lb))
  lb_bvec = lb
  # Update constraint matrices
  A = rbind(A, lb_Amat)
  rhs = c(rhs, lb_bvec)
  sense = c(sense, rep(">=", length(lb_bvec)))
  
  ### Step 3: Update constraints
  # Objective function
  if (modelsense == "max"){
    fcost = - objective_return$obj1 
  } else if (modelsense == "min"){
    fcost = objective_return$obj1
  }
  # Equality constraints
  Emat = A[sense == "=",]
  Fvec = rhs[sense == "="]
  # Inequality constraint >=
  Gmat1 = A[sense == ">=",]
  Hvec1 = rhs[sense == ">="]   
  # Inequality constraint <= 
  Gmat2 = -A[sense == "<=",]
  Hvec2 = -rhs[sense == "<="]
  # Combine G and h matrices
  Gmat = rbind(Gmat1, Gmat2)
  Hvec = as.matrix(c(c(Hvec1), c(Hvec2)), ncol = 1, byrow = TRUE)
  
  ### Step 4 - Solve the model
  # Linear solver is used if obj2 is a zero matrix (i.e. number of zeros equals 
  # the total number of elements) or NULL
  if (is.null(objective_return$obj2) == TRUE | 
      sum(objective_return$obj2 == 0) == length(objective_return$obj2)){
    ### Linear program solver
    solution = linp(E = Emat, F = Fvec, G = Gmat, H = Hvec, Cost = fcost)
    # Obtain objective function, and add back the constant term, negate the 
    # solution if it is a max problem
    if (modelsense == "max"){
      objval = -solution$solutionNorm + objective_return$obj0      
    } else if (modelsense == "min"){
      objval = solution$solutionNorm + objective_return$obj0      
    }
  } else {
    if (modelsense == "min"){
      ### Quadratic program solver
      # Formulate the two matrices
      Amat = Af * sqrt(nf)
      Bvec = bf * sqrt(nf)
      solution = lsei(A = Amat, B = Bvec, E = Emat, F = Fvec, G = Gmat, 
                      H = Hvec)
      # Obtain objective function
      objval = solution$solutionNorm
    } else if (modelsense == "max"){
      stop("This package cannot be used to solve max problems.")
    }
  }
  # Optimal the optimal value of x
  x = solution$X
  return(list(x = as.numeric(x),
              objval = as.numeric(objval)))
}


#' Auxiliary function to return the coefficient terms of the objective 
#' functions
#' 
#' @description This function computes the matrics in the objective functions 
#'    for linear programs. This function takes matrix \eqn{\bm{A}} and 
#'    \eqn{\bm{\beta}} as input and computes the coefficients of the objective 
#'    function. 
#'    
#' @param A The matrix \eqn{\bm{A}},
#' @param beta The column vector \eqn{\bm{\beta}}.
#' @param n The sample size \eqn{n}.
#' 
#' @details 
#' \itemize{
#'   \item{\strong{Quadratic programs} ---
#'      Given inputs \eqn{\bm{A} \in \mathbf{R}^{m\times m}} and 
#'      \eqn{\bm{b} \in \mathbf{R}^m}, the equation of the objective function 
#'      of the quadratic program can be written as
#'      \deqn{n (\bm{A}\bm{x} -\bm{\beta})'(\bm{A}\bm{x}-\bm{\beta})
#'      = n\bm{x}'\bm{A}'\bm{A}\bm{x} - 2n\bm{A}'\bm{\beta}\bm{x} + 
#'      n\bm{\beta}'\bm{\beta}.}}
#'   \item{\strong{Linear programs} ---
#'      For all linear problems that are considered in this code, one of 
#'      \eqn{\bm{A}} and \eqn{\bm{b}} is \code{NULL} or is a zero vector. The 
#'      term that is nonzero and nonnull will be used as \code{obj1}.}
#' }
#' 
#' @return Returns the following three quantities: \code{obj2} is the 
#'   coefficient of the quadratic term, \code{obj1} is the coefficient of the 
#'   linear term and \code{obj0} is the constant term. More explicitly, their
#'   form are given as follows:
#'   \item{obj2}{This is the coefficient for the second-order term. It is 
#'     returned as \code{NULL} for linear programs and \eqn{n\bm{A}'\bm{A}} for
#'     quadratic programs.}
#'   \item{obj1}{This is the coefficient term of the linear term. For quadratic 
#'     programs, it is returned as \eqn{-2n\bm{A}'\bm{\beta}}.}
#'   \item{obj0}{This is the constant term of the linear program. For quadratic
#'     programs, it is returned as \eqn{n\bm{\beta}'\bm{\beta}}.}
#' @export
objective_function <- function(A, b, n){
  # If-else function to determine if it corresponds to a linear or quadratic
  # program. This is identified by whether one of A and b is null or nonzero
  # because it would be the case for linear programs that are considered in
  # this package.
  if (is.null(A) == TRUE | sum(A == 0) == length(A)){
    # Linear program coefficients with nonzero vector b
    obj2 = NULL
    obj1 = b
    obj0 = 0
  } else if (is.null(b) == TRUE | sum(b == 0) == length(b)){
    # Linear program coefficients with nonzero matrix A
    obj2 = NULL
    obj1 = A
    obj0 = 0
  } else {
    # Quadratic program coefficients
    obj2 = t(A) %*% A * n
    obj1 = -2 * t(A) %*% b * n
    obj0 = t(b) %*% b * n
  }
  # Return the above objective functions
  return(list(obj2 = obj2, obj1 = obj1, obj0 = obj0))
}

#' Computes the bootstrap test statistics
#'
#' @description This function computes the bootstrap test statistics.
#'
#' @import modelr
#'
#' @param J The number of distinct nonzero values in vector \eqn{\bm{y}}.
#' @param s_star The value of 
#'    \eqn{\hat{s}^\ast \equiv A_{\mathrm{obs}}\hat{\bm{x}}_n^\ast} 
#'    in the cone-tightening procedure.
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @return Returns the list of estimates from bootstrap:
#'   \item{T_bs}{A list of bootstrap test statistics 
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'  \item{beta_bs_bar_set}{A list of \eqn{\tau_n}-tightened recentered 
#'     bootstrap estimates \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}}
#'
#' @export
beta_bs <- function(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt, func_obs, 
                    beta_obs_hat, beta_tgt, tau, n, solver){
  T_bs = NULL
  beta_bs_bar_set = NULL
  # Loop through all indices in the bootstrap
  for (i in 1:bs_num){
    #### Step 1: Set the seed
    set.seed(bs_seed + i)
    ####  Step 2: Draw the subsample
    #df_bs = as.data.frame(resample_bootstrap(as.data.frame(df)))
    df_bs = as.data.frame(sample_n(df, n, replace = TRUE))
    # Re-index the rows
    rownames(df_bs) = 1:nrow(df_bs)
    ####  Step 3: Compute the bootstrap estimates
    # Compute the value of beta_bs_star using the function func_obs
    beta_bs_star = func_obs(df_bs)
    ####  Step 4: Compute the bootstrap test statistic
    beta_bs_bar = beta_bs_star - beta_obs_hat + s_star
    T_bs_i = prog_cone(A_obs, A_tgt, beta_bs_bar, beta_tgt, tau, "cone", n, 
                       solver)$objval
    T_bs = c(T_bs, T_bs_i)
    beta_bs_bar_set = cbind(beta_bs_bar_set, beta_bs_bar)
  }
  # Return the bootstrap test statistic
  return(list(T_bs = T_bs,
              beta_bs_bar_set = beta_bs_bar_set))
}

#' Auxiliary function to calculate the p-value
#'
#' @description This function computes the \eqn{p}-value of the test based on
#'    the bootstrap estimates.
#'
#' @param T_bs The test statistics obtained from bootstrap.
#' @param T_n The test statistics obtained from quadratic program (5).
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#'
#' @return Returns the \eqn{p}-value:
#'   \item{p_val}{\eqn{p}-value that is corrected to \code{p_sig} decimal
#'      places.}
#'
#' @export
p_eval <- function(T_bs, T_n, p_sig){
  # Initialization
  p_val = NULL
  decision = 1 # decision = 1: rejected, decision = 0: not rejected
  alpha = 0
  while (decision != 0 & alpha <= 1){
    T_quan = as.numeric(quantile(T_bs, probs=c(1 - alpha)))
    if (T_n >= T_quan){
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
#' @return Returns the list of matrices that corresponds to the updated 
#'   constraints:
#'   \item{lp_lhs_tau}{Upated constraint matrix.}
#'   \item{lp_rhs_tau}{Update RHS vector.}
#'   \item{lp_sense_tau}{Update sense for the constraints.}
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

#' Checks and updates the input
#'
#' @description This function checks and updates the input of the user. If 
#'    there is any invalid input, this function will be terminated and 
#'    generates appropriate error messages.
#'
#' @inheritParams dkqs_cone
#' 
#' @return Returns the list of updated parameters as follows:
#'   \item{df}{Upated data in class \code{data.frame}}
#'   \item{A_obs}{Update "observed" matrix in class \code{matrix}.}
#'   \item{A_tgt}{Update "target" matrix in class \code{matrix}.}
#'   \item{beta_obs_tgt}{Obtain \eqn{\widehat{\bm{\beta}}_{\mathrm{tgt}}} 
#'      that is obtained from the function \code{func_obs}.}
#'   \item{solver}{Update name of solver in lower case.}
#' 
#' @export
dkqs_cone_check <- function(df, A_obs, A_tgt, func_obs, beta_tgt, bs_seed, 
                            bs_num, p_sig, tau_input, solver, noisy){
  ### Part 1. Check the dataframe
  if (class(df) %in% c("data.frame", "matrix") == TRUE){
    df = as.data.frame(df)  
  } else {
    stop(gsub("\\s+", " ",
              "The data povided 'df' must either be a data.frame, a data.table, 
               or a matrix."), call. = FALSE)    
  }
  
  ### Part 2. Check the matrices A_obs and A_tgt
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

  ### Part 3. Check the function
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
  
  ### Part 4. Check beta_tgt
  if (!(is.numeric(beta_tgt) == TRUE & length(beta_tgt) == 1)) {
    stop("The argument 'beta_tgt' must be a scalar.", call. = FALSE)
  }  

  ### Part 5. Check bs_seed
  if (!(is.numeric(bs_seed) == TRUE & length(bs_seed) == 1)) {
    stop("The seed to be used in the bootstrap ('bs_seed') must be a scalar.", 
         call. = FALSE)
  }
  
  ### Part 6. Check bs_num
  if ((is.numeric(bs_num) == TRUE & length(bs_num) == 1 & bs_num >= 0 &
       bs_num%%1 == 0) == FALSE){
    stop("The number of bootstrap ('bs_num') must be a positive integer.",
         call. = FALSE)
  }
  
  ### Part 7. Check p_sig
  if ((is.numeric(p_sig) == TRUE & length(p_sig) == 1 & p_sig >= 0 & 
       p_sig%%1 == 0) == FALSE){
    stop("The number of decimal places in the p-value ('p_sig') has to be a
         positive integer.", call. = FALSE)
  }
  
  ### Part 8. Check tau_input
  if ((is.numeric(tau_input) == TRUE & length(tau_input) == 1 & 
       tau_input >= 0 & tau_input <= 1) == FALSE){
    stop("The value of tau ('tau_input') has to be in the interval [0,1].", 
         call. = FALSE)
  }
  
  ### Part 9. Check solvers
  # Check if the user supplied a name of linear or quadratic programming solver
  # that is supported by the function. If the user does not specify any linear
  # programming solver, the function will assign a linear or quadratic
  # programming solver that is supported. 
  
  # Intialize variable
  solver_incomp = FALSE
  # Package recommendation messages
  gurobi_msg = "gurobi (version 8.1-1 or later)"
  cplexapi_msg = "cplexAPI (version 1.3.3 or later)"
  rcplex_msg = "Rcplex (version 0.3-3 or later)"
  limsolve_msg = "limSolve (version 1.5.6 or later)"
  
  # Case 1: If no solver name is provided by the user
  if (is.null(solver) == TRUE){
    if (requireNamespace("gurobi", quietly = TRUE) == TRUE){
      solver = "gurobi"
    } else if (requireNamespace("limSolve", quietly = TRUE) == TRUE){
      solver = "limsolve"
    } else if (requireNamespace("Rcplex", quietly = TRUE) == TRUE){
      solver = "rcplex"
    } else if (requireNamespace("cplexAPI", quietly = TRUE) == TRUE){
      solver = "cplexapi"
    } else {
      stop(gsub("\\s+", " ",
                paste0("Please install one of the following packages to solve 
                       the linear and quadratic programs: ",
                       gurobi_msg, "; ",
                       cplexapi_msg, "; ",
                       rcplex_msg, "; ",
                       limsolve_msg, ".")),
           call. = FALSE)
    }
    # Display message to indicate that no solver is suggested by the user so
    # the module chooses one for the user
    if (noisy == TRUE){
      cat(paste("No linear and quadratic programming solver is suggested by the
                user. The solver '", solver, "' is automatically selected", 
                ".\n", sep = ""))
    }
  } else{
    # Change solver name to lower cases
    solver = tolower(solver)
  # Case 2: If user specifies a package that is not supported, display error
  # message and suggest appropriate solvers
    if ((solver %in% c("gurobi", "cplexapi", "rcplex", "limsolve"))
        == FALSE){
      stop(gsub("\\s+", " ",
                paste0("This function is incompatible with ", solver, 
                       ". Please install one of the following packages to solve 
                       the linear and quadratic programs: ",
                       gurobi_msg, "; ",
                       cplexapi_msg, "; ",
                       rcplex_msg, "; ",
                       limsolve_msg, ".")),
           call. = FALSE)
    }
  }
  
  ### Step 10. Check noisy
  if (!(noisy == TRUE | noisy == FALSE)){
    stop("The argument 'noisy' has to be boolean.")
  }
  
  ### Step 11. Return the upated information
  return(list(df = df, 
              A_obs = A_obs,
              A_tgt = A_tgt,
              beta_obs_hat = beta_obs_hat,
              solver = solver))
}