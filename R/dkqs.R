#' Compute p-value of a quadratic program
#'
#' @description This function runs the procedure suggested by Torgovitsky (2019)
#'    using the cone-tightening procedure proposed by Deb, Kitamura, Quah and
#'    Stoye (2018).
#'
#' @import slam gurobi car
#'
#' @param df The dataframe that contains the sample data.
#' @param A_obs The observed matrix in the linear program.
#' @param A_tgt The "target matrix" in the linear program.
#' @param func_obs The function that generates the required beta_obs.
#' @param beta_tgt The value of t in the null hypothesis.
#' @param bs_seed The starting value of the seed in bootstrap.
#' @param bs_num The total number of bootstraps.
#' @param p_sig The number of decimal places in the \eqn{p}-value.
#' @param tau_input The value of tau chosen by the user.
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
                      bs_num = 100, p_sig = 2, tau_input = .5){

  #### Step 1: Error checks
  dkqs_cone_errormsg(df, A_obs, A_tgt, beta_tgt, bs_seed, bs_num, 
                     p_sig, tau_input)

  #### Step 2: Initialization
  # Initialization
  N = dim(df)[1]
  J = length(unique(df[,"Y"])) - 1
  # Compute beta_obs_hat using the function defined by user
  beta_obs_hat = func_obs(df)

  #### Step 3: Choose the value of tau
  tau_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, 
                         "tau")$objval
  if (tau_input > tau_return){
    tau = tau_return
  } else if (tau_input <= tau_return){
    tau = tau_input
  } else {
    # Error message when the problem is infeasible.
    stop("The problem is infeasible.")
  }

  #### Step 4: Solve QP (5) in Torgovitsky (2019)
  full_return = prog_cone(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, "T")
  x_star = full_return$x
  s_star = A_obs %*% x_star
  # T_star is the test statistic used
  T_star = full_return$objval

  #### Step 5: Compute the bootstrap estimates
  # T_bs is the list of bootstrap test statistics used
  T_bs = beta_bs(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt,
                 func_obs, beta_obs_hat, beta_bs_bar, beta_tgt, tau, N)
  #### Step 6: Compute the p-value
  # decision = 1 refers to rejected, decision = 0 refers to not rejected
  p_val = p_eval(T_bs, T_star, p_sig)

  cat(paste("The p-value is ", p_val, ".\n", sep = ""))
  cat(paste("The value of tau used is ", tau, ".", sep = ""))
  invisible(list(p = p_val, tau = tau, T_star = T_star, T_bs = T_bs))
}

#' Formulating and solving quadratic programs
#'
#' @description This function formulates and solves the quadratic programs (5)
#'    and (6) in Torgovitsky (2019).
#'
#' @param tau The RHS vector fo the linear constraints.
#' @param problem The sense of the linear constraints.
#' @inheritParams dkqs_cone
#'
#' @returns Returns the solution to the quadratic program.
#'
#' @export
prog_cone <- function(A_obs, A_tgt, beta_obs_hat, beta_tgt, tau, problem){
  #### Step 1: Formulation of the objective function for (5)
  rn = dim(A_tgt)[1]
  cn = dim(A_tgt)[2]
  obj2 = t(A_obs) %*% A_obs * rn
  obj1 = -2 * t(A_obs) %*% beta_obs_hat * rn
  obj0 = t(beta_obs_hat) %*% beta_obs_hat * rn

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
  rhs_down = (theta_up$objval - beta_tgt)*tau / length(c(ind_0, ind_down))
  rhs_0 = (1 - rhs_up * length(ind_up) - rhs_down * length(ind_down)) *
            tau / length(ind_0)
  # LP rhs and sense
  lp_rhs = c(beta_tgt, 1)
  lp_sense = c("=", "=")

  #### Step 3: Solve the QP
  # If problem == "T", it is used to solve linear program (5)
  if (problem == "T"){
  # Update lb
    lb_new = lb
    lb_new[ind_up] = rhs_up
    lb_new[ind_down] = rhs_down
    lb_new[ind_0] = rhs_0
    ### Use Gurobi to find the optimum
    ans = gurobi_optim(obj2, obj1, obj0, rbind(A_tgt, ones), lp_rhs, lp_sense,
                       "min", lb_new)
  }
  # If problem == "tau", it is used to solve linear program (6)
  else if (problem == "tau"){
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

#' Gurobi solver for quadratic or linear programs
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


#' Computes the test statistics via bootstrapping
#'
#' @description This function computes the test statistics via bootstrapping.
#'
#' @require modelr
#'
#' @param J The number of distinct nonzero values in vector \eqn{\bm{y}}.
#' @param s_star The value of s_star in the cone-tightening procedure.
#' @param beat_obs_hat The value of beta_obs_hat using the full data.
#' @param beta_bs_bar The tau-tightened recentered bootstrap estimate.
#' @inheritParams dkqs_cone
#' @inheritParams prog_cone
#'
#' @returns Returns the list of test statistics obtained from bootstrap, i.e.
#'    \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.
#'
#' @export
beta_bs <- function(df, bs_seed, bs_num, J, s_star, A_obs, A_tgt,
                    func_obs, beta_obs_hat, beta_bs_bar, beta_tgt, tau, N){
  T_bs = NULL
  # Loop through all indices in the bootstrap
  for (i in 1:bs_num){
    #### Step 1: Set the seed
    set.seed(bs_seed + i)
    ####  Step 2: Draw the subsample
    # Drop unobserved values, i.e. when D = 0
    df = df[df[,"D"] == 1,]
    df_bs = as.data.frame(resample_bootstrap(as.data.frame(df)))
    # Re-index the rows
    rownames(df_bs) = 1:nrow(df_bs)
    ####  Step 3: Compute the bootstrap estimates
    beta_bs_star = NULL
    # Compute the value of beta_bs_star using the function func_obs
    beta_bs_star = func_obs(df_bs)
    ####  Step 4: Compute the bootstrap test statistic
    beta_bs_bar = beta_bs_star - beta_obs + s_star
    T_bs_i = prog_cone(A_obs, A_tgt, beta_bs_bar, beta_tgt, tau, "T")$objval
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
  lp_sense_tau = c(lp_sense_tau, "<=")
  return(list(lp_lhs_tau = lp_lhs_tau,
              lp_rhs_tau = lp_rhs_tau,
              lp_sense_tau = lp_sense_tau))
}

#' Error messages generator
#'
#' @description This function generates the error messages and stops the program
#'    if there is any invalid input.
#'
#' @inheritParams dkqs_cone
#'
#' @export
dkqs_cone_errormsg <- function(df, A_obs, A_tgt, beta_tgt, bs_seed, bs_num, 
                               p_sig, tau_input){

  if (bs_num <= 0 | bs_num %%1 != 0){
    stop("The number of bootstrap (bs_num) has to be a positive integer.")
  }
  if (p_sig <= 0 | p_sig %%1 != 0){
    stop("The number of decimal places in the p-value (p_sig) has to be a
         positive integer.")
  }
  if (tau_input < 0 | tau_input > 1){
    stop("The value of tau has to be in [0,1].")
  }
}

