################################################################################
##
##  Example file for dkqs_cone function
##  
##  This is an example code for applying the R module dkqs_cone on the missing
##  data problem using the sample data by Torgovitsky (2019). This file 
##  illustrates how the module can be used to obtain the p-values using the 
##  full-information and two moments approach, and different linear and 
##  quadratic programming solvers in R. Currently, the solvers supported are
##  as follows:
##    - gurobi (LP & QP)
##    - cplexAPI (LP, QP in progress)
##    - Rcplex (LP & QP)
##    - lpSolveAPI (LP)
##    - osqp (QP)
##
################################################################################

### Part 1: Load required packages
library(modelr)
library(gurobi)
library(e1071)
library(lpSolveAPI)
library(cplexAPI)
library(quadprog)
library(osqp)
library(Rcplex)
library(ddpcr)

### Part 2: Data preparation
# Read data
df = read.csv("sample-data.csv")
# Compute parameters required
N = dim(df)[1]
J = length(unique(df[,"Y"])) - 1
J1 = J + 1
pi = 1 - mean(df[,"D"])
# Compute matrices required
yp = seq(0,1,1/J)
A_tgt = matrix(c(yp, yp), nrow = 1)
beta_obs = matrix(rep((1 - pi)/J1,J1), nrow = J1)
num_bs = 10
# Define the value of tau to be used
tau = sqrt(log(N)/N)

### Part 3: Define functions to compute beta_obs_hat
# Full information approach
func_full_info <- function(df){
  # Initialize beta
  beta = NULL
  # Find the unique elements of Y, sorting in ascending order
  y_list = sort(unique(df[,"Y"]))
  # Count total number of rows of df and y_list
  n = dim(df)[1]
  yn = length(y_list)
  # Generate each entry of beta_obs
  for (i in 1:yn){
    beta_i = sum((df[,"Y"] == y_list[i]) * (df[,"D"] == 1))/n
    beta = c(beta,c(beta_i))
  }
  beta = as.matrix(beta)
  return(beta)
}
# Two moments approach
func_two_moment <- function(df){
  # Initialize beta
  beta = matrix(c(0,0), nrow = 2)
  # Count total number of rows of df and y_list
  n = dim(df)[1]
  # Moment 1 E[YD]
  beta[1] = sum(df[,"Y"] * df[,"D"])/n
  # Moment 2 E[D]
  beta[2] = sum(df[,"D"])/n
  return(beta)
}

### Part 4: Compute A_obs based on the two functions in part 3
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                byrow = TRUE)

### Part 5: Run the dkqs_cone module to compute p-values
# Define the value of beta_tgt
beta_tgt = .365
# Example 1 - Using full information approach, gurobi for LP and QP
full_g_g = dkqs_cone(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 100, 4, 
                     tau, "gurobi", "gurobi")

# Example 2 - Using full information approach, gurobi for LP and osqp for QP
full_g_o = dkqs_cone(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 100, 4, 
                     tau, "gurobi", "osqp")

# Example 3 - Using two moments approach, Rcplex for LP and QP
twom_r_r = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 4, 
                     tau, "Rcplex", "Rcplex")

# Example 4 - Using two moments approach, lpSolveAPI for LP and Rcplex for QP
twom_l_r = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 4, 
                     tau, "lpsolveapi", "Rcplex")

# Example 5 - Using two moments approach, cplexAPI for LP and Rcplex for QP
twom_c_r = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 4, 
                     tau, "cplexapi", "Rcplex")

