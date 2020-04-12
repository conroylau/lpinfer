################################################################################
##
##  Example file for the estbounds function
##  
##  This is an example code for applying the R module estbounds on the missing
##  data problem using the sample data by Torgovitsky (2019). This file 
##  illustrates how the module can be used to obtain the upper and lower bounds 
##  subject to the shape constraints. Currently, only the Gurobi optimzer
##  is supported to compute the bounds.
##
################################################################################

### Part 1: Load required packages
library(modelr)
library(gurobi)
library(e1071)
library(cplexAPI)
library(Rcplex)
library(lpSolveAPI)
library(ddpcr)
library(Momocs)
library(limSolve)
library(foreach)
library(doMC)
library(parallel)

### Part 2: Data preparation
# Read data
data = read.csv("./data/sampledata.csv")
# Compute parameters required
N = dim(data)[1]
J = length(unique(data[,"Y"])) - 1
J1 = J + 1
pi = 1 - mean(data[,"D"])
# Compute matrices required
yp = seq(0,1,1/J)
A_tgt = matrix(c(yp, yp), nrow = 1)
# Define the value of tau to be used
tau = sqrt(log(N)/N)

### Part 3: Define functions to compute beta_obs_hat
# Full information approach
func_full_info <- function(data){
  # Initialize beta
  beta = NULL
  # Find the unique elements of Y, sorting in ascending order
  y_list = sort(unique(data[,"Y"]))
  # Count total number of rows of data and y_list
  n = dim(data)[1]
  yn = length(y_list)
  # Generate each entry of beta_obs
  for (i in 1:yn){
    beta_i = sum((data[,"Y"] == y_list[i]) * (data[,"D"] == 1))/n
    beta = c(beta,c(beta_i))
  }
  beta = as.matrix(beta)
  return(beta)
}
# Two moments approach
func_two_moment <- function(data){
  # Initialize beta
  beta = matrix(c(0,0), nrow = 2)
  # Count total number of rows of data and y_list
  n = dim(data)[1]
  # Moment 1 E[YD]
  beta[1] = sum(data[,"Y"] * data[,"D"])/n
  # Moment 2 E[D]
  beta[2] = sum(data[,"D"])/n
  return(beta)
}

### Part 4: Compute A_obs based on the two functions in part 3
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)

### Part 5 Set the parameters and shape constraints
# Equality constraint
A_shp_dkqs = matrix(rep(1, ncol(A_obs_full)), nrow = 1)
beta_shp_dkqs = c(1)

### Part 5(a) Demonstration of first stage of the procedure
# L2-norm and full info approach
min1a = mincriterion(data = data,
                     func_obs = func_full_info,
                     A_obs = A_obs_full,
                     A_tgt = A_tgt,
                     A_shp = A_shp_dkqs,
                     beta_shp = beta_shp_dkqs,
                     norm = 2,
                     solver = "gurobi")

# L1-norm and full info approach
min1b = mincriterion(data = data,
                     func_obs = func_full_info,
                     A_obs = A_obs_full,
                     A_tgt = A_tgt,
                     A_shp = A_shp_dkqs,
                     beta_shp = beta_shp_dkqs,
                     norm = 1,
                     solver = "gurobi")

# L2-norm and two moments approach
min2a = mincriterion(data = data,
                     func_obs = func_two_moment,
                     A_obs = A_obs_twom,
                     A_tgt = A_tgt,
                     A_shp = A_shp_dkqs,
                     beta_shp = beta_shp_dkqs,
                     norm = 2,
                     solver = "gurobi")

# L1-norm and two moments approach
min2b = mincriterion(data = data,
                     func_obs = func_two_moment,
                     A_obs = A_obs_twom,
                     A_tgt = A_tgt,
                     A_shp = A_shp_dkqs,
                     beta_shp = beta_shp_dkqs,
                     norm = 1,
                     solver = "gurobi")

### Part 5(b) Demonstration of the full procedure 
# Set function arguments
farg = list(A_obs = A_obs_full,
            A_tgt = A_tgt,
            func_obs = func_full_info,
            A_shp = A_shp_dkqs,
            beta_shp = beta_shp_dkqs,
            kappa = 1e-20,
            norm = 2,
            solver = "gurobi",
            estimate = FALSE,
            progress = TRUE)

### True bounds
est_ans1 = estbounds(data = data,
                     A_obs = A_obs_full,
                     A_tgt = A_tgt,
                     func_obs = func_full_info,
                     A_shp = A_shp_dkqs,
                     beta_shp = beta_shp_dkqs,
                     kappa = 1e-20,
                     norm = 1,
                     solver = "gurobi",
                     estimate = FALSE,
                     progress = TRUE)


## Full-information approach
# L2-norm
est_ans2a = estbounds(data = data,
                      A_obs = A_obs_full,
                      A_tgt = A_tgt,
                      func_obs = func_full_info,
                      A_shp = A_shp_dkqs,
                      beta_shp = beta_shp_dkqs,
                      kappa = 1e-20,
                      norm = 2,
                      solver = "gurobi",
                      estimate = TRUE,
                      progress = TRUE)
# L1-norm
est_ans2b = estbounds(data = data,
                      A_obs = A_obs_full,
                      A_tgt = A_tgt,
                      func_obs = func_full_info,
                      A_shp = A_shp_dkqs,
                      beta_shp = beta_shp_dkqs,
                      kappa = 1e-20,
                      norm = 1,
                      solver = "gurobi",
                      estimate = TRUE,
                      progress = TRUE)

## Two-moments approach
# L2-norm
est_ans3a = estbounds(data = data,
                      A_obs = A_obs_twom,
                      A_tgt = A_tgt,
                      func_obs = func_two_moment,
                      A_shp = A_shp_dkqs,
                      beta_shp = beta_shp_dkqs,
                      kappa = 1e-20,
                      norm = 2,
                      solver = "gurobi",
                      estimate = TRUE,
                      progress = TRUE)
# L1-norm
est_ans3b = estbounds(data = data,
                      A_obs = A_obs_twom,
                      A_tgt = A_tgt,
                      func_obs = func_two_moment,
                      A_shp = A_shp_dkqs,
                      beta_shp = beta_shp_dkqs,
                      kappa = 1e-20,
                      norm = 1,
                      solver = "gurobi",
                      estimate = TRUE,
                      progress = TRUE)


