################################################################################
##
##  Example file for subsample function
##  
##  This is an example code for applying the R module subsample on the missing
##  data problem using the sample data by Torgovitsky (2019). This file 
##  illustrates how the module can be used to obtain the p-values using the 
##  full-information and two moments approach.
##
################################################################################

# = = = = = = 
# Part 1: Load required packages
# = = = = = = 
library(modelr)
library(gurobi)
library(e1071)
library(cplexAPI)
library(Rcplex)
library(ddpcr)
library(Momocs)
library(limSolve)
library(foreach)
library(doMC)
library(parallel)
library(PtProcess)
library(doSNOW)

# = = = = = = 
# Part 2: Data preparation
# = = = = = = 
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

# = = = = = = 
# Part 3: Define functions to compute beta_obs_hat and Omega_hat
# = = = = = = 
### Full information approach
# Function for beta_hat
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
# Function for Omega_hat
var_full_info <- function(data){
  len = length(unique(data[,"Y"]))
  return(diag(len))
}

### Two moments approach
# Function for beta_hat
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
# Function for Omega_hat
var_two_moment <- function(data){
  return(diag(2))
}


# = = = = = = 
# Part 4: Compute A_obs based on the two functions in part 3
# = = = = = = 
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)

# = = = = = = 
# Part 5: Introduce shape constraints
# = = = = = = 
A_shp_full = matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom = matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp = 1

# = = = = = = 
# Part 6: Run the subsample module to compute p-values
# = = = = = = 
# Define the parameters to be used
beta_tgt = .375
p_sig = 4
phi_predefine = 2/3
## Full information approach
# Example 1 - Using full information approach and gurobi solver (1 core)
full_gur = subsample(data = data, 
                     A_obs = A_obs_full, 
                     func_obs = func_full_info,
                     func_var = var_full_info,
                     A_shp = A_shp_full,
                     beta_shp = beta_shp,
                     A_tgt = A_tgt,
                     beta_tgt = beta_tgt,
                     R = 100,
                     bs_seed = 1,
                     p_sig = 3,
                     solver = "gurobi",
                     cores = 1,
                     lnorm = 2,
                     phi = phi_predefine,
                     progress = FALSE)

# Example 2 - Using the full information approach and gurobi solver and print
# the result after the operation is completed (8 cores)
subsample(data = data, 
          A_obs = A_obs_full, 
          func_obs = func_full_info,
          func_var = var_full_info,
          A_shp = A_shp_full,
          beta_shp = beta_shp,
          A_tgt = A_tgt,
          beta_tgt = beta_tgt,
          R = 100,
          bs_seed = 1,
          p_sig = 3,
          solver = "gurobi",
          cores = 8,
          lnorm = 2,
          phi = phi_predefine,
          progress = FALSE)


# Example 3 - Using two moments approach and gurobi solver (1 core)
twom_gur = subsample(data = data, 
                     A_obs = A_obs_twom, 
                     func_obs = func_two_moment,
                     func_var = var_two_moment,
                     A_shp = A_shp_full,
                     beta_shp = beta_shp,
                     A_tgt = A_tgt,
                     beta_tgt = beta_tgt,
                     R = 100,
                     bs_seed = 1,
                     p_sig = 3,
                     solver = "gurobi",
                     cores = 1,
                     lnorm = 2,
                     phi = phi_predefine,
                     progress = FALSE)