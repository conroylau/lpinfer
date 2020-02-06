################################################################################
##
##  Example file to compare the time for different number of cores
##
################################################################################

### Part 1: Load required packages
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

### Part 2: Data preparation
# Read data
df = read.csv("./data/sampledata.csv")
# Compute parameters required
N = dim(df)[1]
J = length(unique(df[,"Y"])) - 1
J1 = J + 1
pi = 1 - mean(df[,"D"])
# Compute matrices required
yp = seq(0,1,1/J)
A_tgt = matrix(c(yp, yp), nrow = 1)
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

### Part 5: Run the dkqs module to compute p-values
# Define the value of beta_tgt and significant figures needed
beta_tgt = .365
p_sig = 4

# Example 1 - Using full information approach and gurobi solver (1 core)
t10 = Sys.time()
full_gur = dkqs(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 1000, 
                     p_sig, tau, "gurobi", 1, FALSE)
t11 = Sys.time()
time1 = t11 - t10

# Example 2 - Using full information approach and gurobi solver (n cores)
cores = 8
t80 = Sys.time()
full_gur = dkqs(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 1000, 
                     p_sig, tau, "gurobi", cores, FALSE)
t81 = Sys.time()
time8 = t81 - t80

# Print the time used
print(sprintf("Time used with 1 core: %s", time1))
print(sprintf("Time used with %s cores: %s", cores, time8))