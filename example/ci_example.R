################################################################################
##
##  Example file for qpci function
##  
##  This is an example code for applying the R module qpci on the missing
##  data problem using the sample data by Torgovitsky (2019). This file 
##  illustrates how the module can be used to obtain the confidence interval
##  for the cone-tightening procedure problem.
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
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)

### Part 5: Run the dkqs_cone module to compute p-values
# Define the value of beta_tgt and significant figures needed
beta_tgt = .365
p_sig = 4

### Part 6: Arguments for the dkqs_cone function without beta_tgt
farg = list(df = df,
            A_obs = A_obs_full,
            A_tgt = A_tgt,
            func_obs = func_full_info,
            bs_seed = 1,
            bs_num = 100,
            p_sig = p_sig,
            tau_input = tau,
            solver = "gurobi",
            noisy = FALSE)

### Demonstration 1: Construction of confidence interval
qpci(dkqs_cone, farg, 0.05, 0, 0.4, 1, .6, 0.001, 10, NULL, TRUE)

### Demonstration 2: Use wrapper function to construct a list of 
### confidence intervals
many_qpci(dkqs_cone, farg, c(0.01, 0.02, 0.05, 0.1), 0, 0.4, .6, 1, 
          0.0000001, 10, NULL, FALSE, TRUE)
