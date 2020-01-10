context("Missing data problem")

# Load packages
library(modelr)
library(gurobi)
library(cplexAPI)
library(Rcplex)
library(Momocs)
library(limSolve)

##------------------------------------------------------------------------------
## Define functions to match the moments
##------------------------------------------------------------------------------
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

##------------------------------------------------------------------------------
## Data preparation and declaration
##------------------------------------------------------------------------------
# Read data
df = sampledata
# Declare parameters
N = dim(sampledata)[1]
J = length(unique(sampledata[,"Y"])) - 1
J1 = J + 1
pi = 1 - mean(sampledata[,"D"])
# Compute matrices required
yp = seq(0,1,1/J)
A_tgt = matrix(c(yp, yp), nrow = 1)
# Define the value of tau to be used
tau = sqrt(log(N)/N)
# Define the observed matrix for each appraoch
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)

##------------------------------------------------------------------------------
## Obtain results for each solver
##------------------------------------------------------------------------------
# Parameters to test
beta_tgt = .365
p_sig = 4

### Gurobi solver
# Full information
full_g = dkqs_cone(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 100, 
                   p_sig, tau, "gurobi")
# Two moments
twom_g = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 
                   p_sig, tau, "gurobi")
### Rcplex solver
# Full information
full_r = dkqs_cone(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 100, 
                   p_sig, tau, "rcplex")
# Two moments
twom_r = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 
                   p_sig, tau, "rcplex")

### limSolve solver
# Full information
full_l = dkqs_cone(df, A_obs_full, A_tgt, func_full_info, beta_tgt, 1, 100, 
                   p_sig, tau, "limsolve")
# Two moments
twom_l = dkqs_cone(df, A_obs_twom, A_tgt, func_two_moment, beta_tgt, 1, 100, 
                   p_sig, tau, "limsolve")

##------------------------------------------------------------------------------
## Test 1: Test equivalence of two moments approach and full information 
## appraoch for each optimizer
##
## * I only test whether the p-value and the value of tau used are the same 
##   because the value of bootstrap betas and test statistic could vary across
##   the two approaches because the two appraoches use different moment 
##   information
##------------------------------------------------------------------------------

# Gurobi
test_that("Gurobi solver",{
  expect_equal(full_g$tau, twom_g$tau)
  expect_equal(full_g$p_val, twom_g$p_val)
})

# Rcplex solver
test_that("Rcplex solver",{
  expect_equal(full_r$tau, twom_r$tau)
  expect_equal(full_r$p_val, twom_r$p_val)
})

# limSolve solver
test_that("limSolve solver",{
  expect_equal(full_l$tau, twom_l$tau)
  expect_equal(full_l$p_val, twom_l$p_val)
})

##------------------------------------------------------------------------------
## Test 2: Test equivalence of results across different optimizers for each 
## approach
##
## * In this section, I test only the output of tau, p-value and test statistic
##   from each solver for each approach
##------------------------------------------------------------------------------

# Full information approach
test_that("Full information approach",{
  # Tau
  expect_equal(full_g$tau, full_r$tau)
  expect_equal(full_r$tau, full_l$tau)
  # p-value
  expect_equal(full_g$p_val, full_r$p_val)
  expect_equal(full_r$p_val, full_l$p_val)
  # Test statistic
  expect_equal(full_g$T_n, full_r$T_n)
  expect_equal(full_r$T_n, full_l$T_n)
})

# Two moments approach
test_that("Two moments approach",{
  # Tau
  expect_equal(twom_g$tau, twom_r$tau)
  expect_equal(twom_r$tau, twom_l$tau)
  # p-value
  expect_equal(twom_g$p_val, twom_r$p_val)
  expect_equal(twom_r$p_val, twom_l$p_val)
  # Test statistic
  expect_equal(twom_g$T_n, twom_r$T_n)
  expect_equal(twom_r$T_n, twom_l$T_n)
})


