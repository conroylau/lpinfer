## Test file for invert CI
context("Constructing confidence interval")

# Load packages
library(modelr)
library(gurobi)
library(cplexAPI)
library(Rcplex)
library(Momocs)
library(limSolve)
library(foreach)
library(doMC)

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
## Construct arguments for dkqs command for each solver and approach
##------------------------------------------------------------------------------
# Parameters to test
beta_tgt = .365
p_sig = 4

### Define arguments for dkqs
farg = list(df = df,
            A_tgt = A_tgt,
            bs_seed = 1,
            bs_num = 100,
            p_sig = p_sig,
            tau_input = tau,
            cores = 8,
            progress = FALSE)

### Arguments for full-information approach
# Matrices
farg$A_obs = A_obs_full
farg$func_obs = func_full_info

# Gurobi
farg$solver = "gurobi"
farg_g_full = farg

# Rcplex
farg$solver = "rcplex"
farg_r_full = farg

# limSolve
farg$solver = "limsolve"
farg_l_full = farg


### Generate output for two-moments appraoch
# Matrices
farg$A_obs = A_obs_twom
farg$func_obs = func_two_moment

# Gurobi
farg$solver = "gurobi"
farg_g_twom = farg

# Rcplex
farg$solver = "rcplex"
farg_r_twom = farg

# limSolve
farg$solver = "limsolve"
farg_l_twom = farg

##------------------------------------------------------------------------------
## Construct arguments for constructing confidence interval
##------------------------------------------------------------------------------
invertci_arg = list(f = dkqs,
                alpha = 0.05, 
                lb0 = 0, 
                lb1 = 0.4, 
                ub0 = 1, 
                ub1 = .6, 
                tol = 0.001, 
                max_iter = 50, 
                df_ci = NULL, 
                progress = TRUE)

##------------------------------------------------------------------------------
## Construct the output
##------------------------------------------------------------------------------
### Gurobi & Full info
invertci_arg$farg = farg_g_full
invertci_g_full = do.call(invertci, invertci_arg)

### Gurobi & Two moments
invertci_arg$farg = farg_g_twom
invertci_g_twom = do.call(invertci, invertci_arg)

### Rcplex & Full info
invertci_arg$farg = farg_r_full
invertci_r_full = do.call(invertci, invertci_arg)

### Rcplex & Two moments
invertci_arg$farg = farg_r_twom
invertci_r_twom = do.call(invertci, invertci_arg)

### limSolve & Full info
invertci_arg$farg = farg_l_full
invertci_l_full = do.call(invertci, invertci_arg)

### limSolve & Two moments
invertci_arg$farg = farg_l_twom
invertci_l_twom = do.call(invertci, invertci_arg)

##------------------------------------------------------------------------------
## Test 1: Test equivalence of two moments approach and full information 
## appraoch for each optimizer
##------------------------------------------------------------------------------
dp = 5
# Function to test the output with same solver and two different approaches
# dkqs_return1 = function 1
# dkqs_return1 = function 2
# dp = number of decimal places
ci_test <- function(dkqs_return1, dkqs_return2, dp){
  expect_equal(dkqs_return1$up, dkqs_return2$up)
  expect_equal(dkqs_return1$down, dkqs_return2$down)
  expect_equal(dkqs_return1$df_ci, dkqs_return2$df_ci)
  expect_equal(dkqs_return1$df_down, dkqs_return2$df_down)
  expect_equal(dkqs_return1$df_up, dkqs_return2$df_up)
  expect_equal(dkqs_return1$tol, dkqs_return2$tol)
  expect_equal(dkqs_return1$iter, dkqs_return2$iter)
}

# Gurobi
test_that("Gurobi solver",{
  ci_test(invertci_g_full, invertci_g_twom, dp)
})

# Rcplex solver
test_that("Rcplex solver",{
  ci_test(invertci_r_full, invertci_r_twom, dp)
})

# limSolve solver
test_that("limSolve solver",{
  ci_test(invertci_l_full, invertci_l_twom, dp)
})

##------------------------------------------------------------------------------
## Test 2: Test equivalence of results across different optimizers for each 
## approach
##------------------------------------------------------------------------------

# Full information approach
test_that("Full information approach - Gurobi vs Rcplex",{
  ci_test(invertci_g_full, invertci_r_full, dp)
})

test_that("Full information approach - Rcplex vs limSolve",{
  ci_test(invertci_r_full, invertci_l_full, dp)
})

# Two moments approach
test_that("Two moments approach - Gurobi vs Rcplex",{
  ci_test(invertci_g_twom, invertci_r_twom, dp)
})

test_that("Two moments approach - Rcplex vs limSolve",{
  ci_test(invertci_r_twom, invertci_l_twom, dp)
})


