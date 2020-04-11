## Test file for subsample
context("Tests for subsample")

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
### Full information approach
# Beta
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
# Omega_hat
var_full_info <- function(data){
  len = length(unique(data[,"Y"]))
  return(diag(len))
}

### Two moments approach
# Beta
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
# Omega_hat
var_two_moment <- function(data){
  return(diag(2))
}

##------------------------------------------------------------------------------
## Data preparation and declaration
##------------------------------------------------------------------------------
# Read data
data = sampledata
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
## Shape constraints
##------------------------------------------------------------------------------
A_shp_full = matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom = matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp = 1

##------------------------------------------------------------------------------
### Define arguments and produce output
##------------------------------------------------------------------------------
# Parameters to test
beta_tgt = .365
phi_predefine = 2/3

### Define arguments
farg = list(data = data,
            A_tgt = A_tgt,
            A_shp = A_shp_full,
            beta_shp = beta_shp,
            R = 100,
            beta_tgt = beta_tgt,
            cores = 8,
            norm = 2,
            alpha = .05,
            phi = 2/3,
            progress = TRUE)

### Generate output for full-information appraoch
# Append full-information arguments
farg$A_obs = A_obs_full
farg$func_obs = func_full_info
farg$func_var = var_full_info

# Gurobi
farg$solver = "gurobi"
full_g2 = do.call(subsample, farg)

# Rcplex
farg$solver = "rcplex"
full_r2 = do.call(subsample, farg)

# limSolve
farg$solver = "limsolve"
full_l2 = do.call(subsample, farg)

### L1 norm
farg$norm = 1
# Gurobi
farg$solver = "gurobi"
full_g1 = do.call(subsample, farg)

# Rcplex
farg$solver = "rcplex"
full_r1 = do.call(subsample, farg)

# limSolve
farg$solver = "limsolve"
full_l1 = do.call(subsample, farg)

### Generate output for two-moments appraoch
# Append full-information arguments
farg$A_obs = A_obs_twom
farg$func_obs = func_two_moment
farg$func_var = var_two_moment

# Gurobi
farg$solver = "gurobi"
twom_g2 = do.call(subsample, farg)

# Rcplex
farg$solver = "rcplex"
twom_r2 = do.call(subsample, farg)

# limSolve
farg$solver = "limsolve"
twom_l2 = do.call(subsample, farg)

### L1 norm
farg$norm = 1
# Gurobi
farg$solver = "gurobi"
twom_g1 = do.call(subsample, farg)

# Rcplex
farg$solver = "rcplex"
twom_r1 = do.call(subsample, farg)

# limSolve
farg$solver = "limsolve"
twom_l1 = do.call(subsample, farg)

##------------------------------------------------------------------------------
## Test 1: Test equivalence of two moments approach and full information 
## appraoch for each optimizer with L2 norm
##------------------------------------------------------------------------------

# Function to test the output with same solver and two different approaches
# subsample_return1 = function 1
# subsample_return1 = function 2
# dp = number of decimal places
subsample_test_output_samesolver <- function(subsample_return1, 
                                             subsample_return2, 
                                             dp){
  expect_equal(subsample_return1$p_val, subsample_return2$p_val)
  expect_equal(subsample_return1$solver, subsample_return2$solver)
  expect_equal(subsample_return1$cores, subsample_return2$cores)
  expect_equal(subsample_return1$norm, subsample_return2$norm)
  expect_equal(subsample_return1$decision, subsample_return2$decision)
}

### L1 norm
# Gurobi
test_that("Gurobi solver with L1-norm",{
  subsample_test_output_samesolver(full_g1, twom_g1)
})

# Rcplex solver
test_that("Rcplex solver with L1-norm",{
  subsample_test_output_samesolver(full_r1, twom_r1)
})

# limSolve solver
test_that("limSolve solver with L1-norm",{
  subsample_test_output_samesolver(full_l1, twom_l1)
})

### L2 norm
# Gurobi
test_that("Gurobi solver with L2-norm",{
  subsample_test_output_samesolver(full_g2, twom_g2)
})

# Rcplex solver
test_that("Rcplex solver with L2-norm",{
  subsample_test_output_samesolver(full_r2, twom_r2)
})

# limSolve solver
test_that("limSolve solver with L2-norm",{
  subsample_test_output_samesolver(full_l2, twom_l2)
})

##------------------------------------------------------------------------------
## Test 2: Test equivalence of results across different optimizers for each 
## approach - L2 norm
##------------------------------------------------------------------------------

# Function to test the output with same solver with different approach
# subsample_return1 = function 1
# subsample_return1 = function 2
# dp = number of decimal places
subsample_test_output_approach <- function(subsample_return1, 
                                           subsample_return2, 
                                           dp){
  expect_equal(subsample_return1$p_val, subsample_return2$p_val)
  expect_equal(round(subsample_return1$T_n, digits = dp), 
               round(subsample_return2$T_n, digits = dp))
  expect_equal(subsample_return1$cores, subsample_return2$cores)
  expect_equal(subsample_return1$norm, subsample_return2$norm)
  expect_equal(subsample_return1$decision, subsample_return2$decision)
}

# Set decimal places
dp = 5

### L1 norm
# Full information approach
test_that("Full information approach - Gurobi vs Rcplex with L1-norm",{
  subsample_test_output_approach(full_g1, full_r1, dp)
})

test_that("Full information approach - Rcplex vs limSolve with L1-norm",{
  subsample_test_output_approach(full_r1, full_l1, dp)
})

# Two moments approach
test_that("Two moments approach - Gurobi vs Rcplex with L1-norm",{
  subsample_test_output_approach(twom_g1, twom_r1, dp)
})

test_that("Two moments approach - Rcplex vs limSolve with L1-norm",{
  subsample_test_output_approach(twom_r1, twom_l1, dp)
})


### L2 norm
# Full information approach
test_that("Full information approach - Gurobi vs Rcplex with L2-norm",{
  subsample_test_output_approach(full_g2, full_r2, dp)
})

test_that("Full information approach - Rcplex vs limSolve with L2-norm",{
  subsample_test_output_approach(full_r2, full_l2, dp)
})

# Two moments approach
test_that("Two moments approach - Gurobi vs Rcplex with L2-norm",{
  subsample_test_output_approach(twom_g2, twom_r2, dp)
})

test_that("Two moments approach - Rcplex vs limSolve with L2-norm",{
  subsample_test_output_approach(twom_r2, twom_l2, dp)
})


