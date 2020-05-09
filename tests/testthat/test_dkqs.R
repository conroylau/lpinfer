## Test file for dkqs
context("Tests for dkqs")

# Load packages
library(modelr)
library(gurobi)
library(cplexAPI)
library(Rcplex)
library(limSolve)
library(foreach)
library(doMC)

# ---------------- #
# Define functions to match the moments
# ---------------- #
# Full information approach
func_full_info <- function(data){
  # Initialize beta
  beta <- NULL
  # Find the unique elements of Y, sorting in ascending order
  y_list <- sort(unique(data[,"Y"]))
  # Count total number of rows of data and y_list
  n <- dim(data)[1]
  yn <- length(y_list)
  # Generate each entry of beta_obs
  for (i in 1:yn){
    beta_i <- sum((data[,"Y"] == y_list[i]) * (data[,"D"] == 1))/n
    beta <- c(beta,c(beta_i))
  }
  beta <- as.matrix(beta)
  return(beta)
}

# Two moments approach
func_two_moment <- function(data){
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)
  # Count total number of rows of data and y_list
  n <- dim(data)[1]
  # Moment 1 E[YD]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  # Moment 2 E[D]
  beta[2] <- sum(data[,"D"])/n
  return(beta)
}

# ---------------- #
# Data preparation and declaration
# ---------------- #
# Read data
data <- sampledata
# Declare parameters
N <- dim(sampledata)[1]
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])
# Compute matrices required
yp <- seq(0,1,1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)
# Define the value of tau to be used
tau <- sqrt(log(N)/N)
# Define the observed matrix for each appraoch
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)

# ---------------- #
# Define arguments and produce output
# ---------------- #
# Parameters to test
beta.tgt <- .365

# Define the value of full information method
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        beta.obs = func_full_info)

# Define the value of full two moments method
lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        beta.obs = func_two_moment)

# Define arguments
farg <- list(data = data,
             beta.tgt = beta.tgt,
             R = 100,
             tau = tau,
             cores = 8,
             progress = TRUE)

### Generate output for full-information appraoch
# Append full-information arguments
farg$lpmodel <- lpmodel.full

# Gurobi
farg$solver = "gurobi"
set.seed(1)
full_g = do.call(dkqs, farg)

# Rcplex
farg$solver = "rcplex"
set.seed(1)
full_r = do.call(dkqs, farg)

# limSolve
farg$solver = "limsolve"
set.seed(1)
full_l = do.call(dkqs, farg)

### Generate output for two-moments appraoch
# Append full-information arguments
farg$lpmodel <- lpmodel.twom

# Gurobi
farg$solver = "gurobi"
set.seed(1)
twom_g = do.call(dkqs, farg)

# Rcplex
farg$solver = "rcplex"
set.seed(1)
twom_r = do.call(dkqs, farg)

# limSolve
farg$solver = "limsolve"
set.seed(1)
twom_l = do.call(dkqs, farg)

# ---------------- #
# Test 1: Test equivalence of two moments approach and full information
# appraoch for each optimizer
# ---------------- #

# Function to test the output with same solver and two different approaches
# dkqs_return1 = function 1
# dkqs_return1 = function 2
# dp = number of decimal places
dkqs_test_output_samesolver <- function(dkqs_return1, dkqs_return2, dp){
  expect_equal(dkqs_return1$pval, dkqs_return2$pval)
  expect_equal(dkqs_return1$tau, dkqs_return2$tau)
  expect_equal(dkqs_return1$lb0, dkqs_return2$lb0)
  expect_equal(dkqs_return1$ub0, dkqs_return2$ub0)
  expect_equal(dkqs_return1$solver, dkqs_return2$solver)
  expect_equal(dkqs_return1$cores, dkqs_return2$cores)
}

# Gurobi
test_that("Gurobi solver",{
  dkqs_test_output_samesolver(full_g, twom_g)
})

# Rcplex solver
test_that("Rcplex solver",{
  dkqs_test_output_samesolver(full_r, twom_r)
})

# limSolve solver
test_that("limSolve solver",{
  dkqs_test_output_samesolver(full_l, twom_l)
})

# ---------------- #
# Test 2: Test equivalence of results across different optimizers for each
# approach
# ---------------- #

# Function to test the output with same solver with different approach
# dkqs_return1 = function 1
# dkqs_return1 = function 2
# dp = number of decimal places
dkqs_test_output_approach <- function(dkqs_return1, dkqs_return2, dp){
  expect_equal(dkqs_return1$pval, dkqs_return2$pval)
  expect_equal(dkqs_return1$tau, dkqs_return2$tau)
  expect_equal(round(dkqs_return1$T.n, digits = dp),
               round(dkqs_return2$T.n, digits = dp))
  expect_equal(dkqs_return1$lb0, dkqs_return2$lb0)
  expect_equal(dkqs_return1$ub0, dkqs_return2$ub0)
  expect_equal(dkqs_return1$cores, dkqs_return2$cores)
}

# Set decimal places
dp <- 5

# Full information approach
test_that("Full information approach - Gurobi vs Rcplex",{
  dkqs_test_output_approach(full_g, full_r, dp)
})

test_that("Full information approach - Rcplex vs limSolve",{
  dkqs_test_output_approach(full_r, full_l, dp)
})

# Two moments approach
test_that("Two moments approach - Gurobi vs Rcplex",{
  dkqs_test_output_approach(twom_g, twom_r, dp)
})

test_that("Two moments approach - Rcplex vs limSolve",{
  dkqs_test_output_approach(twom_r, twom_l, dp)
})
