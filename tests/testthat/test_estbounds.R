## Test file for the missing data file
context("Estimation of bounds subject to shape constraints")

# ---------------- #
# Load required packages
# ---------------- #
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
  yn <-length(y_list)
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

# Define the observed matrix for each appraoch
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)

# ---------------- #
# Define shape constraints
# ---------------- #
beta_shp_dkqs <- c(1)
A_shp_dkqs <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)

# ---------------- #
# Define arguments and produce output for L2-norm
# ---------------- #
L2norm <- 2

# Define the value of full information method
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_dkqs,
                        beta.obs = func_full_info,
                        beta.shp = beta_shp_dkqs)

# Define the value of full two moments method
lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_dkqs,
                        beta.obs = func_two_moment,
                        beta.shp = beta_shp_dkqs)

# Define arguments
farg <- list(data = data,
             kappa = 1e-20,
             norm = L2norm,
             solver = "gurobi",
             estimate = FALSE,
             progress = TRUE)

# ---------------- #
# True bounds
# ---------------- #
## Generate true bounds for full-information appraoch
farg$lpmodel <- lpmodel.full
bd_full_true <- do.call(estbounds, farg)

## Generate true bounds for two-moments appraoch
farg$lpmodel <- lpmodel.twom
bd_twom_true <- do.call(estbounds, farg)

### Estimated bounds with kappa = 1e-10
farg$kappa <- 1e-10
farg$estimate <- TRUE
## Generate true bounds for full-information appraoch using kappa = 1e-10
farg$lpmodel <- lpmodel.full
bd_full_est10 <- do.call(estbounds, farg)

## Generate true bounds for two-moments appraoch using kappa = 1e-10
farg$lpmodel <- lpmodel.twom
bd_twom_est10 <- do.call(estbounds, farg)

## Estimated bounds with kappa = 1e-30
farg$kappa <- 1e-30
# Generate true bounds for full-information appraoch using kappa = 1e-30
farg$lpmodel <- lpmodel.full
bd_full_est30 <- do.call(estbounds, farg)

# Generate true bounds for two-moments appraoch using kappa = 1e-30
farg$lpmodel <- lpmodel.twom
bd_twom_est30 <- do.call(estbounds, farg)

# ---------------- #
# Test 1: Test whether the output from computing the true bounds in two
# different approaches yield the same results
# ---------------- #
# Function to test the output with same solver with different approach
# return1 = returned from function 1
# return2 = returned from function 2
# dp = number of decimal places
estbounds_test_1 <- function(return1, return2, dp){
  expect_equal(return1$ub, return2$ub)
  expect_equal(return1$lb, return2$lb)
  expect_equal(return1$est, return2$est)
}

test_that("Full info vs Two moments for true bounds",{
  estbounds_test_1(bd_full_true, bd_twom_true, 10)
})

# ---------------- #
# Test 2: Test whether the output from estimated bounds are same as the
# true bounds corrected 3 decimal places
# ---------------- #
# Function to test the output with same solver with different approach
# return1 = returned from function 1
# return2 = returned from function 2
# dp = number of decimal places
estbounds_test_2 <- function(return1, return2, dp){
  expect_equal(round(return1$ub, digits = dp), round(return2$ub, digits = dp))
  expect_equal(round(return1$lb, digits = dp), round(return2$lb, digits = dp))
}

test_that("Full info vs Two moments for estimated bounds and kappa = 1e-10",{
  estbounds_test_2(bd_full_est10, bd_twom_est10, 3)
})

test_that("Full info vs Two moments for estimated bounds and kappa = 1e-30",{
  estbounds_test_2(bd_full_est30, bd_twom_est30, 3)
})

test_that("Estimated bounds vs true bounds for kappa = 1e-10",{
  estbounds_test_2(bd_full_true, bd_full_est10, 3)
})

test_that("Estimated bounds vs true bounds for kappa = 1e-30",{
  estbounds_test_2(bd_full_true, bd_full_est30, 3)
})
