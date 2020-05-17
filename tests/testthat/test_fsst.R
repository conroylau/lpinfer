## Test file for dkqs
context("Tests for fsst")

# ---------------- #
# Load packages
# ---------------- #
library(modelr)
library(gurobi)
library(cplexAPI)
library(Rcplex)
library(Momocs)
library(limSolve)
library(foreach)
library(doMC)
library(lpinfer)

# ---------------- #
# Define functions to match the moments
# ---------------- #
### Full information approach
# Function for Omega_hat
var_full_info <- function(data){
  len <- length(unique(data[,"Y"]))
  return(diag(len))
}

# Function for beta_hat
func_full_info <- function(data){
  # Initialize beta
  beta <- NULL
  # Find the unique elements of Y, sorting in ascending order
  y_list <- sort(unique(data[,"Y"]))
  # Count total number of rows of data and y_list
  n <- dim(data)[1]
  yn <- length(y_list)
  # Generate each entry of beta_obs
  for (i in 1:yn) {
    beta_i <- sum((data[,"Y"] == y_list[i]) * (data[,"D"] == 1))/n
    beta <- c(beta,c(beta_i))
  }
  beta <- as.matrix(beta)
  # Variance
  var <- var_full_info(data)
  return(list(beta = beta,
              var = var))
}

### Two moments approach
# Function for Omega_hat
var_two_moment <- function(data){
  return(diag(2))
}

# Function for beta_hat
func_two_moment <- function(data){
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)
  # Count total number of rows of data and y_list
  n <- dim(data)[1]
  # Moment 1 E[YD]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  # Moment 2 E[D]
  beta[2] <- sum(data[,"D"])/n
  # Variance
  var <- var_two_moment(data)
  return(list(beta = beta,
              var = var))
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
# Shape constraints
# ---------------- #
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom <- matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp <- c(1)

# ---------------- #
# Define arguments and produce output
# ---------------- #
# Parameters to test
beta.tgt <- .365
phi_predefine <- 2/3

# Define the lpmodels
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_full,
                        beta.obs = func_full_info,
                        beta.shp = beta_shp)

lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_full,
                        beta.obs = func_two_moment,
                        beta.shp = beta_shp)

### Define arguments
farg <- list(data = data,
             beta.tgt = beta.tgt,
             cores = 8,
             alpha = .05,
             lambda = .5,
             rho = 1e-4,
             n = NULL,
             R = 10,
             weight.matrix = "diag",
             solver = "gurobi",
             progress = TRUE)

# ---------------- #
# Test whether the function works when `beta.obs' is a function
# ---------------- #
### Case 1: A single lambda
# Append full-information arguments
farg$lpmodel <- lpmodel.full
set.seed(1)
full_g2 <- do.call(fsst, farg)
beta.obs.bs.full <- full_g2$beta.obs.bs

# Append two-moments arguments
farg$lpmodel <- lpmodel.twom
set.seed(1)
twom_g2 <- do.call(fsst, farg)
beta.obs.bs.twom <- twom_g2$beta.obs.bs

# Check if the p-values are equal in the two approaches
test_that("p-values for full information and two moments approach",{
  expect_equal(full_g2$pval, twom_g2$pval)
})

### Case 2: Single lambda with different cores
# Cores = 1
farg$cores <- 1
farg$lpmodel <- lpmodel.full
set.seed(1)
full_core1 <- do.call(fsst, farg)

# Cores = 8
farg$cores <- 8
set.seed(1)
full_core8 <- do.call(fsst, farg)

# Check if the p-values are equal in parallel and non-parallel approaches
test_that("beta.obs is a function (two m vs full)",{
  expect_equal(full_core1$pval, full_core8$pval)
})

### Case 3: Multiple lambda
farg$lambda <- c(.5, .9)

# Append full-information arguments
farg$lpmodel <- lpmodel.full
set.seed(1)
full_g2_multiple <- do.call(fsst, farg)

# Append two-moments arguments
farg$lpmodel <- lpmodel.twom
set.seed(1)
twom_g2_multiple <- do.call(fsst, farg)

# Check if the p-values are equal in the two approaches
test_that("Multiple lambdas",{
  expect_equal(full_g2_multiple$pval, twom_g2_multiple$pval)
})

# ---------------- #
# Test whether the function works when `beta.obs' is a list
# ---------------- #
# Consolidate the `beta.obs` into a list - full information
lpmodel.full.list <- lpmodel.full
beta.obs.full.estimator <- lpmodel.beta.eval(data,
                                             lpmodel.full$beta.obs, 1)[[1]]
lpmodel.full.list$beta.obs <- c(list(beta.obs.full.estimator),
                                beta.obs.bs.full)

# Consolidate the `beta.obs` into a list - two moments
lpmodel.twom.list <- lpmodel.twom
beta.obs.twom.estimator <- lpmodel.beta.eval(data,
                                             lpmodel.twom$beta.obs, 1)[[1]]
lpmodel.twom.list$beta.obs <- c(list(beta.obs.twom.estimator),
                                beta.obs.bs.twom)

# Two moments approach
set.seed(1)
farg$lpmodel <- lpmodel.full.list
farg$lambda <- .5
full_g2_list <- do.call(fsst, farg)

set.seed(1)
farg$lpmodel <- lpmodel.twom.list
twom_g2_list <- do.call(fsst, farg)

test_that("beta.obs is a list (full vs two m)",{
  expect_equal(full_g2_list$pval, twom_g2_list$pval)
})

test_that("Full information appraoch: function vs list",{
  expect_equal(full_g2$pval, full_g2_list$pval)
})

test_that("Two moments appraoch: function vs list",{
  expect_equal(twom_g2$pval, twom_g2_list$pval)
})
