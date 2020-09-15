## ========================================================================= ##
##
##  Example for the estbounds function
##
##  This followings illustrate how the function can be used to estimate the
##  bounds subject to shape constraints using the missing data problem. 
##  Currently, only Gurobi is supported for the estbounds function with a
##  2-norm.
##
## ========================================================================= ##
rm(list = ls())

# ---------------- #
# Part 1: Load packages
# ---------------- #
library(lpinfer)
library(future)

# ---------------- #
# Part 2: Data and lpmodel preparation
# ---------------- #
source("./inst/example/dgp_missingdata.R")
J <- 5
N <- 1000
data <- missingdata_draw(J = J, n = N, seed = 1, prob.obs = .5)
lpmodel.full <- missingdata_lpm(J = J, info = "full", data = data)
lpmodel.twom <- missingdata_lpm(J = J, info = "mean", data = data)

tau <- sqrt(log(N)/N)
beta.tgt <- .2
kappa <- 1e-5

# ---------------- #
# Step 3: Run estbounds
# ---------------- #
# Example 1 - Compute the true bounds by setting estimate = FALSE
estb1 <- estbounds(data = data,
                   lpmodel = lpmodel.full,
                   kappa = kappa,
                   norm = 1,
                   solver = "gurobi",
                   estimate = FALSE)

# Example 2 - Estimated bounds with full-information approach and 1-norm
estb2 <- estbounds(data = data,
                   lpmodel = lpmodel.full,
                   kappa = kappa,
                   norm = 1,
                   solver = "gurobi",
                   estimate = FALSE)

# Example 3 - Estimated bounds with two-moments approach and 2-norm
estb3 <- estbounds(data = data,
                   lpmodel = lpmodel.twom,
                   kappa = kappa,
                   norm = 2,
                   solver = "gurobi",
                   estimate = FALSE)
