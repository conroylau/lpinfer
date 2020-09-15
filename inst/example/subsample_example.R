## ========================================================================= ##
##
##  Example for the subsample function
##
##  This followings illustrate how the function can be used to compute
##  p-values using the missing data problem.
##
## ========================================================================= ##
rm(list = ls())

# ---------------- #
# Part 1: Load required packages
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
reps <- 100
phi <- 2/3

# ---------------- #
# Part 3: Run the subsampling procedure
# ---------------- #
## Full information approach
# Example 1.1 - Using full information approach and gurobi solver (1 core)
set.seed(1)
plan(multisession, workers = 1)
ss1a <- subsample(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  solver = "gurobi",
                  norm = 2,
                  phi = phi,
                  replace = FALSE,
                  progress = FALSE)

# Example 1.2 - Using full information approach and gurobi solver (8 cores)
set.seed(1)
plan(multisession, workers = 8)
ss1b <- subsample(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  solver = "gurobi",
                  norm = 2,
                  phi = phi,
                  replace = FALSE,
                  progress = FALSE)

# Example 1.3 - Using two moments approach and gurobi solver
set.seed(1)
ss1c <- subsample(data = data,
                  lpmodel = lpmodel.twom,
                  beta.tgt = beta.tgt,
                  R = 100,
                  solver = "gurobi",
                  norm = 2,
                  phi = phi,
                  replace = FALSE,
                  progress = FALSE)

## Two moments approach
# Example 2.1 - Using the two moments approach with bootstrap
# (i.e. replace and phi = 1)
set.seed(1)
ss2a <- subsample(data = data,
                  lpmodel = lpmodel.twom,
                  beta.tgt = beta.tgt,
                  R = 100,
                  solver = "gurobi",
                  norm = 2,
                  phi = 1,
                  replace = TRUE,
                  progress = FALSE)

# Example 3.1 - Using the two moments approach with m out of n bootstrap
# (i.e. replace and phi != 1)
set.seed(1)
ss3a <- subsample(data = data,
                  lpmodel = lpmodel.twom,
                  beta.tgt = beta.tgt,
                  R = 100,
                  solver = "gurobi",
                  norm = 2,
                  phi = .5,
                  replace = TRUE,
                  progress = FALSE)
