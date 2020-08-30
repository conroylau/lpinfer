## ========================================================================= ##
##
##  Example file for dkqs function
##
##  This is an example code for applying dkqs procedure on the missing
##  data problem using the sample data. This file illustrates how the module
##  can be used to obtain the p-values using the ull-information and two
##  moments approach.
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
source("./example/dgp_missingdata.R")
J <- 5
N <- 1000
data <- missingdata_draw(J = J, n = N, seed = 1, prob.obs = .5)
lpmodel.full <- missingdata_lpm(J = J, info = "full", data = data)
lpmodel.twom <- missingdata_lpm(J = J, info = "mean", data = data)

tau <- sqrt(log(N)/N)
beta.tgt <- .2
reps <- 100

# ---------------- #
# Part 3: Run the `dkqs` procedure
# ---------------- #
## Full information approach
# Example 1.1a - Using full information approach and gurobi solver (1 core)
set.seed(1)
plan(multisession, workers = 1)
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "gurobi",
                 progress = TRUE)

# Example 1.1b - Using full information approach and gurobi solver (8 cores)
set.seed(1)
plan(multisession, workers = 8)
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "gurobi",
                 progress = TRUE)

# Example 1.2 - Using full information approach and Rcplex solver
set.seed(1)
full_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "Rcplex",
                 progress = TRUE)

# Example 1.3 - Using full information approach and limSolve solver
set.seed(1)
full_lim <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "limSolve",
                 progress = TRUE)

## Two moments approach
# Example 2.1 - Using two moments approach and gurobi solver
set.seed(1)
twom_gur <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "gurobi",
                 progress = TRUE)

# Example 2.2 - Using two moments approach and Rcplex solver
set.seed(1)
twom_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "Rcplex",
                 progress = TRUE)

# Example 2.3 - Using two moments approach and limSolve solver
set.seed(1)
twom_lim <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = reps,
                 tau = tau,
                 solver = "limSolve",
                 progress = TRUE)
