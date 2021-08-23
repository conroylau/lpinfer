## ========================================================================= ##
##
##  Example for the fsst function
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
lambda <- .5
rho <- 1e-4
reps <- 100

# ---------------- #
# Part 3: Run the fsst procedure
# ---------------- #
# Example 1 - Using full information approach and gurobi solver (1 core)
set.seed(1)
plan(multisession, workers = 1)
fsst.full1 <- fsst(data = sampledata,
                   lpmodel = lpmodel.full,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = lambda,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "identity",
                   solver = "gurobi",
                   progress = TRUE)

# Example 2 - Using two moments approach and gurobi solver (1 core)
set.seed(1)
plan(multisession, workers = 1)
fsst.twom1 <- fsst(data = sampledata,
                   lpmodel = lpmodel.twom,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = lambda,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "identity",
                   solver = "gurobi",
                   progress = TRUE)

# Example 3 - Using two moments approach and gurobi solver (1 core) with
# weight.matrix = "diag"
set.seed(1)
plan(multisession, workers = 1)
fsst.twom2 <- fsst(data = sampledata,
                   lpmodel = lpmodel.twom,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = lambda,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "diag",
                   solver = "gurobi",
                   progress = TRUE)

# Example 4 - Using two moments approach and gurobi solver (1 core) with
# weight.matrix = "avar"
set.seed(1)
plan(multisession, workers = 1)
fsst.twom3 <- fsst(data = sampledata,
                   lpmodel = lpmodel.twom,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = lambda,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "avar",
                   solver = "gurobi",
                   progress = TRUE)

# Example 5 - Using full information approach and gurobi solver (1 core)
# with multiple lambdas
set.seed(1)
fsst.full2 <- fsst(data = sampledata,
                   lpmodel = lpmodel.full,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = c(.1, .2, .5),
                   rho = rho,
                   n = nrow(sampledata),
                   weight.matrix = "identity",
                   solver = "gurobi",
                   progress = TRUE)

# Example 6 - Using full information approach and gurobi solver (1 core)
# with data-driven lambda
set.seed(1)
fsst.full3 <- fsst(data = sampledata,
                   lpmodel = lpmodel.full,
                   beta.tgt = beta.tgt,
                   R = reps,
                   lambda = NA,
                   rho = rho,
                   n = nrow(sampledata),
                   weight.matrix = "identity",
                   solver = "gurobi",
                   progress = TRUE)
