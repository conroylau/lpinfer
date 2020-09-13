## ========================================================================= ##
##
##  Example file for the mincriterion function
##
##  This is an example code for applying the mincriterion function on
##  the missing data problem using the sample data by Torgovitsky (2019).
##  This file illustrates how the module can be used to obtain the solution
##  to the first stage problems in estbounds.
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
kappa <- 1e-5

# ---------------- #
# Step 3: Run mincriterion
# ---------------- #
# Example 1 - Running mincriterion with 1-norm and full-information approach
minc1 <- mincriterion(data = data,
                      lpmodel = lpmodel.full,
                      norm = 1,
                      solver = "gurobi")

# Example 2 - Running mincriterion with 2-norm and two-moments approach
minc2 <- mincriterion(data = data,
                      lpmodel = lpmodel.twom,
                      norm = 2,
                      solver = "gurobi")
