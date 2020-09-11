## ========================================================================= ##
##
##  Example file for invertci function
##
##  This is an example code for applying the invertci procedure on the missing
##  data problem using the sample data. This file illustrates how the module
##  can be used to obtain the confidence interval
##  for the subsampling procedure problem.
##
## ========================================================================= ##
rm(list = ls())

# ---------------- #
# Part 1: Load required packages
# ---------------- #
library(lpinfer)
library(future)

# ---------------- #
# Part 2: Data, lpmodel preparation and arguments for estbounds
# ---------------- #
source("./example/dgp_missingdata.R")
J <- 5
N <- 1000
data <- missingdata_draw(J = J, n = N, seed = 1, prob.obs = .5)
lpmodel.full <- missingdata_lpm(J = J, info = "full", data = data)

tau <- sqrt(log(N)/N)
beta.tgt <- .2
reps <- 100
phi <- 2/3


# Define the arguments
farg <- list(data = data,
             lpmodel = lpmodel.full,
             R = reps,
             phi = phi,
             solver = "gurobi",
             progress = FALSE)

# Example 1: Construction of one confidence interval
set.seed(1)
invertci1 <- invertci(f = subsample,
                      farg = farg,
                      alpha = 0.05,
                      lb0 = 0,
                      lb1 = .4,
                      ub0 = 1,
                      ub1 = .6,
                      tol = 0.001,
                      max.iter = 50,
                      df_ci = NULL,
                      progress = FALSE)
print(invertci1)
summary(invertci1)

# Example 2: Construct a list of multiple confidence intervals
set.seed(1)
invertci2 <- invertci(f = subsample,
                      farg = farg,
                      alpha = c(0.05, 0.1, 0.2),
                      lb0 = 0,
                      lb1 = .4,
                      ub0 = 1,
                      ub1 = .6,
                      tol = 0.001,
                      max.iter = 5,
                      df_ci = NULL,
                      progress = FALSE)
print(invertci2)
summary(invertci2)

# Example 3: Print only the list of selected output
summary(invertci2, alphas = .05)
