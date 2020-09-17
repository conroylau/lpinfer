## ========================================================================= ##
##
##  Example for the invertci function
##
##  This followings illustrate how the function can be used to construct
##  confidence intervals for the target parameter under different procedures.
##  The missing data problem is used with the subsampling procedure in
##  constructing the confidence intervals.
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
source("./inst/example/dgp_missingdata.R")
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
                      init.lb = c(0, .4),
                      init.ub = c(.6, 1),
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
                      init.lb = c(0, .4),
                      init.ub = c(.6, 1),
                      tol = 0.001,
                      max.iter = 5,
                      df_ci = NULL,
                      progress = FALSE)
print(invertci2)
summary(invertci2)

# Example 3: Print only the list of selected output
summary(invertci2, alphas = .05)

# Example 4: Construction of one confidence interval by specifying the 
# logical lower and upper bounds
set.seed(1)
invertci4 <- invertci(f = subsample,
                      farg = farg,
                      alpha = 0.05,
                      init.lb = 0,
                      init.ub = 1,
                      tol = 0.001,
                      max.iter = 50,
                      df_ci = NULL,
                      progress = FALSE)
print(invertci4)
summary(invertci4)

# Example 5: Construction of one confidence interval without specifying the 
# logical lower and upper bounds - they are computed automatically inside
# the intertci function
set.seed(1)
invertci5 <- invertci(f = subsample,
                      farg = farg,
                      alpha = 0.05,
                      tol = 0.001,
                      max.iter = 50,
                      df_ci = NULL,
                      progress = FALSE)
print(invertci5)
summary(invertci5)