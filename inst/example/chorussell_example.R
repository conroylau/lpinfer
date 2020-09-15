## ========================================================================= ##
##
##  Example for the chorussell function
##
##  This followings illustrate how the function can be used to compute
##  p-values and construct confidence intervals using the missing data
##  problem.
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
kappa <- 1e-5

# ---------------- #
# Step 3: Run chorussell
# ---------------- #
# Example 1: By default, the function is used to find the p-values,
# i.e. ci = FALSE
set.seed(1)
cr1 <- chorussell(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = kappa,
                  solver = "gurobi")

# Example 2: Sample as example 1 except that the two-moments approach is used
set.seed(1)
cr2 <- chorussell(data = data,
                  lpmodel = lpmodel.twom,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = kappa,
                  solver = "gurobi")

# Example 3: Find the 95% confidence interval with full-information approach
# The default is to set the signficance level as 0.05
set.seed(1)
cr3 <- chorussell(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = kappa,
                  ci = TRUE,
                  solver = "gurobi")

# Example 4: Find the 90% confidence interval with full-information approach
set.seed(1)
cr4 <- chorussell(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = kappa,
                  ci = TRUE,
                  alpha = .1,
                  solver = "gurobi")

# Example 5: Multiple lambda can be specified
set.seed(1)
cr5 <- chorussell(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = c(1e-10, 1e-5, 1e-1),
                  ci = TRUE,
                  alpha = .1,
                  solver = "gurobi")


# Example 6: Uses parallel programming by the future package
plan(multisession, workers = 8)
set.seed(1)
cr6 <- chorussell(data = data,
                  lpmodel = lpmodel.full,
                  beta.tgt = beta.tgt,
                  R = 100,
                  norm = 1,
                  kappa = c(1e-10, 1e-5),
                  solver = "gurobi")

# Example 7: The confidence interval can also be obtained via the `invertci`
# command although it can be done without using the `invertci` command.
set.seed(1)
farg <- list(data = data,
             lpmodel = lpmodel.full,
             beta.tgt = beta.tgt,
             R = 100,
             norm = 1,
             kappa = kappa,
             solver = "gurobi",
             ci = TRUE,
             progress = FALSE)
cr7 <- invertci(f = chorussell,
                farg = farg,
                alpha = .05,
                max.iter = 100,
                tol = 1e-4)

# Example 8: To construct confidence interval for multiple alphas, this can be
# done by specifying multiple alpha in the `alpha` option
set.seed(1)
cr8 <- invertci(f = chorussell,
                farg = farg,
                alpha = c(.01, .1),
                max.iter = 100,
                tol = 1e-5)
