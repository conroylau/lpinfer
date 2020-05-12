## ========================================================================= ##
##
##  Example file for dkqs function
##
##  This is an example code for applying the R module dkqs on the missing
##  data problem using the sample data by Torgovitsky (2019). This file
##  illustrates how the module can be used to obtain the p-values using the
##  full-information and two moments approach, and different solvers in R that
##  can solve both linear and quadratic programs. Currently, the solvers
##  supported are as follows:
##
##    - gurobi
##    - cplexAPI (QP in progress)
##    - Rcplex
##    - limsolve
##
## ========================================================================= ##

# ---------------- #
# Step 1: Load required packages
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
library(PtProcess)
library(doSNOW)

# ---------------- #
# Step 2: Data preparation
# ---------------- #
# Read data
data <- read.csv("./data/sampledata.csv")

# Compute parameters required
N <- dim(data)[1]
J <- length(unique(data[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(data[,"D"])

# Compute matrices required
yp <- seq(0,1,1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)

# Define the value of tau to be used
tau <- sqrt(log(N)/N)

# ---------------- #
# Step 3: Define functions to compute beta_obs_hat
# ---------------- #
# Full information approach
func_full_info <- function(data){
  # Initialize beta
  beta <- NULL
  # Find the unique elements of Y, sorting in ascending order
  y_list <- sort(unique(data[,"Y"]))
  # Count total number of rows of data and y_list
  n <- dim(data)[1]
  yn <- length(y_list)
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
# Step 4: Compute A_obs based on the two functions in part 3
# ---------------- #
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)

# ---------------- #
# Step 5: Run the dkqs module to compute p-values
# ---------------- #
### Define the value of beta_tgt
beta.tgt <- .375

# Define the value of full information method
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        beta.obs = func_full_info)

# Define the value of full two moments method
lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        beta.obs = func_two_moment)

### Full information approach
# Example 1.1a - Using full information approach and gurobi solver (1 core)
set.seed(1)
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "gurobi",
                 cores = 1,
                 progress = TRUE)

# Example 1.1b - Using full information approach and gurobi solver (8 cores)
set.seed(1)
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "gurobi",
                 cores = 8,
                 progress = TRUE)

# Example 2.1a - Using full information approach and Rcplex solver (1 core)
set.seed(1)
full_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "Rcplex",
                 cores = 1,
                 progress = TRUE)

# Example 1.2b - Using full information approach and Rcplex solver (8 cores)
set.seed(1)
full_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "Rcplex",
                 cores = 8,
                 progress = TRUE)

# Example 1.3a - Using full information approach and limSolve solver (1 core)
set.seed(1)
full_lim <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "limSolve",
                 cores = 1,
                 progress = TRUE)

# Example 1.3b - Using full information approach and limSolve solver (8 cores)
set.seed(1)
full_lim <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "limSolve",
                 cores = 8,
                 progress = TRUE)

### Two moments approach
# Example 2.1a - Using two moments approach and gurobi solver (1 core)
set.seed(1)
twom_gur <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "gurobi",
                 cores = 1,
                 progress = TRUE)

# Example 2.1b - Using two moments approach and gurobi solver (8 cores)
set.seed(1)
twom_gur <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "gurobi",
                 cores = 8,
                 progress = TRUE)

# Example 2.2a - Using two moments approach and Rcplex solver (1 core)
set.seed(1)
twom_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "Rcplex",
                 cores = 1,
                 progress = TRUE)

# Example 2.2b - Using two moments approach and Rcplex solver (8 cores)
set.seed(1)
twom_rcp <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "Rcplex",
                 cores = 8,
                 progress = TRUE)

# Example 2.3a - Using two moments approach and limSolve solver (1 core)
set.seed(1)
twom_lim <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "limSolve",
                 cores = 1,
                 progress = TRUE)

# Example 2.3b - Using two moments approach and limSolve solver (8 cores)
set.seed(1)
twom_lim <- dkqs(data = data,
                 lpmodel = lpmodel.twom,
                 beta.tgt = beta.tgt,
                 R = 100,
                 tau = tau,
                 solver = "limSolve",
                 cores = 8,
                 progress = TRUE)
