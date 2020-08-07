## ========================================================================= ##
##
##  Example file for invertci function
##
##  This is an example code for applying the R module invertci on the missing
##  data problem using the sample data by Torgovitsky (2019). This file
##  illustrates how the module can be used to obtain the confidence interval
##  for the subsampling procedure problem.
##
## ========================================================================= ##

# ---------------- #
# Part 1: Load required packages
# ---------------- #
library(modelr)
library(gurobi)
library(e1071)
library(cplexAPI)
library(Rcplex)
library(ddpcr)
library(limSolve)
library(foreach)
library(doMC)
library(parallel)
library(PtProcess)
library(doSNOW)

# ---------------- #
# Part 2: Data preparation
# ---------------- #
# Read data
data <- read.csv("./data/sampledata.csv")

# Compute parameters required
N <- nrow(data)
J <- length(unique(data[,"Y"])) - 1
J1 <- J + 1

# Compute matrices required
yp <- seq(0,1,1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)

# Define the value of tau to be used
tau <- sqrt(log(N)/N)

# ---------------- #
# Part 3: Define functions to compute beta_obs_hat
# ---------------- #
### Full information approach
# Function for Omega_hat
var_full_info <- function(data) {
  len <- length(unique(data[,"Y"]))
  return(diag(len))
}

# Function for beta_hat
func_full_info <- function(data) {
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
var_two_moment <- function(data) {
  return(diag(2))
}

# Function for beta_hat
func_two_moment <- function(data) {
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
# Part 4: Compute A_obs based on the two functions in part 3
# ---------------- #
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)

# ---------------- #
# Part 5: Introduce shape constraints
# ---------------- #
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom <- matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp <- c(1)

# ---------------- #
# Part 6: Arguments for the subsample function without beta.tgt
# ---------------- #
phi.predefine <- .75

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

# Define the arguments
farg <- list(data = data,
             lpmodel = lpmodel.full,
             R = 100,
             phi = phi.predefine,
             solver = "gurobi",
             cores = 1,
             progress = FALSE)

# Demonstration 1: Construction of confidence interval
set.seed(1)
invertci_output <- invertci(f = subsample,
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
print(invertci_output)
summary(invertci_output)

# Demonstration 2: Construct a list of multiple confidence intervals
set.seed(1)
invertci_output_many1 <- invertci(f = subsample,
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
print(invertci_output_many1)
summary(invertci_output_many1)

# Print only the list of selected output
summary(invertci_output_many1, alphas = .05)
