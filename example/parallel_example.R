## ========================================================================= ##
##
##  Example file to compare the time for different number of workers for
##  parallelization
##
## ========================================================================= ##

# ---------------- #
# Part 1: Load required packages
# ---------------- #
library(lpinfer)
library(future)

# ---------------- #
# Part 2: Data preparation
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
# Part 3: Define functions to compute beta_obs_hat
# ---------------- #
# Full information approach
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
  return(beta)
}

# Two moments approach
func_two_moment <- function(data) {
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)
  # Count total number of rows of data and y_list
  n = dim(data)[1]
  # Moment 1 E[YD]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  # Moment 2 E[D]
  beta[2] <- sum(data[,"D"])/n
  return(beta)
}

# ---------------- #
# Part 4: Compute A_obs based on the two functions in part 3
# ---------------- #
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))

# ---------------- #
# Part 5: Run the dkqs module to compute p-values
# ---------------- #
# Define the value of full information method
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        beta.obs = func_full_info)

# Define the value of beta_tgt and significant figures needed
beta_tgt <- .365

# Example 1 - Using full information approach and gurobi solver (1 core)
t10 <- Sys.time()
set.seed(1)
plan(multisession, workers = 1)
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 5000,
                 tau = tau,
                 solver = "gurobi",
                 progress = TRUE)

t11 <- Sys.time()
time1 <- t11 - t10

# Example 2 - Using full information approach and gurobi solver (8 cores)
set.seed(1)
workers <- 8
plan(multisession, workers = workers)
t80 <- Sys.time()
full_gur <- dkqs(data = data,
                 lpmodel = lpmodel.full,
                 beta.tgt = beta.tgt,
                 R = 5000,
                 tau = tau,
                 solver = "gurobi",
                 progress = TRUE)
t81 <- Sys.time()
time8 <- t81 - t80

# Print the time used
print(sprintf("Time used with 1 core: %s", time1))
print(sprintf("Time used with %s cores: %s", workers, time8))
