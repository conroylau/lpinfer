## ========================================================================= ##
##
##  Example file for fsst function
##
##  This is an example code for applying the R module fsst on the missing
##  data problem using the sample data by Torgovitsky (2019). This file
##  illustrates how the module can be used to obtain the p-values.
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
# Part 3: Construct matrices in lpmodel
# ---------------- #
# Extract relevant information from data
N <- nrow(sampledata)
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1

# Construct A.obs
Aobs.full <- cbind(matrix(rep(0, J1 * J1), nrow = J1), diag(1, J1))

# Construct A.tgt
yp <- seq(0, 1, 1/J)
Atgt <- matrix(c(yp, yp), nrow = 1)

# Construct A.shp
Ashp <- matrix(rep(1, ncol(Aobs.full)), nrow = 1)

# ---------------- #
# Part 4: Define functions to compute beta_obs_hat and Omega_hat
# ---------------- #
# Construct beta.obs (a function)
betaobs.fullinfo <- function(data){
  beta <- NULL
  y.list <- sort(unique(data[,"Y"]))
  n <- dim(data)[1]
  yn <- length(y.list)
  for (i in 1:yn) {
    beta.i <- sum((data[,"Y"] == y.list[i]) * (data[,"D"] == 1))/n
    beta <- c(beta, c(beta.i))
  }
  beta <- as.matrix(beta)
  return(list(beta = beta,
              var = diag(yn)))
}

# ---------------- #
# Part 5: Define the lpmodel object
# ---------------- #
lpm.full <- lpmodel(A.obs    = Aobs.full,
                    A.tgt    = Atgt,
                    A.shp    = Ashp,
                    beta.obs = betaobs.fullinfo,
                    beta.shp = 1)

# ---------------- #
# Part 6: Run the fsst module to compute p-value
# ---------------- #
set.seed(1)
fsst.full1 <- fsst(data = sampledata, 
                   lpmodel = lpm.full,
                   beta.tgt = 0.375,
                   R = 100,
                   lambda = 0.5,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "diag",
                   solver = "gurobi",
                   cores = 1,
                   progress = TRUE)
print(fsst.full1)
summary(fsst.full1)

# ---------------- #
# Part 7: Run the fsst module to compute multiple p-values
# ---------------- #
set.seed(1)
fsst.full2 <- fsst(data = sampledata, 
                   lpmodel = lpm.full,
                   beta.tgt = 0.375,
                   R = 100,
                   lambda = c(0.1, 0.2, 0.5),
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "diag",
                   solver = "gurobi",
                   cores = 1,
                   progress = TRUE)
print(fsst.full2)
summary(fsst.full2)
