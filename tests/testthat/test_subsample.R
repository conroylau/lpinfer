context("Tests for subsample")

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)

# ---------------- #
# Define functions to match the moments
# ---------------- #
## 1. Full information approach
# Function for Omega_hat
var_full_info <- function(data){
  len <- length(unique(data[,"Y"]))
  return(diag(len))
}

# Function for beta_hat
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
  # Variance
  var <- var_full_info(data)
  return(list(beta = beta,
              var = var))
}

## 2. Two moments approach
# Function for Omega_hat
var_two_moment <- function(data){
  return(diag(2))
}

# Function for beta_hat
func_two_moment <- function(data){
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
# Data preparation and arguments for the subsample function
# ---------------- #
# Declare parameters
N <- dim(sampledata)[1]
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])
reps <- 100
# Compute matrices
yp <- seq(0, 1, 1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)
A_obs_full <- cbind(matrix(rep(0, J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0, J1), yp, rep(0, J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)
# Define the target beta and the value of tau
tau <- sqrt(log(N)/N)
beta.tgt <- .365

# ---------------- #
# Data preparation and arguments for the subsample function
# ---------------- #
# Declare parameters
N <- dim(sampledata)[1]
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])
reps <- 100
# Compute matrices
yp <- seq(0, 1, 1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)
A_obs_full <- cbind(matrix(rep(0, J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0, J1), yp, rep(0, J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)
# Shape constraints
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom <- matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp <- c(1)
# Define the target beta and the value of tau
tau <- sqrt(log(N)/N)
beta.tgt <- .365
phi <- 2/3
# Define arguments for the `subsample` function
farg <- list(data = sampledata,
             R = 100,
             beta.tgt = beta.tgt,
             cores = 1,
             norm = 2,
             phi = phi,
             replace = FALSE,
             progress = TRUE)

# ---------------- #
# Define the lpmodel objects
# ---------------- #
# 1. Full information approach
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_full,
                        beta.obs = func_full_info,
                        beta.shp = beta_shp)

# 2. Two moments approach
lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_twom,
                        beta.obs = func_two_moment,
                        beta.shp = beta_shp)

# ---------------- #
# Run results by different solvers with 1 core
# ---------------- #
### 1. Full information approach
farg$lpmodel <- lpmodel.full
## L2-norm
farg$norm <- 2

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g2 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r2 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l2 <- do.call(subsample, farg)

## L1-norm
farg$norm <- 1

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g1 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r1 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l1 <- do.call(subsample, farg)

### 2. Two moments approach
farg$lpmodel <- lpmodel.twom
## L2-norm
farg$norm <- 2

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g2 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r2 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l2 <- do.call(subsample, farg)

## L1-norm
farg$norm <- 1

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g1 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r1 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l1 <- do.call(subsample, farg)

# ---------------- #
# Run results by different solvers with 8 cores
# ---------------- #
farg$cores <- 8
### 1. Full information approach
farg$lpmodel <- lpmodel.full
## L2-norm
farg$norm <- 2

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g28 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r28 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l28 <- do.call(subsample, farg)

## L1-norm
farg$norm <- 1

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g18 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r18 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l18 <- do.call(subsample, farg)

### 2. Two moments approach
farg$lpmodel <- lpmodel.twom
## L2-norm
farg$norm <- 2

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g28 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r28 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l28 <- do.call(subsample, farg)

## L1-norm
farg$norm <- 1

# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g18 <- do.call(subsample, farg)
# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r18 <- do.call(subsample, farg)
# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l18 <- do.call(subsample, farg)

# ---------------- #
# Consolidate the list of outputs
# ---------------- #
ss.output <- list()
ss.output[[1]] <- full_g2
ss.output[[2]] <- full_r2
ss.output[[3]] <- full_l2
ss.output[[4]] <- twom_g2
ss.output[[5]] <- twom_r2
ss.output[[6]] <- twom_l2
ss.output[[7]] <- full_g1
ss.output[[8]] <- full_r1
ss.output[[9]] <- full_l1
ss.output[[10]] <- twom_g1
ss.output[[11]] <- twom_r1
ss.output[[12]] <- twom_l1
ss.output[[13]] <- full_g28
ss.output[[14]] <- full_r28
ss.output[[15]] <- full_l28
ss.output[[16]] <- twom_g28
ss.output[[17]] <- twom_r28
ss.output[[18]] <- twom_l28
ss.output[[19]] <- full_g18
ss.output[[20]] <- full_r18
ss.output[[21]] <- full_l18
ss.output[[22]] <- twom_g18
ss.output[[23]] <- twom_r18
ss.output[[24]] <- twom_l18

# ---------------- #
# Create Gurobi solver
# ---------------- #
gurobi.qlp <- function(Q = NULL, obj, objcon, A, rhs, quadcon = NULL, sense,
                       modelsense, lb) {
  # Set up the model
  model <- list()
  
  # Objective function
  model$Q <- Q
  model$obj <- obj
  model$objcon <- objcon
  
  # Linear constraints
  model$A <- A
  model$rhs <- rhs
  
  # Quadrtaic constraints
  model$quadcon <- quadcon
  
  # Model sense and lower bound
  model$sense <- sense
  model$modelsense <- modelsense
  model$lb <- lb
  
  # Obtain the results to the optimization problem
  params <- list(OutputFlag = 0, FeasibilityTol = 1e-9)
  result <- gurobi::gurobi(model, params)
  
  return(result)
}

# ---------------- #
# Construct the answer by the programs
# ---------------- #
# 1. Obtain the parameters
m <- floor(N^phi)
nx <- ncol(A_tgt)
ng.full <- nrow(A_obs_full)
ng.twom <- nrow(A_obs_twom)

# 2. Define the LP for the L1-problem and the QP for the L2-problem
## L1-problem
subsample.l1.arg <- function(lpmodel, beta.tgt, data, ng, n) {
  # Obtain beta and omega
  diag.omega <- diag(lpmodel$beta.obs(data)$var)
  g <- 1/diag.omega
  g[diag.omega == 0] <- 0
  G <- diag(g)
  bobs <- lpmodel$beta.obs(data)$beta
  # Consolidate the constraints
  Acomb <- rbind(lpmodel$A.shp, lpmodel$A.tgt)
  Afull1 <- cbind(Acomb, matrix(rep(0, ng*2*nrow(Acomb)), nrow = nrow(Acomb)))
  Afull2 <- cbind(G %*% lpmodel$A.obs, -diag(ng), diag(ng))
  # Formulate the list of arguments
  Tn.arg <- list(Q = NULL,
                 obj = c(rep(0, nx), rep(1, ng), rep(-1, ng)) * sqrt(n),
                 objcon = 0,
                 A = rbind(Afull1, Afull2),
                 rhs = c(lpmodel$beta.shp, beta.tgt, G %*% bobs),
                 sense = "=",
                 modelsense = "min",
                 lb = rep(0, nx + ng*2))
  return(Tn.arg)
}
## L2-problem
subsample.l2.arg <- function(lpmodel, beta.tgt, data, ng, n) {
  # Obtain beta and omega
  diag.omega <- diag(lpmodel$beta.obs(data)$var)
  g <- 1/diag.omega
  g[diag.omega == 0] <- 0
  G <- diag(g)
  bobs <- lpmodel$beta.obs(data)$beta
  # Consolidate the constraints
  GA <- G %*% lpmodel$A.obs
  Gb <- G %*% bobs
  # Formulate the list of arguments
  Tn.arg <- list(Q = n * t(GA) %*% GA,
                 obj = -2 * n %*% t(Gb) %*% GA,
                 objcon = n * t(Gb) %*% Gb,
                 A = rbind(lpmodel$A.shp, lpmodel$A.tgt),
                 rhs = c(lpmodel$beta.shp, beta.tgt),
                 sense = "=",
                 modelsense = "min",
                 lb = rep(0, nx))
  return(Tn.arg)
}

# 3. Solve the L1 and L2 problem with the full data
ss.full.Tn <- list()
ss.twom.Tn <- list()
## L1-problem - Full information
ss.l1.full.arg <- subsample.l1.arg(lpmodel.full, beta.tgt, sampledata,
                                   ng.full, N)
ss.full.Tn[[1]] <- do.call(gurobi.qlp, ss.l1.full.arg)
## L1-problem - Two moments
ss.l1.twom.arg <- subsample.l1.arg(lpmodel.twom, beta.tgt, sampledata,
                                   ng.twom, N)
ss.twom.Tn[[1]]  <- do.call(gurobi.qlp, ss.l1.twom.arg)
## L2-problem - Full information
ss.l2.full.arg <- subsample.l2.arg(lpmodel.full, beta.tgt, sampledata,
                                   ng.full, N)
ss.full.Tn[[2]] <- do.call(gurobi.qlp, ss.l2.full.arg)
## L2-problem - Two moments
ss.l2.twom.arg <- subsample.l2.arg(lpmodel.twom, beta.tgt, sampledata,
                                   ng.twom, N)
ss.twom.Tn[[2]] <- do.call(gurobi.qlp, ss.l2.twom.arg)

# 4. Bootstrap procedure
ss.full.Tbs <- list()
ss.full.Tbs[[1]] <- list()
ss.full.Tbs[[2]] <- list()
ss.twom.Tbs <- list()
ss.twom.Tbs[[1]] <- list()
ss.twom.Tbs[[2]] <- list()
ss.full.Tn.temp <- list()
ss.twom.Tn.temp <- list()

set.seed(1)
for (i in 1:reps) {
  ## Construct tau-tightened recentered bootstrap estimates
  data.bs <- as.data.frame(sampledata[sample(1:nrow(sampledata),
                                             size = m,
                                             replace = FALSE),])
  ## L1-problem - Full information
  ss.l1.full.arg <- subsample.l1.arg(lpmodel.full, beta.tgt, data.bs,
                                     ng.full, m)
  ss.full.Tn.temp[[1]] <- do.call(gurobi.qlp, ss.l1.full.arg)
  ## L1-problem - Two moments
  ss.l1.twom.arg <- subsample.l1.arg(lpmodel.twom, beta.tgt, data.bs,
                                     ng.twom, m)
  ss.twom.Tn.temp[[1]]  <- do.call(gurobi.qlp, ss.l1.twom.arg)
  ## L2-problem - Full information
  ss.l2.full.arg <- subsample.l2.arg(lpmodel.full, beta.tgt, data.bs,
                                     ng.full, m)
  ss.full.Tn.temp[[2]] <- do.call(gurobi.qlp, ss.l2.full.arg)
  ## L2-problem - Two moments
  ss.l2.twom.arg <- subsample.l2.arg(lpmodel.twom, beta.tgt, data.bs,
                                     ng.twom, m)
  ss.twom.Tn.temp[[2]] <- do.call(gurobi.qlp, ss.l2.twom.arg)
  ## Consolidate them into a list
  for (j in 1:2) {
    ### Full information
    ss.full.Tbs[[j]][[i]] <- ss.full.Tn.temp[[j]]$objval
    ### Two moments
    ss.twom.Tbs[[j]][[i]] <- ss.twom.Tn.temp[[j]]$objval
  }
}

# 5. Evaluate p-value
pval.full <- list()
pval.twom <- list()
for (j in 1:2) {
  pval.full[[j]] <- mean(unlist(ss.full.Tbs[[j]]) > ss.full.Tn[[j]]$objval)
  pval.twom[[j]] <- mean(unlist(ss.twom.Tbs[[j]]) > ss.twom.Tn[[j]]$objval)
}

# ---------------- #
# Test if the output are equal
# ---------------- #
# Define the indices
full.l2.c1.ind <- 1:3
full.l1.c1.ind <- 7:9
twom.l2.c1.ind <- 4:6
twom.l1.c1.ind <- 10:12
full.l2.c8.ind <- full.l2.c1.ind + 12
full.l1.c8.ind <- full.l1.c1.ind + 12
twom.l2.c8.ind <- twom.l2.c1.ind + 12
twom.l1.c8.ind <- twom.l1.c1.ind + 12

full.l2.ind <- c(full.l2.c1.ind, full.l2.c8.ind)
full.l1.ind <- c(full.l1.c1.ind, full.l1.c8.ind)
twom.l2.ind <- c(twom.l2.c1.ind, twom.l2.c8.ind)
twom.l1.ind <- c(twom.l1.c1.ind, twom.l1.c8.ind)

# 1. Full information approach p-values
test_that("Full information approach p-values",{
  for (i in full.l2.ind) {
    expect_equal(pval.full[[2]], ss.output[[i]]$pval)
  }
  for (i in full.l1.ind) {
    expect_equal(pval.full[[1]], ss.output[[i]]$pval)
  }
})

# 2. Two moments approach p-values
test_that("Two moments approach p-values",{
  for (i in twom.l2.ind) {
    expect_equal(pval.twom[[2]], ss.output[[i]]$pval)
  }
  for (i in twom.l1.ind) {
    expect_equal(pval.twom[[1]], ss.output[[i]]$pval)
  }
})

# 3. Full information approach test statistics
test_that("Full information approach test statistics",{
  for (i in full.l2.ind) {
    expect_lte(abs(ss.full.Tn[[2]]$objval - ss.output[[i]]$T.n), 1e-6)
  }
  for (i in full.l1.ind) {
    expect_lte(abs(ss.full.Tn[[1]]$objval - ss.output[[i]]$T.n), 1e-6)
  }
})


# 4. Two moments approach test statistics
test_that("Two moments approach test statistics",{
  for (i in twom.l2.ind) {
    expect_lte(abs(ss.twom.Tn[[2]]$objval - ss.output[[i]]$T.n), 1e-6)
  }
  for (i in twom.l1.ind) {
    expect_lte(abs(ss.twom.Tn[[1]]$objval - ss.output[[i]]$T.n), 1e-6)
  }
})

# 5. Solver name
test_that("Solver name",{
  for (i in seq(1, 24, 3)) {
    expect_equal("gurobi", ss.output[[i]]$solver)
  }
  for (i in seq(2, 24, 3)) {
    expect_equal("Rcplex", ss.output[[i]]$solver)
  }
  for (i in seq(3, 24, 3)) {
    expect_equal("limSolve", ss.output[[i]]$solver)
  }
})

# 6. Cores
test_that("Cores",{
  for (i in 1:12) {
    expect_equal(1, ss.output[[i]]$cores)
  }
  for (i in 13:24) {
    expect_equal(8, ss.output[[i]]$cores)
  }
})

# 7. cv.table
test_that("cv.table",{
  k <- length(ss.full.Tbs[[1]])
  full99 <- list()
  full95 <- list()
  full90 <- list()
  twom99 <- list()
  twom95 <- list()
  twom90 <- list()
  for (j in 1:2) {
    full99[[j]] <- sort(unlist(ss.full.Tbs[[j]]))[ceiling(.99*k)]
    full95[[j]] <- sort(unlist(ss.full.Tbs[[j]]))[ceiling(.95*k)]
    full90[[j]] <- sort(unlist(ss.full.Tbs[[j]]))[ceiling(.90*k)]
    twom99[[j]] <- sort(unlist(ss.twom.Tbs[[j]]))[ceiling(.99*k)]
    twom95[[j]] <- sort(unlist(ss.twom.Tbs[[j]]))[ceiling(.95*k)]
    twom90[[j]] <- sort(unlist(ss.twom.Tbs[[j]]))[ceiling(.90*k)]
  }
  for (i in full.l2.ind) {
    j <- 2
    expect_lte(abs(full99[[j]] - ss.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(full95[[j]] - ss.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(full90[[j]] - ss.output[[i]]$cv.table[4,2]), 1e-5)
  }
  for (i in full.l1.ind) {
    j <- 1
    expect_lte(abs(full99[[j]] - ss.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(full95[[j]] - ss.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(full90[[j]] - ss.output[[i]]$cv.table[4,2]), 1e-5)
  }
  for (i in twom.l2.ind) {
    j <- 2
    expect_lte(abs(twom99[[j]] - ss.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(twom95[[j]] - ss.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(twom90[[j]] - ss.output[[i]]$cv.table[4,2]), 1e-5)
  }
  for (i in twom.l1.ind) {
    j <- 1
    expect_lte(abs(twom99[[j]] - ss.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(twom95[[j]] - ss.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(twom90[[j]] - ss.output[[i]]$cv.table[4,2]), 1e-5)
  }
})

# 8. Phi
test_that("Phi",{
  for (i in 1:24) {
    expect_equal(2/3, ss.output[[i]]$phi)
  }
})

# 9. Norm
test_that("Norm",{
  for (i in c(full.l2.ind, twom.l2.ind)) {
    expect_equal(2, ss.output[[i]]$norm)
  }
  for (i in c(full.l1.ind, twom.l1.ind)) {
    expect_equal(1, ss.output[[i]]$norm)
  }
})

# 10. Subsample size
test_that("subsample size",{
  for (i in 1:24) {
    expect_equal(floor(nrow(sampledata)^phi), ss.output[[i]]$subsample.size)
  }
})

# 11. Test logical
test_that("test logical",{
  for (i in 1:24) {
    expect_equal(1, ss.output[[i]]$test.logical)
  }
})

# 12. df.error
test_that("Table for problematic bootstrap replications",{
  for (i in 1:24) {
    expect_equal(0, nrow(ss.output[[i]]$df.error))
  }
})

# 13. Number of successful bootstrap replications
test_that("Number of successful bootstrap replications",{
  for (i in 1:24) {
    expect_equal(100, ss.output[[i]]$R.succ)
  }
})
