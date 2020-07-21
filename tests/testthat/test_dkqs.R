context("Tests for dkqs")

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)

# ---------------- #
# Define functions to match the moments
# ---------------- #
# 1. Full information approach
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

# 2. Two moments approach
func_two_moment <- function(data) { 
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
# Data preparation and arguments for the dkqs function
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
# Define arguments for the `dkqs` function
farg <- list(data = sampledata,
             beta.tgt = beta.tgt,
             R = reps,
             tau = tau,
             cores = 1,
             progress = TRUE)

# ---------------- #
# Define the lpmodel objects
# ---------------- #
# 1. Full information approach
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        beta.obs = func_full_info)

# 2. Two moments approach
lpmodel.twom <- lpmodel(A.obs    = A_obs_twom,
                        A.tgt    = A_tgt,
                        beta.obs = func_two_moment)

# ---------------- #
# Run results by different solvers with 1 core
# ---------------- #
## 1. Full information approach
farg$lpmodel <- lpmodel.full
# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g <- do.call(dkqs, farg)

# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r <- do.call(dkqs, farg)

# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l <- do.call(dkqs, farg)

## 2. Two moments approach
farg$lpmodel <- lpmodel.twom
# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g = do.call(dkqs, farg)

# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r <- do.call(dkqs, farg)

# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l <- do.call(dkqs, farg)

# ---------------- #
# Run results by different solvers with 8 core
# ---------------- #
farg$cores <- 8
## 1. Full information approach
farg$lpmodel <- lpmodel.full
# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
full_g8 <- do.call(dkqs, farg)

# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
full_r8 <- do.call(dkqs, farg)

# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
full_l8 <- do.call(dkqs, farg)

## 2. Two moments approach
farg$lpmodel <- lpmodel.twom
# (a) Gurobi
farg$solver <- "gurobi"
set.seed(1)
twom_g8 = do.call(dkqs, farg)

# (b) Rcplex
farg$solver <- "rcplex"
set.seed(1)
twom_r8 <- do.call(dkqs, farg)

# (c) limSolve
farg$solver <- "limsolve"
set.seed(1)
twom_l8 <- do.call(dkqs, farg)

# ---------------- #
# Consolidate the list of outputs
# ---------------- #
dkqs.output <- list()
dkqs.output[[1]] <- full_g
dkqs.output[[2]] <- full_r
dkqs.output[[3]] <- full_l
dkqs.output[[4]] <- twom_g
dkqs.output[[5]] <- twom_r
dkqs.output[[6]] <- twom_l
dkqs.output[[7]] <- full_g8
dkqs.output[[8]] <- full_r8
dkqs.output[[9]] <- full_l8
dkqs.output[[10]] <- twom_g8
dkqs.output[[11]] <- twom_r8
dkqs.output[[12]] <- twom_l8

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
# Construct the answer by the linear programs themselves
# ---------------- #
# 1. Obtain the parameters and the index sets
nx <- ncol(A_tgt)
ones <- matrix(rep(1, nx), nrow = 1)
## Obtain theta.up and theta.down
theta.arg <- list(obj = A_tgt,
                  objcon = 0,
                  A = ones,
                  rhs = c(1),
                  sense = "=",
                  lb = rep(0, nx))
theta.arg$modelsense = "max"
theta.up <- do.call(gurobi.qlp, theta.arg)
theta.arg$modelsense = "min"
theta.down <- do.call(gurobi.qlp, theta.arg)
## Construct the index sets
ind.up <- which(A_tgt %in% theta.up$objval)
ind.down <- which(A_tgt %in% theta.down$objval)
ind.0 <- (1:ncol(A_tgt))[-c(ind.up, ind.down)]

# 2. Solve problem (4)
## Define the general list
dkqs.4.arg <- function(Aobs, betaobs) {
  Tn.arg <- list(Q = N * t(Aobs) %*% Aobs,
                 obj = -2 * N * t(Aobs) %*% betaobs,
                 objcon = N * t(betaobs) %*% betaobs,
                 A = rbind(A_tgt, ones),
                 rhs = c(beta.tgt, 1),
                 sense = "=",
                 modelsense = "min",
                 lb = rep(0, nx))
  return(Tn.arg)
}
## Define the general list of problem (5)
dkqs.5.arg <- function(Aobs, betaobs) {
  rhs <- c(beta.tgt, 1)
  sense <- rep("=", length(rhs))
  rhs.up <- (beta.tgt - theta.down$objval) * tau / (length(c(ind.up, ind.0)))
  rhs.down <- (theta.up$objval - beta.tgt) * tau / (length(c(ind.down, ind.0)))
  rhs.0 <- (1 - 
              rhs.up * length(ind.up) / tau -
              rhs.down * length(ind.down) / tau) * tau / length(ind.0)
  lb.new <- rep(0, nx)
  lb.new[ind.up] <- rhs.up
  lb.new[ind.down] <- rhs.down
  lb.new[ind.0] <- rhs.0
  Tn.arg <- list(Q = N * t(Aobs) %*% Aobs,
                 obj = -2 * N * t(Aobs) %*% betaobs,
                 objcon = N * t(betaobs) %*% betaobs,
                 A = rbind(A_tgt, ones),
                 rhs = rhs,
                 sense = sense,
                 modelsense = "min",
                 lb = lb.new)
  return(Tn.arg)
}
## Full information approach
beta.obs.full <- func_full_info(sampledata)
Tn.full.arg <- dkqs.4.arg(Aobs = A_obs_full, betaobs = beta.obs.full)
Tn.return.full <- do.call(gurobi.qlp, Tn.full.arg)
Tn.full <- Tn.return.full$objval
tau.full.arg <- dkqs.5.arg(Aobs = A_obs_full, betaobs = beta.obs.full)
tau.return.full <- do.call(gurobi.qlp, tau.full.arg)
x.full.star <- tau.return.full$x
s.full.star <- A_obs_full %*% x.full.star
## Two moments approach
beta.obs.twom <- func_two_moment(sampledata)
Tn.twom.arg <- dkqs.4.arg(Aobs = A_obs_twom, betaobs = beta.obs.twom)
Tn.return.twom <- do.call(gurobi.qlp, Tn.twom.arg)
Tn.twom <- Tn.return.twom$objval
tau.twom.arg <- dkqs.5.arg(Aobs = A_obs_twom, betaobs = beta.obs.twom)
tau.return.twom <- do.call(gurobi.qlp, tau.twom.arg)
x.twom.star <- tau.return.twom$x
s.twom.star <- A_obs_twom %*% x.twom.star

# 3. Compute the bootstrap estimators using the cone-tightening procedure
# i.e. solving problem (5)
Tbs.full <- NULL
Tbs.twom <- NULL
## Bootstrap replications
set.seed(1)
beta.full.bar <- NULL
beta.twom.bar <- NULL
for (i in 1:reps) {
  ## Construct tau-tightened recentered bootstrap estimates
  data.bs <- as.data.frame(sampledata[sample(1:nrow(sampledata), replace = TRUE),])
  ## Full information approach
  beta.obs.full.bs <- func_full_info(data.bs)
  beta.bar.full.bs <- beta.obs.full.bs - beta.obs.full + s.full.star
  bs.full.arg <-  dkqs.5.arg(Aobs = A_obs_full, betaobs = beta.bar.full.bs)
  Tn.return.full <- do.call(gurobi.qlp, bs.full.arg)
  Tbs.full <- c(Tbs.full, Tn.return.full$objval)
  ## Two moments approach
  beta.obs.twom.bs <- func_two_moment(data.bs)
  beta.bar.twom.bs <- beta.obs.twom.bs - beta.obs.twom + s.twom.star
  bs.twom.arg <-  dkqs.5.arg(Aobs = A_obs_twom, betaobs = beta.bar.twom.bs)
  Tn.return.twom <- do.call(gurobi.qlp, bs.twom.arg)
  Tbs.twom <- c(Tbs.twom, Tn.return.twom$objval)
  ## Append the bootstrap betas
  beta.full.bar <- cbind(beta.full.bar, beta.bar.full.bs)
  beta.twom.bar <- cbind(beta.twom.bar, beta.bar.twom.bs)
}

# 4. Compute the p-values
pval.full <- mean(Tbs.full > Tn.full)
pval.twom <- mean(Tbs.twom > Tn.twom)

# 5. Compute the maximum feasible tau
tauobj <- c(1, rep(0, nx))
tauA <- NULL
for (i in 1:ncol(A_tgt)) {
  tauA <- rbind(tauA, rep(0, length(tauobj)))
  tauA[i, i + 1] <- 1
  if (i %in% ind.up) {
    tauA[i, 1] <- -(beta.tgt - theta.down$objval) / (length(c(ind.up, ind.0)))
  } else if (i %in% ind.down) {
    tauA[i, 1] <- -(theta.up$objval - beta.tgt) / (length(c(ind.down, ind.0)))
  } else {
    tauA[i, 1] <- -(1 - 
                      (theta.up$objval - beta.tgt) / 
                      (length(c(ind.down, ind.0))) * length(ind.down) -
                      (beta.tgt - theta.down$objval) / 
                      (length(c(ind.up, ind.0))) * length(ind.up)) / 
      length(ind.0)
  }
}
tauA <- rbind(tauA, c(0, A_tgt), c(0, ones))
taurhs <- c(rep(0, ncol(A_tgt)), beta.tgt, 1)
tausense <- c(rep(">=", ncol(A_tgt)), rep("=", 2))
taureturn <- gurobi.qlp(Q = NULL, 
                        obj = tauobj,
                        objcon = 0,
                        A = tauA, 
                        rhs = taurhs,
                        quadcon = NULL,
                        sense = tausense,
                        modelsense = "max",
                        lb = rep(0, length(tauobj)))
taumax <- taureturn$objval

# ---------------- #
# Test if the output are equal
# ---------------- #
# 1. Full information approach p-values
test_that("Full information approach",{
  for (i in c(1:3, 7:9)) {
    expect_equal(pval.full, dkqs.output[[i]]$pval$`p-value`)
  }
})

# 2. Two moments p-values
test_that("Two moments approach",{
  for (i in c(4:6, 10,12)) {
    expect_equal(pval.twom, dkqs.output[[i]]$pval$`p-value`)
  }
})

# 3. The list of feasible taus
test_that("Feasible taus",{
  for (i in 1:12) {
    expect_equal(tau, dkqs.output[[i]]$tau.feasible)
  }
})

# 4. The list of infeasible taus
test_that("Infeasible taus",{
  for (i in 1:12) {
    expect_equal(NULL, dkqs.output[[i]]$tau.infeasible)
  }
})

# 5. Maximum feasible tau
test_that("Maximum feasible tau",{
  for (i in 1:12) {
    expect_equal(taumax, dkqs.output[[i]]$tau.max)
  }
})

# 6. Test statistics
test_that("Test statistics",{
  for (i in c(1:3, 7:9)) {
    expect_lte(abs(Tn.full - dkqs.output[[i]]$T.n), 1e-6)
  }
  for (i in c(4:6, 10:12)) {
    expect_lte(abs(Tn.twom - dkqs.output[[i]]$T.n), 1e-6)
  }
})

# 7. Test logical lower bound
test_that("Logical lower bound",{
  for (i in 1:12) {
    expect_equal(theta.down$objval, dkqs.output[[i]]$lb0[1,2])
  }
})

# 8. Test logical upper bound
test_that("Logical upper bound",{
  for (i in 1:12) {
    expect_equal(theta.up$objval, dkqs.output[[i]]$ub0[1,2])
  }
})

# 9. Solver name
test_that("Solver name",{
  for (i in c(1, 4, 7, 10)) {
    expect_equal("gurobi", dkqs.output[[i]]$solver)
  }
  for (i in c(2, 5, 8, 11)) {
    expect_equal("Rcplex", dkqs.output[[i]]$solver)
  }
  for (i in c(3, 6, 9, 12)) {
    expect_equal("limSolve", dkqs.output[[i]]$solver)
  }
})

# 10. Cores
test_that("Cores",{
  for (i in 1:6) {
    expect_equal(1, dkqs.output[[i]]$cores)
  }
  for (i in 7:12) {
    expect_equal(8, dkqs.output[[i]]$cores)
  }
})

# 11. cv.table
test_that("cv.table",{
  full99 <- sort(Tbs.full)[ceiling(.99*length(Tbs.full))]
  full95 <- sort(Tbs.full)[ceiling(.95*length(Tbs.full))]
  full90 <- sort(Tbs.full)[ceiling(.90*length(Tbs.full))]
  twom99 <- sort(Tbs.twom)[ceiling(.99*length(Tbs.twom))]
  twom95 <- sort(Tbs.twom)[ceiling(.95*length(Tbs.twom))]
  twom90 <- sort(Tbs.twom)[ceiling(.90*length(Tbs.twom))]
  for (i in c(1:3, 7:9)) {
    expect_lte(abs(full99 - dkqs.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(full95 - dkqs.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(full90 - dkqs.output[[i]]$cv.table[4,2]), 1e-5)
  }
  for (i in c(4:6, 10:12)) {
    expect_lte(abs(twom99 - dkqs.output[[i]]$cv.table[2,2]), 1e-5)
    expect_lte(abs(twom95 - dkqs.output[[i]]$cv.table[3,2]), 1e-5)
    expect_lte(abs(twom90 - dkqs.output[[i]]$cv.table[4,2]), 1e-5)
  }
})

# 12. Test logical
test_that("test logical",{
  for (i in 1:12) {
    expect_equal(1, dkqs.output[[i]]$test.logical)
  }
})

# 13. df.error
test_that("Table for problematic bootstrap replications",{
  for (i in 1:12) {
    expect_equal(0, nrow(dkqs.output[[i]]$df.error))
  }
})

# 14. Number of successful bootstrap replications
test_that("Number of successful bootstrap replications",{
  for (i in 1:12) {
    expect_equal(100, dkqs.output[[i]]$R.succ)
  }
})
