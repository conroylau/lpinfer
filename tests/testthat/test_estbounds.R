context("Tests for estbounds and mincriterion")

# ---------------- #
# Load packages
# ---------------- #
library(lpinfer)

# ---------------- #
# Define functions to match the moments
# ---------------- #
### Full information approach
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
  n <- dim(sampledata)[1]
  yn <- length(y_list)
  # Generate each entry of beta_obs
  for (i in 1:yn) {
    beta_i <- sum((data[,"Y"] == y_list[i]) * (data[,"D"] == 1))/n
    beta <- c(beta,c(beta_i))
  }
  beta <- as.matrix(beta)
  # Variance
  var <- var_full_info(sampledata)
  return(list(beta = beta,
              var = var))
}

### Two moments approach
# Function for Omega_hat
var_two_moment <- function(data){
  return(diag(2))
}

# Function for beta_hat
func_two_moment <- function(data){
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)
  # Count total number of rows of data and y_list
  n <- dim(sampledata)[1]
  # Moment 1 E[YD]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  # Moment 2 E[D]
  beta[2] <- sum(data[,"D"])/n
  # Variance
  var <- var_two_moment(sampledata)
  return(list(beta = beta,
              var = var))
}

# ---------------- #
# Data preparation and declaration
# ---------------- #
# Declare parameters
N <- dim(sampledata)[1]
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])

# Compute matrices required
yp <- seq(0,1,1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)

# Define the observed matrix for each appraoch
A_obs_full <- cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
A_obs_twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)

# ---------------- #
# Shape constraints
# ---------------- #
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
A_shp_twom <- matrix(rep(1, ncol(A_obs_twom)), nrow = 1)
beta_shp <- c(1)

# ---------------- #
# Define arguments and produce output
# ---------------- #
# Parameters to test
beta.tgt <- .365
kap <- 1e-5

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

i.lpmodel <- list(lpmodel.full, lpmodel.twom)
ni <- length(i.lpmodel)

# Define arguments
farg <- list(data = sampledata,
             solver = "gurobi")

# ---------------- #
# Generate output from MINCRITERION
# i: lpmodel/approach, j: norm
# ---------------- #
mincriterion.out <- list()
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel[[i]]
  mincriterion.out[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    mincriterion.out[[i]][[j]] <- do.call(mincriterion, farg)
  }
}

# ---------------- #
# Generate output from ESTBOUNDS
# i: lpmodel/approach, j: norm
# ---------------- #
farg$kappa <- kap
farg$estimate <- TRUE
farg$progress <- TRUE
estbounds.out <- list()
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel[[i]]
  estbounds.out[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    estbounds.out[[i]][[j]] <- do.call(estbounds, farg)
  }
}

# ---------------- #
# Create Gurobi solver
# ---------------- #
gurobi.qlp <- function(Q = NULL, obj, objcon, A, rhs,
                       qc = NULL, sense, modelsense, lb) {
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
  model$quadcon <- qc
  
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
# 1. Create arguments 
## Arguments function for step 1 of the 1-norm problem
estb.11.args <- function(lpmodel) {
  Aobs <- lpmodel$A.obs
  Ashp <- lpmodel$A.shp
  bobs <- lpmodel$beta.obs(sampledata)$beta
  bshp <- lpmodel$beta.shp
  p <- length(bobs)
  args <- list(obj = c(rep(0, ncol(Aobs)), rep(1, p * 2)),
               objcon = 0,
               A = rbind(cbind(Ashp,
                               matrix(0, nrow = nrow(Ashp), ncol = 2 * p)),
                         cbind(Aobs, -diag(p), diag(p))),
               rhs = c(bshp, bobs),
               sense = "=",
               modelsense = "min",
               lb = rep(0, ncol(Aobs) + 2 * p))
  return(args)
}
## Arguments function for step 2 of the 1-norm problem
estb.12.args <- function(lpmodel, Qhat, kappa, modelsense) {
  Aobs <- lpmodel$A.obs
  Ashp <- lpmodel$A.shp
  Atgt <- lpmodel$A.tgt
  bobs <- lpmodel$beta.obs(sampledata)$beta
  bshp <- lpmodel$beta.shp
  p <- length(bobs)
  args <- list(obj = cbind(Atgt, matrix(0, nrow = nrow(Ashp), ncol = 2 * p)),
               objcon = 0,
               A = rbind(cbind(Ashp,
                               matrix(0, nrow = nrow(Ashp), ncol = 2 * p)),
                         cbind(Aobs, -diag(p), diag(p)),
                         c(rep(0, ncol(Aobs)), rep(1, p * 2))),
               rhs = c(bshp, bobs, Qhat * (1 + kappa)),
               sense = c(rep("=", p + length(bshp)), "<="),
               modelsense = modelsense,
               lb = rep(0, ncol(Aobs) + 2 * p))
  return(args)
}
## Arguments function for step 1 of the 2-norm problem
estb.21.args <- function(lpmodel) {
  Aobs <- lpmodel$A.obs
  Ashp <- lpmodel$A.shp
  bobs <- lpmodel$beta.obs(sampledata)$beta
  bshp <- lpmodel$beta.shp
  args <- list(Q = t(Aobs) %*% Aobs,
               obj = -2 * t(Aobs) %*% bobs,
               objcon = t(bobs) %*% bobs,
               A = Ashp,
               rhs = bshp,
               sense = "=",
               modelsense = "min",
               lb = rep(0, ncol(Aobs)))
  return(args)
}
## Arguments function for step 2 of the 2-norm problem
estb.22.args <- function(lpmodel, Qhat, kappa, modelsense) {
  Aobs <- lpmodel$A.obs
  Ashp <- lpmodel$A.shp
  Atgt <- lpmodel$A.tgt
  bobs <- lpmodel$beta.obs(sampledata)$beta
  bshp <- lpmodel$beta.shp
  args <- list(Q = NULL,
               obj = Atgt,
               objcon = 0,
               A = Ashp,
               sense = "=",
               rhs = bshp,
               qc = list(list(Qc = t(Aobs) %*% Aobs,
                         q = -2 * t(Aobs) %*% bobs,
                         rhs = Qhat * (1 + kappa) - t(bobs) %*% bobs,
                         sense = "<=")),
               modelsense = modelsense,
               lb = rep(0, ncol(Aobs)))
  return(args)
}
## Consolidate by steps
estb.step1.args <- function(lpmodel, norm) {
  if (norm == 1) {
    args <- estb.11.args(lpmodel)
  } else if (norm == 2) {
    args <- estb.21.args(lpmodel)
  }
  return(args)
}
## Consolidate by steps
estb.step2.args <- function(lpmodel, Qhat, kappa, modelsense, norm) {
  if (norm == 1) {
    args <- estb.12.args(lpmodel, Qhat, kappa, modelsense)
  } else if (norm == 2) {
    args <- estb.22.args(lpmodel, Qhat, kappa, modelsense)
  }
  return(args)
}

# ---------------- #
# Construct the answer by the programs
# ---------------- #
lb <- list()
ub <- list()
Qhat <- list()
minx <- list()
for (i in 1:ni) {
  lb[[i]] <- list()
  ub[[i]] <- list()
  Qhat[[i]] <- list()
  minx[[i]] <- list()
  for (j in 1:2) {
    arg1 <- estb.step1.args(i.lpmodel[[i]], j)
    step1.return <- do.call(gurobi.qlp, arg1)
    Qhat[[i]][[j]] <- step1.return$objval
    minx[[i]][[j]] <- step1.return$x
    arg2ub <- estb.step2.args(i.lpmodel[[i]], Qhat[[i]][[j]], kap, "max", j)
    arg2lb <- estb.step2.args(i.lpmodel[[i]], Qhat[[i]][[j]], kap, "min", j)
    ub[[i]][[j]] <- do.call(gurobi.qlp, arg2ub)$objval
    lb[[i]][[j]] <- do.call(gurobi.qlp, arg2lb)$objval
  }
}

# ---------------- #
# Test if the output are equal in MINCRITERION
# ---------------- #
# 1. Qhat
test_that("'mincriterion': Q-hat",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_lte(abs(Qhat[[i]][[j]] - mincriterion.out[[i]][[j]]$objval),
                 1e-10)
    }
  }
})

# 2. x
test_that("'mincriterion': X",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(minx[[i]][[j]], mincriterion.out[[i]][[j]]$x)
    }
  }
})

# 3. Norm
test_that("'mincriterion': Norm",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(j, mincriterion.out[[i]][[j]]$norm)
    }
  }
})

# 4. Solver
test_that("'mincriterion': Solver",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal("gurobi", mincriterion.out[[i]][[j]]$solver)
    }
  }
})

# ---------------- #
# Test if the output are equal in ESTBOUNDS
# ---------------- #
# 1. Lower bounds (lb)
test_that("'estbounds': Lower bounds",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(lb[[i]][[j]], estbounds.out[[i]][[j]]$lb)
    }
  }
})

# 2. Upper bounds (ub)
test_that("'estbounds': Upper bounds",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(ub[[i]][[j]], estbounds.out[[i]][[j]]$ub)
    }
  }
})

# 3. Estimate
test_that("'estbounds': Estimate",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(TRUE, estbounds.out[[i]][[j]]$est)
    }
  }
})

# 4. Norm
test_that("'estbounds': Norm",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal(j, estbounds.out[[i]][[j]]$norm)
    }
  }
})

# 5. Solver
test_that("'estbounds': Solver",{
  for (i in 1:ni) {
    for (j in 1:2) {
      expect_equal("gurobi", estbounds.out[[i]][[j]]$solver)
    }
  }
})
