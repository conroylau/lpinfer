context("Tests for dkqs")

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)
library(future)
library(future.apply)

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
# Output 1: beta.obs is a function
# ---------------- #
# List of cores, lpmodel and norm objects to be used
i.cores <- list(1)
j.lpmodel <- list(lpmodel.full, lpmodel.twom)
k.solver <- list("gurobi", "Rcplex", "limSolve")

# Generate output
dkqs.out <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  dkqs.out[[i]] <- list()
  for (j in seq_along(j.lpmodel)) {
    farg$lpmodel <- j.lpmodel[[j]]
    dkqs.out[[i]][[j]] <- list()
    for (k in seq_along(k.solver)) {
      set.seed(1)
      farg$solver <- k.solver[[k]]
      dkqs.out[[i]][[j]][[k]] <- do.call(dkqs, farg)
    }
  }
}

# ---------------- #
# Output 2: beta.obs is a list that contains the sample and bootstrap estimates
# ---------------- #
# Function to draw bootstrap data
draw.bs.data <- function(x, f, data) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  bobs.bs <- f(data.bs)
  return(bobs.bs)
}

# Draw bootstrap data for the full information and two moments method
set.seed(1)
bobs.bs.full.list <- future.apply::future_lapply(1:reps,
                                                 FUN = draw.bs.data,
                                                 future.seed = TRUE,
                                                 f = func_full_info,
                                                 data = sampledata)
set.seed(1)
bobs.bs.twom.list <- future.apply::future_lapply(1:reps,
                                                 FUN = draw.bs.data,
                                                 future.seed = TRUE,
                                                 f = func_two_moment,
                                                 data = sampledata)

bobs.full.list <- c(list(func_full_info(sampledata)), bobs.bs.full.list)
bobs.twom.list <- c(list(func_two_moment(sampledata)), bobs.bs.twom.list)

# Redefine the 'lpmodel' object with 'beta.obs' being a list
lpmodel.full.list <- lpmodel.full
lpmodel.twom.list <- lpmodel.twom
lpmodel.full.list$beta.obs <- bobs.full.list
lpmodel.twom.list$beta.obs <- bobs.twom.list

# Define the new lpmodel object and the arguments to be passed to the function
j.lpmodel2 <- list(lpmodel.full.list, lpmodel.twom.list)
farg2 <- list(beta.tgt = beta.tgt,
              R = reps,
              tau = tau,
              n = nrow(sampledata),
              progress = TRUE)

# Compute the dkqs output again
dkqs.out2 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  dkqs.out2[[i]] <- list()
  for (j in seq_along(j.lpmodel2)) {
    farg2$lpmodel <- j.lpmodel2[[j]]
    dkqs.out2[[i]][[j]] <- list()
    for (k in seq_along(k.solver)) {
      set.seed(1)
      farg2$solver <- k.solver[[k]]
      dkqs.out2[[i]][[j]][[k]] <- do.call(dkqs, farg2)
    }
  }
}

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

  # Obtain the results to  the optimization problem
  params <- list(OutputFlag = 0, FeasibilityTol = 1e-9)
  result <- gurobi::gurobi(model, params)

  return(result)
}

# ---------------- #
# Obtain the results without using the dkqs function
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
s.star.list <- list(s.full.star, s.twom.star)

Tn <- list(Tn.full, Tn.twom)

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
  data.bs <- as.data.frame(sampledata[sample(1:nrow(sampledata),
                                             replace = TRUE),])
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

## Function to compute one bootstrap replication
dkqs.fn <- function(x, data, lpmodel, beta.obs, s.star) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  beta.obs.bs <- lpmodel$beta.obs(data.bs)
  beta.bar.bs <- beta.obs.bs - beta.obs + s.star
  bs.arg <-  dkqs.5.arg(Aobs = lpmodel$A.obs, betaobs = beta.bar.bs)
  Tn.return <- do.call(gurobi.qlp, bs.arg)
  Ts <- Tn.return$objval
  return(Ts)
}

## Compute the test statistics by 'future'
T.bs <- list()
for (j in seq_along(j.lpmodel)) {
  set.seed(1)
  lpm <- j.lpmodel[[j]]
  T.bs[[j]] <-
    unlist(future.apply::future_lapply(1:reps,
                                       FUN = dkqs.fn,
                                       future.seed = TRUE,
                                       data = sampledata,
                                       lpmodel = lpm,
                                       beta.obs = lpm$beta.obs(sampledata),
                                       s.star = s.star.list[[j]]))
}

# 4. Compute the p-values
pval <- list()
for (j in seq_along(j.lpmodel)) {
  pval[[j]] <- mean(T.bs[[j]] > Tn[[j]])
}

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
# A list of unit tests
# ---------------- #
tests.dkqs <- function(dkqs.out, test.name) {
  # Assign the name
  test.name <- sprintf("'beta.obs' as %s:", test.name)
  
  # 1. Full information approach p-values
  test_that(sprintf("%s Full information approach", test.name), {
    for (i in seq_along(i.cores)) {
      j <- 1
      for (k in seq_along(k.solver)) {
        expect_equal(pval[[j]], dkqs.out[[i]][[j]][[k]]$pval[1, 2])
      }
    }
  })
  
  # 2. Two moments p-values
  test_that(sprintf("%s Two moments approach", test.name), {
    for (i in seq_along(i.cores)) {
      j <- 2
      for (k in seq_along(k.solver)) {
        expect_equal(pval[[j]], dkqs.out[[i]][[j]][[k]]$pval[1, 2])
      }
    }
  })
  
  # 3. The list of feasible taus
  test_that(sprintf("%s Feasible taus", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(tau, dkqs.out[[i]][[j]][[k]]$tau.feasible)
        }
      }
    }
  })
  
  # 4. The list of infeasible taus
  test_that(sprintf("%s Infeasible taus", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(NULL, dkqs.out[[i]][[j]][[k]]$tau.infeasible)
        }
      }
    }
  })
  
  # 5. Maximum feasible tau
  test_that(sprintf("%s Maximum feasible tau", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(taumax, dkqs.out[[i]][[j]][[k]]$tau.max)
        }
      }
    }
  })
  
  # 6. Test statistics
  test_that(sprintf("%s Test statistics", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_lte(abs(Tn[[j]] - dkqs.out[[i]][[j]][[k]]$T.n), 1e-6)
        }
      }
    }
  })
  
  # 7. Test logical lower bound
  test_that(sprintf("%s Logical lower bound", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(theta.down$objval, dkqs.out[[i]][[j]][[k]]$lb0[1,2])
        }
      }
    }
  })
  
  # 8. Test logical upper bound
  test_that(sprintf("%s Logical upper bound", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(theta.up$objval, dkqs.out[[i]][[j]][[k]]$ub0[1,2])
        }
      }
    }
  })
  
  # 9. Solver name
  test_that(sprintf("%s Solver name", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(k.solver[[k]], dkqs.out[[i]][[j]][[k]]$solver)
        }
      }
    }
  })
  
  # 10. cv.table
  test_that(sprintf("%s cv.table", test.name), {
    cv <- list()
    n99 <- ceiling(.99 * reps)
    n95 <- ceiling(.95 * reps)
    n90 <- ceiling(.90 * reps)
    # Compute the critical values
    for (j in seq_along(j.lpmodel)) {
      cv[[j]] <- list()
      cv[[j]][[1]] <- Tn[[j]]
      cv[[j]][[2]] <- sort(T.bs[[j]])[n99]
      cv[[j]][[3]] <- sort(T.bs[[j]])[n95]
      cv[[j]][[4]] <- sort(T.bs[[j]])[n90]
    }
    
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_lte(abs(cv[[j]][[1]] - dkqs.out[[i]][[j]][[k]]$cv.table[1,2]),
                     1e-5)
          expect_lte(abs(cv[[j]][[2]] - dkqs.out[[i]][[j]][[k]]$cv.table[2,2]),
                     1e-5)
          expect_lte(abs(cv[[j]][[3]] - dkqs.out[[i]][[j]][[k]]$cv.table[3,2]),
                     1e-5)
          expect_lte(abs(cv[[j]][[4]] - dkqs.out[[i]][[j]][[k]]$cv.table[4,2]),
                     1e-5)
        }
      }
    }
  })
  
  # 11. Test logical
  test_that(sprintf("%s Test logical", test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(1, dkqs.out[[i]][[j]][[k]]$test.logical)
        }
      }
    }
  })
  
  # 12. df.error
  test_that(sprintf("%s Table for problematic bootstrap replications",
                    test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(NULL, dkqs.out[[i]][[j]][[k]]$df.error)
        }
      }
    }
  })
  
  # 13. Number of successful bootstrap replications
  test_that(sprintf("%s Number of successful bootstrap replications",
                    test.name), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.solver)) {
          expect_equal(reps, dkqs.out[[i]][[j]][[k]]$R.succ)
        }
      }
    }
  })
}

# ---------------- #
# Run the tests
# ---------------- #
# beta.obs is a function
tests.dkqs(dkqs.out, "function")

# beta.obs is a list of bootstrap estimates
tests.dkqs(dkqs.out2, "list")
