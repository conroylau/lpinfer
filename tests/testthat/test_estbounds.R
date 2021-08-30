context("Tests for estbounds and mincriterion")
rm(list = ls())

# ---------------- #
# Load packages
# ---------------- #
library(lpinfer)
library(future)
library(future.apply)

# ---------------- #
# Define functions to match the moments
# ---------------- #
# Full information approach
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
  var <- diag(length(unique(data[,"Y"])))
  return(list(beta = beta,
              var = var))
}

# Two moments approach
func_two_moment <- function(data){
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)

  # Count total number of rows of data and y_list
  n <- dim(sampledata)[1]

  # Computes the two moments E[YD] and E[D]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  beta[2] <- sum(data[,"D"])/n

  # Variance
  var <- diag(2)
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

# Define the lpmodels when beta.obs is a function
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

# Define the lpmodels when beta.obs is deterministic
lpmodel.full2 <- lpmodel.full
lpmodel.full2$beta.obs <- lpmodel.full$beta.obs(sampledata)
lpmodel.twom2 <- lpmodel.twom
lpmodel.twom2$beta.obs <- lpmodel.twom$beta.obs(sampledata)
i.lpmodel2 <- list(lpmodel.full2, lpmodel.twom2)

# Define other testing parameters
ni <- length(i.lpmodel)
k.solver <- list("gurobi", "Rcplex", "limSolve")
nk <- length(k.solver)

# Define arguments
farg <- list(data = sampledata,
             solver = "gurobi")

# ---------------- #
# Output 1a: Generate output from MINCRITERION and beta.obs is a function
# i: lpmodel/approach, j: norm, k: solver
# ---------------- #
mincriterion.out <- list()
plan(multisession, workers = 8)
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel[[i]]
  mincriterion.out[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    mincriterion.out[[i]][[j]] <- list()
    for (k in 1:nk) {
      farg$solver <- k.solver[[k]]
      mincriterion.out[[i]][[j]][[k]] <- do.call(mincriterion, farg)
    }
  }
}

# ---------------- #
# Output 1b: Generate output from ESTBOUNDS and beta.obs is a function
# i: lpmodel/approach, j: norm, k: solver
# ---------------- #
farg$kappa <- kap
farg$estimate <- TRUE
estbounds.out <- list()
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel[[i]]
  estbounds.out[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    estbounds.out[[i]][[j]] <- list()
    if (j == 1) {
      for (k in 1:nk) {
        farg$solver <- k.solver[[k]]
        estbounds.out[[i]][[j]][[k]] <- do.call(estbounds, farg)
      }
    } else if (j == 2) {
      k <- 1
      farg$solver <- k.solver[[k]]
      estbounds.out[[i]][[j]][[k]] <- do.call(estbounds, farg)
    }
  }
}

# ---------------- #
# Output 2a: Generate output from MINCRITERION and beta.obs is deterministic
# i: lpmodel/approach, j: norm, k: solver
# ---------------- #
farg <- list(data = sampledata,
             solver = "gurobi")

mincriterion.out2 <- list()
plan(multisession, workers = 8)
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel2[[i]]
  mincriterion.out2[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    mincriterion.out2[[i]][[j]] <- list()
    for (k in 1:nk) {
      farg$solver <- k.solver[[k]]
      mincriterion.out2[[i]][[j]][[k]] <- do.call(mincriterion, farg)
    }
  }
}

# ---------------- #
# Output 2b: Generate output from ESTBOUNDS and beta.obs is deterministic
# i: lpmodel/approach, j: norm, k: solver
# ---------------- #
farg$kappa <- kap
farg$estimate <- TRUE
estbounds.out2 <- list()
for (i in 1:ni) {
  farg$lpmodel <- i.lpmodel2[[i]]
  estbounds.out2[[i]] <- list()
  for (j in 1:2) {
    farg$norm <- j
    estbounds.out2[[i]][[j]] <- list()
    if (j == 1) {
      for (k in 1:nk) {
        farg$solver <- k.solver[[k]]
        estbounds.out2[[i]][[j]][[k]] <- do.call(estbounds, farg)
      }
    } else if (j == 2) {
      k <- 1
      farg$solver <- k.solver[[k]]
      estbounds.out2[[i]][[j]][[k]] <- do.call(estbounds, farg)
    }
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
                         q = as.vector(-2 * t(Aobs) %*% bobs),
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
# A list of unit tests for MINCRITERION
# ---------------- #
tests.minc <- function(estbounds.out, test.name) {
  # Assign the name
  test.name <- sprintf("'beta.obs' as %s:", test.name)

  # 1. Qhat
  test_that(sprintf("%s 'mincriterion': Q-hat", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          expect_lte(abs(Qhat[[i]][[j]] -
                           mincriterion.out[[i]][[j]][[k]]$objval),
                     1e-10)
        }
      }
    }
  })

  # 2. Norm
  test_that(sprintf("%s 'mincriterion': Norm", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          expect_equal(j, mincriterion.out[[i]][[j]][[k]]$norm)
        }
      }
    }
  })

  # 3. Solver
  test_that(sprintf("%s 'mincriterion': Solver", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          expect_equal(k.solver[[k]], mincriterion.out[[i]][[j]][[k]]$solver)
        }
      }
    }
  })
}

# ---------------- #
# A list of unit tests for ESTBUONDS
# ---------------- #
tests.estb <- function(estbounds.out, test.name) {
  # Assign the name
  test.name <- sprintf("'beta.obs' as %s:", test.name)

  # 1. Lower bounds (lb)
  test_that(sprintf("%s 'estbounds': Lower bounds", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(lb[[i]][[j]], estbounds.out[[i]][[j]][[k]]$lb)
          } else if (j == 2) {
            expect_equal(lb[[i]][[j]], estbounds.out[[i]][[j]][[1]]$lb)
          }
        }
      }
    }
  })

  # 2. Upper bounds (ub)
  test_that(sprintf("%s 'estbounds': Upper bounds", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(ub[[i]][[j]], estbounds.out[[i]][[j]][[k]]$ub)
          } else if (j == 2) {
            expect_equal(ub[[i]][[j]], estbounds.out[[i]][[j]][[1]]$ub)
          }
        }
      }
    }
  })

  # 3. Checking the value of the mincriterion
  test_that(sprintf("%s 'estbounds': Checking the value of the mincriterion",
                    test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(Qhat[[i]][[j]],
                         estbounds.out[[i]][[j]][[k]]$mincriterion)
          } else if (j == 2) {
            expect_equal(Qhat[[i]][[j]],
                         estbounds.out[[i]][[j]][[1]]$mincriterion)
          }
        }
      }
    }
  })

  # 4. Estimate
  test_that(sprintf("%s 'estbounds': Estimate", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(TRUE, estbounds.out[[i]][[j]][[k]]$est)
          } else if (j == 2) {
            expect_equal(TRUE, estbounds.out[[i]][[j]][[1]]$est)
          }
        }
      }
    }
  })

  # 5. Norm
  test_that(sprintf("%s 'estbounds': Norm", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(j, estbounds.out[[i]][[j]][[k]]$norm)
          } else if (j == 2) {
            expect_equal(j, estbounds.out[[i]][[j]][[1]]$norm)
          }
        }
      }
    }
  })

  # 6. Solver
  test_that(sprintf("%s 'estbounds': Solver", test.name), {
    for (i in 1:ni) {
      for (j in 1:2) {
        for (k in 1:nk) {
          if (j == 1) {
            expect_equal(k.solver[[k]], estbounds.out[[i]][[j]][[k]]$solver)
          } else if (j == 2) {
            expect_equal(k.solver[[1]], estbounds.out[[i]][[j]][[1]]$solver)
          }
        }
      }
    }
  })
}

# ---------------- #
# Run the tests for the optimal cases in ESTBOUNDS
# ---------------- #
# beta.obs is a function
tests.minc(mincriterion.out, "function")

# beta.obs is a list of bootstrap estimates
tests.minc(mincriterion.out2, "list")

# ---------------- #
# Run the tests for the optimal cases in ESTBOUNDS
# ---------------- #
# beta.obs is a function
tests.estb(estbounds.out, "function")

# beta.obs is a list of bootstrap estimates
tests.estb(estbounds.out2, "list")

# ---------------- #
# Unit tests on nonoptimal cases
# ---------------- #
Amat <- matrix(c(1,0), nrow = 1)
lpm.inf <- lpmodel(A.obs = Amat,
                   A.shp = Amat,
                   A.tgt = matrix(c(1,1), nrow = 1),
                   beta.obs = 1,
                   beta.shp = 1)

# This should give [1, Inf]
estinf1 <- estbounds(lpmodel = lpm.inf,
                     solver = "gurobi",
                     estimate = FALSE)

test_that("'estinf1' bounds are [1, Inf]", {
  expect_equal(estinf1$lb, 1)
  expect_equal(estinf1$ub, Inf)
})

test_that("Status for 'estinf1' bounds", {
  expect_equal(estinf1$lb.status, "OPTIMAL")
  expect_equal(estinf1$ub.status, "UNBOUNDED")
})

# This setup is infeasible
lpm.inf2 <- lpm.inf
lpm.inf2$beta.shp <- 2

estinf2 <- estbounds(lpmodel = lpm.inf2,
                     solver = "gurobi",
                     estimate = FALSE)

test_that("'estinf2' bounds are [Inf, -Inf]", {
  expect_equal(estinf2$lb, Inf)
  expect_equal(estinf2$ub, -Inf)
})

test_that("Status for 'estinf2' bounds", {
  expect_equal(estinf2$lb.status, "INFEASIBLE")
  expect_equal(estinf2$ub.status, "INFEASIBLE")
})
