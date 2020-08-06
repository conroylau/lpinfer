context("Tests for subsample")
rm(list = ls())

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)
library(future)
library(future.apply)

# ---------------- #
# Define functions to match the moments
# ---------------- #
## 1. Full information approach
func_full_info <- function(data){
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
  var <- diag(length(unique(data[,"Y"])))
  return(list(beta = beta,
              var = var))
}

## 2. Two moments approach
func_two_moment <- function(data){
  # Initialize beta
  beta <- matrix(c(0,0), nrow = 2)

  # Count total number of rows of data and y_list
  n <- dim(data)[1]

  # Computes the two moments E[YD] and E[D]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  beta[2] <- sum(data[,"D"])/n

  # Variance
  var <- diag(length(beta))
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
# Run results by different solvers and cores
# ---------------- #
# List of cores, lpmodel and norm objects to be used
i.cores <- list(1, 8)
j.lpmodel <- list(lpmodel.full, lpmodel.twom)
k.norm <- list(1, 2)
l.solver <- list("gurobi", "Rcplex", "limSolve")

# Generate output
ss.out <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  ss.out[[i]] <- list()
  for (j in seq_along(j.lpmodel)) {
    farg$lpmodel <- j.lpmodel[[j]]
    ss.out[[i]][[j]] <- list()
    for (k in seq_along(k.norm)) {
      farg$norm <- k.norm[[k]]
      ss.out[[i]][[j]][[k]] <- list()
      for (l in seq_along(l.solver)) {
        RNGkind(kind = "L'Ecuyer-CMRG")
        set.seed(1)
        farg$solver <- l.solver[[l]]
        ss.out[[i]][[j]][[k]][[l]] <- do.call(subsample, farg)
      }
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
ngs <- list(ng.full, ng.twom)

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
ss.ts <- list()
for (j in seq_along(j.lpmodel)) {
  ss.ts[[j]] <- list()
  # 1-norm
  ss.l1.rg <- subsample.l1.arg(j.lpmodel[[j]], beta.tgt, sampledata, ngs[[j]],
                               N)
  ss.ts[[j]][[1]] <- do.call(gurobi.qlp, ss.l1.rg)
  # 2-norm
  ss.l2.rg <- subsample.l2.arg(j.lpmodel[[j]], beta.tgt, sampledata, ngs[[j]],
                               N)
  ss.ts[[j]][[2]] <- do.call(gurobi.qlp, ss.l2.rg)
}

# 4. Bootstrap procedure
## Initialize the lists
ss.full.Tbs <- list()
ss.full.Tbs[[1]] <- list()
ss.full.Tbs[[2]] <- list()
ss.twom.Tbs <- list()
ss.twom.Tbs[[1]] <- list()
ss.twom.Tbs[[2]] <- list()
ss.full.Tn.temp <- list()
ss.twom.Tn.temp <- list()

## Draw the bootstrap data
set.seed(1)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- .Random.seed

## Function to compute one bootstrap replication
ss.fn <- function(x, data, lpmodel, m, ng, norm) {
  # Draw bootstrap data
  data.bs <- as.data.frame(sampledata[sample(1:nrow(sampledata),
                                             size = m,
                                             replace = FALSE),])

  # Compute argument for the norm
  if (norm == 1) {
    ss.arg <- subsample.l1.arg(lpmodel, beta.tgt, data.bs, ng, m)
  } else if (norm == 2) {
    ss.arg <- subsample.l2.arg(lpmodel, beta.tgt, data.bs, ng, m)
  }

  # Compute the test statistic
  Treturn <- do.call(gurobi.qlp, ss.arg)
  Ts <- Treturn$objval
  return(Ts)
}

# Compute the bootstrap test statistics
plan(multisession, workers = 8)
ss.bs.ts <- list()
for (j in seq_along(ngs)) {
  ss.bs.ts[[j]] <- list()
  for (k in seq_along(k.norm)) {
    ss.bs.ts[[j]][[k]] <-
      unlist(future.apply::future_lapply(1:reps,
                                         FUN = ss.fn,
                                         future.seed = seed,
                                         data = sampledata,
                                         lpmodel = j.lpmodel[[j]],
                                         m = m,
                                         ng = ngs[[j]],
                                         norm = k.norm[[k]]))
  }
}

# 5. Evaluate p-value
ss.pval <- list()
for (j in seq_along(j.lpmodel)) {
  ss.pval[[j]] <- list()
  for (k in seq_along(k.norm)) {
    ss.pval[[j]][[k]] <- mean(ss.bs.ts[[j]][[k]] > ss.ts[[j]][[k]]$objval)
  }
}

# ---------------- #
# Test if the output are equal
# ---------------- #
# 1. Full information approach p-values
test_that("Full information approach p-values",{
  for (i in seq_along(i.cores)) {
    j <- 1
    for (k in seq_along(k.norm)) {
      for (l in seq_along(l.solver)) {
        expect_lte(abs(ss.pval[[j]][[k]] - ss.out[[i]][[j]][[k]][[l]]$pval),
                   1e-5)
      }
    }
  }
})

# 2. Two moments approach p-values
test_that("Two moments approach p-values",{
  for (i in seq_along(i.cores)) {
    j <- 2
    for (k in seq_along(k.norm)) {
      expect_lte(abs(ss.pval[[j]][[k]] - ss.out[[i]][[j]][[k]][[l]]$pval),
                 1e-5)
    }
  }
})

# 3. Full information approach test statistics
test_that("Full information approach test statistics",{
  for (i in seq_along(i.cores)) {
    j <- 1
    for (k in seq_along(k.norm)) {
      expect_lte(abs(ss.ts[[j]][[k]]$objval - ss.out[[i]][[j]][[k]][[l]]$T.n),
                 1e-5)
    }
  }
})

# 4. Two moments approach test statistics
test_that("Two moments approach test statistics",{
  for (i in seq_along(i.cores)) {
    j <- 2
    for (k in seq_along(k.norm)) {
      expect_lte(abs(ss.ts[[j]][[k]]$objval - ss.out[[i]][[j]][[k]][[l]]$T.n),
                 1e-5)
    }
  }
})

# 5. Solver name
test_that("Solver name",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(l.solver[[l]], ss.out[[i]][[j]][[k]][[l]]$solver)
        }
      }
    }
  }
})

# 6. cv.table
cv <- list()
test_that("cv.table",{
  n99 <- ceiling(.99 * reps)
  n95 <- ceiling(.95 * reps)
  n90 <- ceiling(.90 * reps)
  # Compute the critical values
  for (j in seq_along(j.lpmodel)) {
    cv[[j]] <- list()
    for (k in seq_along(k.norm)) {
      cv[[j]][[k]] <- list()
      cv[[j]][[k]][[1]] <- ss.ts[[j]][[k]]$objval
      cv[[j]][[k]][[2]] <- sort(unlist(ss.bs.ts[[j]][[k]]))[n99]
      cv[[j]][[k]][[3]] <- sort(unlist(ss.bs.ts[[j]][[k]]))[n95]
      cv[[j]][[k]][[4]] <- sort(unlist(ss.bs.ts[[j]][[k]]))[n90]
    }
  }

  # Check the critical values
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_lte(abs(cv[[j]][[k]][[1]] -
                         ss.out[[i]][[j]][[k]][[l]]$cv.table[1,2]),
                     1e-5)
          for (p in 2:4) {
            expect_lte(abs(cv[[j]][[k]][[p]] -
                            ss.out[[i]][[j]][[k]][[l]]$cv.table[p,2]),
                           1e-5)
          }
        }
      }
    }
  }
})

# 7. Phi
test_that("Phi",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(2/3, ss.out[[i]][[j]][[k]][[l]]$phi)
        }
      }
    }
  }
})

# 8. Norm
test_that("Norm",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(k.norm[[k]], ss.out[[i]][[j]][[k]][[l]]$norm)
        }
      }
    }
  }
})

# 9. Subsample size
test_that("subsample size",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(floor(nrow(sampledata)^phi),
                       ss.out[[i]][[j]][[k]][[l]]$subsample.size)
        }
      }
    }
  }
})

# 10. Test logical
test_that("test logical",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(1, ss.out[[i]][[j]][[k]][[l]]$test.logical)
        }
      }
    }
  }
})

# 11. df.error
test_that("Table for problematic bootstrap replications",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(TRUE, is.null(ss.out[[i]][[j]][[k]][[l]]$df.error))
        }
      }
    }
  }
})

# 12. Number of successful bootstrap replications
test_that("Number of successful bootstrap replications",{
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        for (l in seq_along(l.solver)) {
          expect_equal(reps, ss.out[[i]][[j]][[k]][[l]]$R.succ)
        }
      }
    }
  }
})
