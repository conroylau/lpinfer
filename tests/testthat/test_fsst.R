context("Tests for fsst")

# ---------------- #
# Load packages
# ---------------- #
library(lpinfer)

# =========================================================================== #
# Case 1: d > p
# =========================================================================== #
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
lam1 <- .5
lam2 <- c(.1, .5)
rho <- 1e-4
reps <- 100

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

# Define arguments
farg <- list(data = sampledata,
             beta.tgt = beta.tgt,
             cores = 1,
             lambda = lam1,
             rho = rho,
             n = NULL,
             R = reps,
             weight.matrix = "diag",
             solver = "gurobi",
             progress = TRUE)

# ---------------- #
# Generate output from FSST
# ---------------- #
# Define the parameters
i.cores <- list(1, 8)
j.lpmodel <- list(lpmodel.full, lpmodel.twom)
k.lambdas <- list(lam1, lam2)
set.seed(1)
farg$cores = 1
farg$lpmodel = lpmodel.full
farg$lambda = lam1
fsst_g1 = do.call(fsst, farg)

# Generate output
fsst.out <- list()
for (i in 1:2) {
  farg$cores <- i.cores[[i]]
  fsst.out[[i]] <- list()
  for (j in 1:2) {
    farg$lpmodel <- j.lpmodel[[j]]
    fsst.out[[i]][[j]] <- list()
    for (k in 1:2) {
      set.seed(1)
      farg$lambda <- k.lambdas[[k]]
      fsst.out[[i]][[j]][[k]] <- do.call(fsst, farg)
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
# 1. Obtain the parameters and the matrices
nx <- ncol(A_tgt)
ones <- matrix(rep(1, nx), nrow = 1)
# Asymptotic variance of beta.obs
n.beta <- list(length(lpmodel.full$beta.obs(sampledata)$beta) +
                 length(beta_shp) +
                 length(beta.tgt),
               length(lpmodel.twom$beta.obs(sampledata)$beta) +
                 length(beta_shp) +
                 length(beta.tgt))
sigma.beta.obs <- list()
sigma.beta <- list()
for (j in 1:2) {
  sigma.beta.obs[[j]] <- j.lpmodel[[j]]$beta.obs(sampledata)$var
  matrix.temp <- matrix(0L, nrow = n.beta[[j]], ncol = n.beta[[j]])
  length.temp <- nrow(sigma.beta.obs[[j]])
  matrix.temp[1:length.temp, 1:length.temp] <- sigma.beta.obs[[j]]
  sigma.beta[[j]] <- matrix.temp
}
# Compute p and d
d <- list()
p <- list()
for (j in 1:2) {
  p[[j]] <- n.beta[[j]]
  d[[j]] <- ncol(j.lpmodel[[j]]$A.obs)
}

# 2. Compute beta.obs
## Sample beta.obs
beta.obs <- list()
for (j in 1:2) {
  beta.obs[[j]] <- j.lpmodel[[j]]$beta.obs(sampledata)$beta
}
## Estimator of asymptotic variance of beta.obs
set.seed(1)
beta.bs <- list()
beta.bs[[1]] <- list()
beta.bs[[2]] <- list()
for (i in 1:reps) {
  data.bs <- as.data.frame(sampledata[sample(1:nrow(sampledata),
                                             replace = TRUE),])
  for (j in 1:2) {
    beta.bs[[j]][[i]] <- j.lpmodel[[j]]$beta.obs(data.bs)$beta
  }
}

# 3. Solve problem (3) - with the sample estimates and the bootstrap estimates
## Function to obtain the arguments
fsst.3.arg <- function(lpmodel, bobs, beta.tgt, sigma.mat, weight.matrix) {
  Aobs <- lpmodel$A.obs
  # Solve weighting matrix
  if (weight.matrix == "identity") {
    Xi <- diag(length(bobs))
  } else if (weight.matrix == "diag") {
    Xi <- diag(diag(solve(sigma.mat)))
  } else if (weight.matrix == "avar") {
    Xi <- solve(sigma.mat)
  }
  # Formulate the arguments
  args <- list(Q = t(Aobs) %*% Xi %*% Aobs,
               obj = -2 * t(bobs) %*% Xi %*% Aobs,
               objcon = t(bobs) %*% Xi %*% bobs,
               A = rbind(lpmodel$A.shp, lpmodel$A.tgt),
               rhs = c(lpmodel$beta.shp, beta.tgt),
               sense = "=",
               modelsense = "min",
               lb = rep(-Inf, ncol(lpmodel$A.obs)))
  return(list(args = args,
              Xi = Xi))
}
## Solve problem (3) for full data and the bootstrap beta
x.star <- list()
beta.obs.star <- list()
beta.star <- list()
x.star.bs <- list()
beta.obs.star.bs <- list()
beta.star.bs <- list()
Xi <- list()
for (j in 1:2) {
  ### Full data
  fsst.3.args <- fsst.3.arg(j.lpmodel[[j]], beta.obs[[j]], beta.tgt,
                            sigma.beta.obs[[j]], "diag")
  Xi[[j]] <- fsst.3.args$Xi
  fsst.3.return <- do.call(gurobi.qlp, fsst.3.args$args)
  x.star[[j]] <- fsst.3.return$x
  if (d[[j]] >= p[[j]]) {
    beta.obs.star[[j]] <- beta.obs[[j]]
  } else {
    beta.obs.star[[j]] <- j.lpmodel[[j]]$A.obs %*% x.star[[j]]
  }
  beta.star[[j]] <- c(beta.obs.star[[j]], j.lpmodel[[j]]$beta.shp, beta.tgt)
  ### Bootstrap data
  x.star.bs[[j]] <- list()
  beta.obs.star.bs[[j]] <- list()
  beta.star.bs[[j]] <- list()
  for (i in 1:reps) {
    fsst.3.bs.args <- fsst.3.arg(j.lpmodel[[j]], beta.bs[[j]][[i]], beta.tgt,
                                 sigma.beta.obs[[j]], "diag")
    fsst.3.bs.return <- do.call(gurobi.qlp, fsst.3.bs.args$args)
    x.star.bs[[j]][[i]] <- fsst.3.bs.return$x
    if (d[[j]] >= p[[j]]) {
      beta.obs.star.bs[[j]][[i]] <- beta.bs[[j]][[i]]
    } else {
      beta.obs.star.bs[[j]][[i]] <- j.lpmodel[[j]]$A.obs %*% x.star.bs[[j]][[i]]
    }
    beta.star.bs[[j]][[i]] <- c(beta.obs.star.bs[[j]][[i]],
                                j.lpmodel[[j]]$beta.shp,
                                beta.tgt)
  }
}

# 4. Standardization
## Compute studentization matrix
rhobar <- list()
student.matrix <- list()
omega <- list()
for (j in 1:2) {
  if (d[[j]] >= p[[j]]) {
    student.matrix[[j]] <- sigma.beta[[j]]
  } else {
    student.matrix[[j]] <- matrix(0, nrow = n.beta[[j]], ncol = n.beta[[j]])
    for (i in 1:reps) {
      beta.t <- beta.star.bs[[j]][[i]] - beta.star[[j]]
      student.matrix[[j]] <- student.matrix[[j]] + beta.t %*% t(beta.t)
    }
    student.matrix[[j]] <- N * student.matrix[[j]] / reps
  }
  ## Compute regularization parameter
  rhobar[[j]] <- base::norm(student.matrix[[j]], type = "f") * rho
  ## Compute regularization matrix
  omega[[j]] <- expm::sqrtm(student.matrix[[j]] + rhobar[[j]] * diag(p[[j]]))
}

# 5. Test statistic
## Cone program arguments
fsst.56.args <- function(lpmodel, p, d, beta.star, omega) {
  A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
  ns <- ncol(t(A))
  nb <- length(beta.star)
  nx <- ncol(A)
  args <- list(Q = NULL,
               obj = c(beta.star, rep(0, p * 2)),
               objcon = 0,
               A = rbind(cbind(omega, -diag(p), diag(p)),
                         c(rep(0, p), rep(1, p * 2)),
                         cbind(t(A), matrix(0L, nrow = ncol(A), ncol = p * 2))),
               rhs = c(rep(0, p), 1, rep(0, ncol(A))),
               sense = c(rep("=", p), rep("<=", 1 + ncol(A))),
               modelsense = "max",
               lb = c(rep(-Inf, nb), rep(0, p * 2)))
  if (d < p) {
    args1 <- list(Q = NULL,
                  obj = c(args$obj, rep(0, nx)),
                  objcon = 0,
                  A = rbind(cbind(args$A,
                                  matrix(0L,
                                         nrow = nrow(args$A),
                                         ncol = ncol(A))),
                            cbind(-diag(ns),
                                  matrix(0L, nrow = ns, ncol = p*2),
                                  A)),
                  rhs = c(args$rhs, rep(0, ns)),
                  sense = c(args$sense, rep("=", ns)),
                  modelsense = "max",
                  lb = c(args$lb, rep(-Inf, nx)))
    return(args1)
  } else {
    return(args)
  }
}
## Range program
fsst.range.soln <- function(beta.obs.star, beta.obs, Xi, p, d) {
  if (d >= p) {
    range.n <- 0
  } else {
    range.n <- base::norm(sqrt(N) * expm::sqrtm(Xi) %*%
                            (beta.obs - beta.obs.star),
                          type = "I")
  }
  return(range.n)
}
## Compute the statistics
range.n <- list()
cone.n <- list()
ts <- list()
for (j in 1:2) {
  ### Range
  range.n[[j]] <- fsst.range.soln(beta.obs.star[[j]], beta.obs[[j]], Xi,
                                  p[[j]], d[[j]])
  ### Cone
  cone.args <- fsst.56.args(j.lpmodel[[j]], p[[j]], d[[j]], beta.star[[j]],
                            omega[[j]])
  cone.return <- do.call(gurobi.qlp, cone.args)
  cone.n[[j]] <- sqrt(N) * cone.return$objval
  ### Test statistics
  ts[[j]] <- max(cone.n[[j]], range.n[[j]])
}

# 6. Restricted estimator
## Arguments
fsst.89.args <- function(lpmodel, p, d, beta.star, beta, omega, beta.tgt) {
  A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
  nbobs <- p - length(c(lpmodel$beta.shp, beta.tgt))
  args <- list(Q = NULL,
               obj = c(rep(0, 2 * p + 2 * d + nbobs), 1),
               objcon = 0,
               A = rbind(cbind(sqrt(N) * diag(p),
                               matrix(0L, nrow = p, ncol = d + nbobs),
                               -omega,
                               A,
                               matrix(0L, nrow = p, ncol = 1)),
                         cbind(matrix(0L, nrow = p, ncol = p + d + nbobs),
                               diag(p),
                               matrix(0L, nrow = p, ncol = d),
                               matrix(1L, nrow = p, ncol = 1)),
                         cbind(matrix(0L, nrow = p, ncol = p + d + nbobs),
                               -diag(p),
                               matrix(0L, nrow = p, ncol = d),
                               matrix(1L, nrow = p, ncol = 1)),
                         cbind(-diag(p),
                               A,
                               matrix(0L, nrow = p, ncol = p + d + nbobs + 1)),
                         cbind(diag(p),
                               matrix(0L, nrow = p, ncol = d),
                               rbind(-diag(nbobs),
                                     matrix(0L, nrow = (p - nbobs),
                                            ncol = nbobs)),
                               matrix(0L, nrow = p, ncol = p + d + 1))),
               rhs = c(sqrt(N) * beta, rep(0, p * 3 + nbobs),
                       c(lpmodel$beta.shp, beta.tgt)),
               sense = c(rep("=", p), rep(">=", p * 2), rep("=", 2 * p)),
               modelsense = "min",
               lb = c(rep(-Inf, p), rep(0, d), rep(-Inf, nbobs), rep(-Inf, p),
                      rep(0, d + 1)))
  if (d < p) {
    args1 <- list(Q = NULL,
                  obj = args$obj,
                  objcon = 0,
                  A = rbind(t(A) %*% args$A[1:p,],
                            args$A[(p + 1):nrow(args$A),]),
                  rhs = c(sqrt(N) * t(A) %*% beta.star,
                          args$rhs[(p + 1):nrow(args$A)]),
                  sense = c(rep("=", ncol(A)),
                            args$sense[(p + 1):nrow(args$A)]),
                  modelsense = args$modelsense,
                  lb = args$lb)
    return(args1)
  } else {
    return(args)
  }
}
## Obtain the estimators
beta.r <- list()
for (j in 1:2) {
  beta.r.args <- fsst.89.args(j.lpmodel[[j]], p[[j]], d[[j]], beta.star[[j]],
                              c(beta.obs[[j]],
                                j.lpmodel[[j]]$beta.shp,
                                beta.tgt),
                              omega[[j]], beta.tgt)
  beta.r[[j]] <- do.call(gurobi.qlp, beta.r.args)$x[1:p[[j]]]
}

# 7. Compute bootstrap components
range.n.bs <- list()
cone.n.bs <- list()
for (j in 1:2) {
  cone.n.bs[[j]] <- list()
  for (k in 1:2) {
    cone.n.bs[[j]][[k]] <- list()
    for (kk in 1:length(k.lambdas[[k]])) {
      cone.n.bs[[j]][[k]][[kk]] <- list()
    }
  }
}
ts.bs <- cone.n.bs
for (j in 1:2) {
  range.n.bs[[j]] <- list()
  for (i in 1:reps) {
    b1 <- beta.bs[[j]][[i]] - beta.obs[[j]]
    b2 <- beta.obs.star.bs[[j]][[i]] - beta.obs.star[[j]]
    range.n.bs[[j]][[i]] <- fsst.range.soln(b1, b2, Xi[[j]], p[[j]], d[[j]])
    for (k in 1:2) {
      for (kk in 1:length(k.lambdas[[k]])) {
        beta.res <- beta.star.bs[[j]][[i]] - beta.star[[j]] +
          k.lambdas[[k]][kk] * beta.r[[j]]
        cone.args <- fsst.56.args(j.lpmodel[[j]], p[[j]], d[[j]], beta.res,
                                  omega[[j]])
        cone.n.bs[[j]][[k]][[kk]][[i]] <- sqrt(N) * do.call(gurobi.qlp,
                                                            cone.args)$objval
        ts.bs[[j]][[k]][[kk]][[i]] <- max(cone.n.bs[[j]][[k]][[kk]][[i]],
                                          range.n.bs[[j]][[i]])
      }
    }
  }
}

# 8. Compute pvalues
pval <- list()
for (j in 1:2) {
  pval[[j]] <- list()
  for (k in 1:2) {
    pval[[j]][[k]] <- list()
    for (kk in 1:length(k.lambdas[[k]])) {
      pval[[j]][[k]][[kk]] <- mean(unlist(ts.bs[[j]][[k]][[kk]]) > ts[[j]])
    }
  }
}

# 9. Critical values
cv <- list()
for (j in 1:2) {
  cv[[j]] <- list()
  n99 <- ceiling(.99 * reps)
  n95 <- ceiling(.95 * reps)
  n90 <- ceiling(.90 * reps)
  for (k in 1:2) {
    cv[[j]][[k]] <- list()
    for (kk in 1:length(k.lambdas[[k]])) {
      cv[[j]][[k]][[kk]] <- list()
      # Test statistics
      cv[[j]][[k]][[kk]][[1]] <- ts[[j]]
      cv[[j]][[k]][[kk]][[2]] <- sort(unlist(ts.bs[[j]][[k]][[kk]]))[n99]
      cv[[j]][[k]][[kk]][[3]] <- sort(unlist(ts.bs[[j]][[k]][[kk]]))[n95]
      cv[[j]][[k]][[kk]][[4]] <- sort(unlist(ts.bs[[j]][[k]][[kk]]))[n90]
      # Cone
      cv[[j]][[k]][[kk]][[5]] <- cone.n[[j]]
      cv[[j]][[k]][[kk]][[6]] <- sort(unlist(cone.n.bs[[j]][[k]][[kk]]))[n99]
      cv[[j]][[k]][[kk]][[7]] <- sort(unlist(cone.n.bs[[j]][[k]][[kk]]))[n95]
      cv[[j]][[k]][[kk]][[8]] <- sort(unlist(cone.n.bs[[j]][[k]][[kk]]))[n90]
    }
    # Range
    cv[[j]][[k]][[1]][[9]] <- range.n[[j]]
    cv[[j]][[k]][[1]][[10]] <- sort(unlist(range.n.bs[[j]]))[n99]
    cv[[j]][[k]][[1]][[11]] <- sort(unlist(range.n.bs[[j]]))[n95]
    cv[[j]][[k]][[1]][[12]] <- sort(unlist(range.n.bs[[j]]))[n90]
  }
}

# ---------------- #
# Test if the output are equal - For d > p only
# i: cores, j: lpmodel approach, k: lambdas
# ---------------- #
# 1. Full information approach p-values
test_that("'d > p': Full information approach",{
  for (i in 1:2) {
    j <- 1
    for (k in 1:2) {
      for (kk in 1:length(k.lambdas[[k]])) {
        expect_equal(pval[[j]][[k]][[kk]],
                     fsst.out[[i]][[j]][[k]]$pval$`p-value`[kk])
      }
    }
  }
})

# 2. Two moments approach p-values
test_that("'d > p': Full information approach",{
  for (i in 1:2) {
    j <- 2
    for (k in 1:2) {
      for (kk in 1:length(k.lambdas[[k]])) {
        expect_equal(pval[[j]][[k]][[kk]],
                     fsst.out[[i]][[j]][[k]]$pval$`p-value`[kk])
      }
    }
  }
})

# 3. Full information CV table
test_that("'d > p': Full information CV table",{
  for (i in 1:2) {
    j <- 1
    for (k in 1:2) {
      for (kk in 1:length(k.lambdas[[k]])) {
        for (l in 1:8) {
          expect_lte(abs(cv[[j]][[k]][[kk]][[l]] -
                           fsst.out[[i]][[j]][[k]]$cv.table[l, kk + 2]),
                     1e-5)
        }
      }
      for (l in 9:12) {
        expect_lte(abs(cv[[j]][[k]][[1]][[l]] -
                     fsst.out[[i]][[j]][[k]]$cv.table[l, 3]),
                   1e-5)
      }
    }
  }
})

# 4. Two moments CV table
test_that("'d > p': Two moments CV table",{
  for (i in 1:2) {
    j <- 2
    for (k in 1:2) {
      for (kk in 1:length(k.lambdas[[k]])) {
        for (l in 1:8) {
          expect_lte(abs(cv[[j]][[k]][[kk]][[l]] -
                         fsst.out[[i]][[j]][[k]]$cv.table[l, kk + 2]),
                     1e-5)
        }
      }
      for (l in 9:12) {
        expect_lte(abs(cv[[j]][[k]][[1]][[l]] -
                       fsst.out[[i]][[j]][[k]]$cv.table[l, 3]),
                   1e-5)
      }
    }
  }
})

# 5. Cores
test_that("'d > p': Cores",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(i.cores[[i]], fsst.out[[i]][[j]][[k]]$cores)
      }
    }
  }
})

# 6. Range test statistics
test_that("'d > p': Range test statistics",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(range.n[[j]], fsst.out[[i]][[j]][[k]]$range)
      }
    }
  }
})

# 7. Cone test statistics
test_that("'d > p': Cone test statistics",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(cone.n[[j]], fsst.out[[i]][[j]][[k]]$cone$objval)
      }
    }
  }
})

# 8. Test statistics
test_that("'d > p': Test statistics",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(ts[[j]], fsst.out[[i]][[j]][[k]]$test)
      }
    }
  }
})

# 9. Solver name
test_that("'d > p': Solver name",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal("gurobi", fsst.out[[i]][[j]][[k]]$solver.name)
      }
    }
  }
})

# 10. Rho parameter
test_that("'d > p': Rho parameter",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(rho, fsst.out[[i]][[j]][[k]]$rho)
      }
    }
  }
})

# 11. Regularization parameter
test_that("'d > p': Regularization parameter",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(rhobar[[j]], fsst.out[[i]][[j]][[k]]$rhobar.i)
      }
    }
  }
})

# 12. Method of obtaining the beta.var matrix
test_that("'d > p': Method of obtaining beta.var",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal("function", fsst.out[[i]][[j]][[k]]$beta.var.method)
      }
    }
  }
})

# 13. Test logical
test_that("'d > p': Omega.i matrix",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(omega[[j]], fsst.out[[i]][[j]][[k]]$omega.i)
      }
    }
  }
})

# 14. Test logical
test_that("'d > p': Test logical",{
  for (i in 1:2) {
    for (j in 1:2) {
      for (k in 1:2) {
        expect_equal(1, fsst.out[[i]][[j]][[k]]$test.logical)
      }
    }
  }
})

# =========================================================================== #
# Case 2: d <= p
# =========================================================================== #
# ---------------- #
# Extract information for d <= p case
# ---------------- #
# Load data
load("test_lpm_fsst.RData")
# Set parameters
n <- 1000
reps2 <- 100
rho2 <- .0001
lambda2 <- .5
btgt2 <- .21242552380635624

farg2 <- list(lpmodel = lpm,
              beta.tgt = btgt2,
              R = reps,
              lambda = lam2,
              rho = rho2,
              weight.matrix = "avar",
              solver = "gurobi",
              n = n,
              progress = TRUE)

# ---------------- #
# Compute solution using FSST
# ---------------- #
fsst.out2 <- list()
for (i in 1:2) {
  farg2$cores <- i.cores[[i]]
  fsst.out2[[i]] <- do.call(fsst, farg2)
}

# ---------------- #
# Define arguments and produce output
# ---------------- #
# 1. Extract relevant information
bobs2 <- lpm$beta.obs[[1]][[1]]
sigma.bobs2 <- lpm$beta.obs[[1]][[2]]
p2 <- length(c(lpm$beta.obs[[1]][[1]], lpm$beta.shp, btgt2))
d2 <- ncol(lpm$A.obs)

# 2. Solve problem (3) - with the sample estimates and the bootstrap estimates
## With full data
fsst.3.args2 <- fsst.3.arg(lpm, bobs2, btgt2, sigma.bobs2, "avar")
Xi2 <- fsst.3.args2$Xi
fsst.3.return2 <- do.call(gurobi.qlp, fsst.3.args2$args)
x.star2 <- fsst.3.return2$x
if (d2 >= p2) {
  bobs.star2 <- bobs2
} else {
  bobs.star2 <- lpm$A.obs %*% x.star2
}
beta.star2 <- c(bobs.star2, lpm$beta.shp, btgt2)
## Bootstrap components
x.star2.list <- list()
beta.obs.star.bs2 <- list()
beta.star.bs2 <- list()
for (i in 1:reps2) {
  fsst.3.bs.args2 <- fsst.3.arg(lpm, lpm$beta.obs[[i + 1]], btgt2,
                                sigma.bobs2, "avar")
  fsst.3.bs.return2 <- do.call(gurobi.qlp, fsst.3.bs.args2$args)
  x.star2.list[[i]] <- fsst.3.bs.return2$x
  if (d2 >= p2) {
    beta.obs.star.bs2[[i]] <- lpm$beta.obs[[i + 1]]
  } else {
    beta.obs.star.bs2[[i]] <- lpm$A.obs %*% x.star2.list[[i]]
  }
  beta.star.bs2[[i]] <- c(beta.obs.star.bs2[[i]], lpm$beta.shp, btgt2)
}
sigma2 <- matrix(0, nrow = p2, ncol = p2)
n.bobs2 <- length(bobs2)
sigma2[1:n.bobs2, 1:n.bobs2] <- sigma.bobs2

# 3. Standardization
## Compute studentization matrix
if (d2 >= p2) {
  student.matrix2 <- sigma2
} else {
  student.matrix2 <- matrix(0, nrow = p2, ncol = p2)
}
for (i in 1:reps2) {
    beta.t <- beta.star.bs2[[i]] - beta.star2
    student.matrix2 <- student.matrix2 + beta.t %*% t(beta.t)
}
student.matrix2 <- n * student.matrix2 / reps2
## Compute regularization parameter
rhobar2 <- base::norm(student.matrix2, type = "f") * rho2
## Compute regularization matrix
omega2 <- expm::sqrtm(student.matrix2 + rhobar2 * diag(p2))

# 4. Compute test statistics
## Range
range.n2 <- fsst.range.soln(bobs.star2, bobs2, Xi2, p2, d2)
## Cone
cone.args2 <- fsst.56.args(lpm, p2, d2, beta.star2, omega2)
cone.return2 <- do.call(gurobi.qlp, cone.args2)
cone.n2 <- sqrt(n) * cone.return2$objval
## Test statistics
ts2 <- max(cone.n2, range.n2)

# 5. Restricted estimator
beta.r.args2 <- fsst.89.args(lpm, p2, d2, beta.star2,
                             c(bobs2, lpm$beta.shp, btgt2), omega2, btgt2)
beta.r2 <- do.call(gurobi.qlp, beta.r.args2)$x[1:p2]

# 6. Compute bootstrap components
range.n.bs2 <- list()
cone.n.bs2 <- list()
ts.bs2 <- list()
for (i in 1:reps) {
  ## Range
  b1 <- lpm$beta.obs[[i + 1]] - bobs2
  b2 <- beta.obs.star.bs2[[i]] - bobs.star2
  range.n.bs2[[i]] <- fsst.range.soln(b1, b2, Xi2, p2, d2)
  
  ## Cone
  beta.res2 <- beta.star.bs2[[i]] - beta.star2 + lambda2 * beta.r2
  cone.args2 <- fsst.56.args(lpm, p2, d2, beta.res2, omega2)
  cone.n.bs2[[i]] <- sqrt(n) * do.call(gurobi.qlp, cone.args2)$objval
  ts.bs2[[i]] <- max(cone.n.bs2[[i]], range.n.bs2[[i]])
}

# 7. Compute p-value
pval2 <- mean(unlist(ts.bs2) > ts2)

# 8. Critical values
n99.2 <- ceiling(.99 * reps2)
n95.2 <- ceiling(.95 * reps2)
n90.2 <- ceiling(.90 * reps2)
## Fill in the critical values
cv2 <- c(ts2,
         sort(unlist(ts.bs2))[n99.2],
         sort(unlist(ts.bs2))[n95.2],
         sort(unlist(ts.bs2))[n90.2],
         cone.n2,
         sort(unlist(cone.n.bs2))[n99.2],
         sort(unlist(cone.n.bs2))[n95.2],
         sort(unlist(cone.n.bs2))[n90.2],
         range.n2,
         sort(unlist(range.n.bs2))[n99.2],
         sort(unlist(range.n.bs2))[n95.2],
         sort(unlist(range.n.bs2))[n90.2])

# ---------------- #
# Test if the output are equal - For d <= p only
# i: cores
# ---------------- #
# 1. p-value
test_that("'d <= p': p-value",{
  for (i in 1:2) {
    expect_equal(pval2, fsst.out2[[i]]$pval$`p-value`[1])
  }
})

# 2. Full information CV table
test_that("'d <= p': CV table",{
  for (i in 1:2) {
    for (l in 1:12) {
      expect_lte(abs(cv2[l] - fsst.out2[[i]]$cv.table[l, 3]), 1e-5)
    }
  }
})

# 3. Cores
test_that("'d <= p': Cores",{
  for (i in 1:2) {
    expect_equal(i.cores[[i]], fsst.out2[[i]]$cores)
  }
})

# 4. Range test statistics
test_that("'d <= p': Range test statistics",{
  for (i in 1:2) {
    expect_equal(range.n2, fsst.out2[[i]]$range)
  }
})

# 5. Cone test statistics
test_that("'d <= p': Cone test statistics",{
  for (i in 1:2) {
    expect_equal(cone.n2, fsst.out2[[i]]$cone$objval)
  }
})

# 6. Test statistics
test_that("'d <= p': Test statistics",{
  for (i in 1:2) {
    expect_equal(ts2, fsst.out2[[i]]$test)
  }
})

# 7. Solver name
test_that("'d <= p': Solver name",{
  for (i in 1:2) {
    expect_equal("gurobi", fsst.out2[[i]]$solver.name)
  }
})

# 8. Rho parameter
test_that("'d <= p': Rho parameter",{
  for (i in 1:2) {
    expect_equal(rho2, fsst.out2[[i]]$rho)
  }
})

# 9. Regularization parameter
test_that("'d <= p': Regularization parameter",{
  for (i in 1:2) {
    expect_equal(rhobar2, fsst.out2[[i]]$rhobar.i)
  }
})

# 10. Method of obtaining the beta.var matrix
test_that("'d <= p': Method of obtaining beta.var",{
  for (i in 1:2) {
    expect_equal("list", fsst.out2[[i]]$beta.var.method)
  }
})

# 11. Test logical
test_that("'d <= p': Omega.i matrix",{
  for (i in 1:2) {
    expect_equal(omega2, fsst.out2[[i]]$omega.i)
  }
})

# 12. Test logical
test_that("'d <= p': Test logical",{
  for (i in 1:2) {
    expect_equal(1, fsst.out2[[i]]$test.logical)
  }
})
