context("Tests for fsst")
rm(list = ls())

# ---------------- #
# Load packages
# ---------------- #
library(lpinfer)
library(future)
library(future.apply)
library(dplyr)

# =========================================================================== #
# Case 1: d >= p
# =========================================================================== #
# ---------------- #
# Define functions to match the moments
# ---------------- #
## Full information approach
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

## Two moments approach
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
A_obs_full <- cbind(matrix(rep(0, J1*J1), nrow = J1), diag(1, J1))
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
lam1 <- c(.1, NA)
lam2 <- c(.1, .5, NA)
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
             lambda = lam1,
             rho = rho,
             n = NULL,
             R = reps,
             weight.matrix = "diag",
             solver = "gurobi",
             progress = TRUE)

# ---------------- #
# Output 1: beta.obs is a function
# d >= p
# ---------------- #
# Define the parameters
i.cores <- list(1, 8)
j.lpmodel <- list(lpmodel.full, lpmodel.twom)
k.lambdas <- list(lam1, lam2)

# Generate output
fsst.out <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  fsst.out[[i]] <- list()
  for (j in seq_along(j.lpmodel)) {
    farg$lpmodel <- j.lpmodel[[j]]
    fsst.out[[i]][[j]] <- list()
    for (k in seq_along(k.lambdas)) {
      RNGkind(kind = "L'Ecuyer-CMRG")
      set.seed(1)
      farg$lambda <- k.lambdas[[k]]
      fsst.out[[i]][[j]][[k]] <- do.call(fsst, farg)
    }
  }
}

# ---------------- #
# Output 2: beta.obs is a list that contains the sample and bootstrap estimates
# d >= p
# ---------------- #
# Function to draw bootstrap data
draw.bs.data <- function(x, f, data) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  bobs.bs <- f(data.bs)$beta
  return(bobs.bs)
}

# Draw bootstrap data for the full information and two moments method
set.seed(1)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- .Random.seed
bobs.bs.full.list <- future.apply::future_lapply(1:reps,
                                                 FUN = draw.bs.data,
                                                 future.seed = seed,
                                                 f = func_full_info,
                                                 data = sampledata)
bobs.bs.twom.list <- future.apply::future_lapply(1:reps,
                                                 FUN = draw.bs.data,
                                                 future.seed = seed,
                                                 f = func_two_moment,
                                                 data = sampledata)

bobs.full.list <- c(list(func_full_info(sampledata)$beta), bobs.bs.full.list)
bobs.twom.list <- c(list(func_two_moment(sampledata)$beta), bobs.bs.twom.list)
bobs.full.list[[1]] <- func_full_info(sampledata)
bobs.twom.list[[1]] <- func_two_moment(sampledata)

# Redefine the 'lpmodel' object with 'beta.obs' being a list
lpmodel.full.list <- lpmodel.full
lpmodel.twom.list <- lpmodel.twom
lpmodel.full.list$beta.obs <- bobs.full.list
lpmodel.twom.list$beta.obs <- bobs.twom.list

# Define the new lpmodel object and the arguments to be passed to the function
j.lpmodel2 <- list(lpmodel.full.list, lpmodel.twom.list)
farg2 <- list(beta.tgt = beta.tgt,
              lambda = lam1,
              rho = rho,
              n = nrow(sampledata),
              R = reps,
              weight.matrix = "diag",
              solver = "gurobi",
              progress = TRUE)

# Compute the fsst output again
fsst.out2 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  fsst.out2[[i]] <- list()
  for (j in seq_along(j.lpmodel2)) {
    farg2$lpmodel <- j.lpmodel2[[j]]
    fsst.out2[[i]][[j]] <- list()
    for (k in seq_along(k.lambdas)) {
      RNGkind(kind = "L'Ecuyer-CMRG")
      set.seed(1)
      farg2$lambda <- k.lambdas[[k]]
      fsst.out2[[i]][[j]][[k]] <- do.call(fsst, farg2)
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
## Function to draw bootstrap data and the corresponding estimate of 'beta.obs'
fsst.bs.fn <- function(x, data, lpmodel) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  beta <- lpmodel$beta.obs(data.bs)$beta
  return(beta)
}
## Estimator of asymptotic variance of beta.obs
set.seed(1)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- .Random.seed

beta.bs <- list()
for (j in seq_along(j.lpmodel)) {
  beta.bs[[j]] <- future.apply::future_lapply(1:reps,
                                              FUN = fsst.bs.fn,
                                              future.seed = seed,
                                              data = sampledata,
                                              lpmodel = j.lpmodel[[j]])
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
  colnames(A) <- 1:ncol(A)
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
    A1 <- cbind(args$A,
                matrix(0L,
                       nrow = nrow(args$A),
                       ncol = ncol(A)))
    A2 <- cbind(-diag(ns),
                matrix(0L, nrow = ns, ncol = p*2),
                A)
    colnames(A1) <- 1:ncol(A1)
    colnames(A2) <- 1:ncol(A2)
    args1 <- list(Q = NULL,
                  obj = c(args$obj, rep(0, nx)),
                  objcon = 0,
                  A = as.matrix(rbind(A1, A2)),
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

# 6. Obtain the data-driven lambda
lambda.dd1 <- list()
lambda.dd.ts <- list()
alpha.n <- 1/sqrt(log(log(N)))
for (j in seq_along(j.lpmodel)) {
  lambda.dd.ts[[j]] <- list()
  for (i in 1:reps) {
    # Solve the bootstrap cone problem
    cone.dd.args <- fsst.56.args(j.lpmodel[[j]],
                                 p[[j]],
                                 d[[j]],
                                 beta.star.bs[[j]][[i]] - beta.star[[j]],
                                 omega[[j]])
    lambda.dd.ts[[j]][[i]] <- sqrt(N) * do.call(gurobi.qlp,
                                                cone.dd.args)$objval
  }
  nq <- ceiling(reps * (1 - alpha.n))
  q <- sort(unlist(lambda.dd.ts[[j]]))[nq]
  lambda.dd1[[j]] <- min(1/q, 1)
}

# 7. Restricted estimator
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

# 8. Compute bootstrap components
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
        ## Replace NA with data-driven lambda
        k.lambdas[[k]][is.na(k.lambdas[[k]])] <- lambda.dd1[[j]]
        ## Compute the bootstrap estimates
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

# 9. Compute pvalues
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

# 10. Critical values
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
# A list of unit tests for d >= p
# i: cores, j: lpmodel approach, k: lambdas
# ---------------- #
tests.fsst.dgeqp <- function(fsst.out, test.name) {
  # Assign the name
  test.name1 <- sprintf("'beta.obs' as %s:", test.name)

  # 1. Full information approach p-values
  test_that(sprintf("%s 'd >= p': Full information approach", test.name1), {
    for (i in seq_along(i.cores)) {
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
  test_that(sprintf("%s 'd >= p': Full information approach", test.name1), {
    for (i in seq_along(i.cores)) {
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
  test_that(sprintf("%s 'd >= p': Full information CV table", test.name1), {
    for (i in seq_along(i.cores)) {
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
  test_that(sprintf("%s 'd >= p': Two moments CV table", test.name1), {
    for (i in seq_along(i.cores)) {
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

  # 5. Range test statistics
  test_that(sprintf("%s 'd >= p': Range test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(range.n[[j]], fsst.out[[i]][[j]][[k]]$range)
        }
      }
    }
  })

  # 6. Cone test statistics
  test_that(sprintf("%s 'd >= p': Cone test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(cone.n[[j]], fsst.out[[i]][[j]][[k]]$cone)
        }
      }
    }
  })

  # 7. Test statistics
  test_that(sprintf("%s 'd >= p': Test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(ts[[j]], fsst.out[[i]][[j]][[k]]$test)
        }
      }
    }
  })

  # 8. Solver name
  test_that(sprintf("%s 'd >= p': Solver name", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal("gurobi", fsst.out[[i]][[j]][[k]]$solver.name)
        }
      }
    }
  })

  # 9. Rho parameter
  test_that(sprintf("%s 'd >= p': Rho parameter", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(rho, fsst.out[[i]][[j]][[k]]$rho)
        }
      }
    }
  })

  # 10. Regularization parameter
  test_that(sprintf("%s 'd >= p': Regularization parameter", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(rhobar[[j]], fsst.out[[i]][[j]][[k]]$rhobar.i)
        }
      }
    }
  })

  # 11. Data-driven lambda
  test_that(sprintf("%s 'd >= p': Data-driven lambda", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(lambda.dd1[[j]], fsst.out[[i]][[j]][[k]]$lambda.data)
        }
      }
    }
  })

  # 12. Method of obtaining the beta.var matrix
  test_that(sprintf("%s 'd >= p': Method of obtaining beta.var", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(test.name, fsst.out[[i]][[j]][[k]]$var.method)
        }
      }
    }
  })

  # 13. Omega.i matrix
  test_that(sprintf("%s 'd >= p': Omega.i matrix", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(omega[[j]], fsst.out[[i]][[j]][[k]]$omega.i)
        }
      }
    }
  })

  # 14. Test logical
  test_that(sprintf("%s 'd >= p': Test logical", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in 1:2) {
        for (k in 1:2) {
          expect_equal(1, fsst.out[[i]][[j]][[k]]$test.logical)
        }
      }
    }
  })
}

# ---------------- #
# Run the tests for d >= p
# ---------------- #
# beta.obs is a function
tests.fsst.dgeqp(fsst.out, "function")

# beta.obs is a list of bootstrap estimates
tests.fsst.dgeqp(fsst.out2, "list")


# =========================================================================== #
# Case 2: d < p
# =========================================================================== #
# ---------------- #
# Define the 'lpmodel' object
# ---------------- #
# Function for beta.obs
betafunc <- function(data) {
  beta <- matrix(0L, nrow = 2)
  beta[1, 1] <- mean(data$Y > .5)
  beta[2, 1] <- 1 - beta[1, 1]
  return(list(beta = beta,
              var = diag(2)))
}

# Other components in the 'lpmodel' object
Aobs <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE)
Ashp <- matrix(c(.2, .2, .6), nrow = 1, byrow = TRUE)
Atgt <- matrix(c(.3, .4, 1), nrow = 1, byrow = TRUE)
bshp <- c(.8)

# Formulate the 'lpmodel' object
lpm2 <- lpmodel(A.obs = Aobs,
                A.shp = Ashp,
                A.tgt = Atgt,
                beta.obs = betafunc,
                beta.shp = bshp)

# Other paramters - Use the same lambda as in the earlier test
btgt <- 1.5
rho2 <- 1e-3
reps <- 100
n <- nrow(sampledata)

# Define arguments for FSST
farg3 <- list(data = sampledata,
              lpmodel = lpm2,
              beta.tgt = btgt,
              rho = rho2,
              n = NULL,
              R = reps,
              weight.matrix = "avar",
              solver = "gurobi",
              lambda = lam1,
              progress = TRUE)

# ---------------- #
# Output 3: beta.obs is a function
# d < p
# ---------------- #
# Generate output
fsst.out3 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(1)
  fsst.out3[[i]] <- do.call(fsst, farg3)
}

# ---------------- #
# Output 4: beta.obs is a list that contains the sample and bootstrap estimates
# d < p
# ---------------- #
# Draw bootstrap data for the full information and two moments method
set.seed(1)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- .Random.seed
bobs.dlp.list <- future.apply::future_lapply(1:reps,
                                             FUN = draw.bs.data,
                                             future.seed = seed,
                                             f = betafunc,
                                             data = sampledata)

bobs.dlp.list <- c(list(betafunc(sampledata)$beta), bobs.dlp.list)
bobs.dlp.list[[1]] <- betafunc(sampledata)

# Redefine the 'lpmodel' object with 'beta.obs' being a list
lpm2.list <- lpm2
lpm2.list$beta.obs <- bobs.dlp.list

# Define the new lpmodel object and the arguments to be passed to the function
farg4 <- list(lpmodel = lpm2.list,
              beta.tgt = btgt,
              rho = rho2,
              n = nrow(sampledata),
              R = reps,
              weight.matrix = "avar",
              solver = "gurobi",
              lambda = lam1,
              progress = TRUE)

# Generate output
fsst.out4 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(1)
  fsst.out4[[i]] <- do.call(fsst, farg4)
}

# ---------------- #
# Construct the answer by the programs
# ---------------- #
# 1. Obtain the parameters and the matrices
A2 <- rbind(Aobs, Ashp, Atgt)
p2 <- nrow(A2)
d2 <- ncol(A2)
# sigma.bobs2 <- matrix(0L, nrow = p2 + 2, ncol = p2 + 2)
sigma.bobs2 <- lpm2$beta.obs(sampledata)$var

# 2. Compute beta.obs
## Sample beta.obs
bobs2 <- list()
bobs2 <- lpm2$beta.obs(sampledata)$beta
## Estimator of asymptotic variance of beta.obs
set.seed(1)
RNGkind(kind = "L'Ecuyer-CMRG")
seed <- .Random.seed
bobs.bs2 <- list()
bobs.bs2 <- future.apply::future_lapply(1:reps,
                                        FUN = fsst.bs.fn,
                                        future.seed = seed,
                                        data = sampledata,
                                        lpmodel = lpm2)

# 3. Solve problem (3) - with the sample estimates and the bootstrap estimates
fsst.3.args2 <- fsst.3.arg(lpm2, bobs2, btgt, sigma.bobs2, "avar")
Xi2 <- fsst.3.args2$Xi
fsst.3.return2 <- do.call(gurobi.qlp, fsst.3.args2$args)
Xi2 <- fsst.3.args2$Xi
fsst.3.return2 <- do.call(gurobi.qlp, fsst.3.args2$args)
x.star2 <- fsst.3.return2$x
if (d2 >= p2) {
  bobs.star2 <- bobs2
} else {
  bobs.star2 <- lpm2$A.obs %*% x.star2
}
beta.star2 <- c(bobs.star2, lpm2$beta.shp, btgt)
## Bootstrap components
x.star2.list <- list()
beta.obs.star.bs2 <- list()
beta.star.bs2 <- list()
for (i in 1:reps) {
  fsst.3.bs.args2 <- fsst.3.arg(lpm2, bobs.bs2[[i]], btgt,
                                sigma.bobs2, "avar")
  fsst.3.bs.return2 <- do.call(gurobi.qlp, fsst.3.bs.args2$args)
  x.star2.list[[i]] <- fsst.3.bs.return2$x
  if (d2 >= p2) {
    beta.obs.star.bs2[[i]] <- bobs.bs2[[i]]
  } else {
    beta.obs.star.bs2[[i]] <- lpm2$A.obs %*% x.star2.list[[i]]
  }
  beta.star.bs2[[i]] <- c(beta.obs.star.bs2[[i]], lpm2$beta.shp, btgt)
}
sigma2 <- matrix(0, nrow = p2, ncol = p2)
n.bobs2 <- length(bobs2)
sigma2[1:n.bobs2, 1:n.bobs2] <- sigma.bobs2

# 4. Standardization
## Compute studentization matrix
if (d2 >= p2) {
  student.matrix2 <- sigma2
} else {
  student.matrix2 <- matrix(0, nrow = p2, ncol = p2)
}
for (i in 1:reps) {
  beta.t <- beta.star.bs2[[i]] - beta.star2
  student.matrix2 <- student.matrix2 + beta.t %*% t(beta.t)
}
student.matrix2 <- n * student.matrix2 / reps
## Compute regularization parameter
rhobar2 <- base::norm(student.matrix2, type = "f") * rho2
## Compute regularization matrix
omega2 <- expm::sqrtm(student.matrix2 + rhobar2 * diag(p2))

# 5. Test statistics
## Range
range.n2 <- fsst.range.soln(bobs.star2, bobs2, Xi2, p2, d2)
## Cone
cone.args2 <- fsst.56.args(lpm2, p2, d2, beta.star2, omega2)
cone.return2 <- do.call(gurobi.qlp, cone.args2)
cone.n2 <- sqrt(n) * cone.return2$objval
## Test statistics
ts2 <- max(cone.n2, range.n2)

# 6. Obtain the data-driven lambda
alpha.n2 <- 1/sqrt(log(log(n)))
lambda.dd.ts2 <- NULL
for (i in 1:reps) {
  # Solve the bootstrap cone problem
  cone.dd.args <- fsst.56.args(lpm2, p2, d2,
                               beta.star.bs[[j]][[i]] - beta.star[[j]],
                               omega[[j]])
  lambda.dd.ts2 <- c(lambda.dd.ts2,
                     sqrt(n) * do.call(gurobi.qlp, cone.dd.args)$objval)
}
nq2 <- ceiling(reps * (1 - alpha.n2))
q2 <- sort(lambda.dd.ts2)[nq2]
lambda.dd2 <- min(1/q2, 1)

# 7. Restricted estimator
beta.r.args2 <- fsst.89.args(lpm2, p2, d2, beta.star2,
                             c(bobs2, lpm2$beta.shp, btgt), omega2, btgt)
beta.r2 <- do.call(gurobi.qlp, beta.r.args2)$x[1:p2]

# 8. Compute bootstrap components
range.n.bs2 <- list()
cone.n.bs2 <- list()
ts.bs2 <- list()
## Replace NA with data-driven lambda
lam2[is.na(lam2)] <- lambda.dd2
for (i in 1:reps) {
  ## Range
  b1 <- bobs.bs2[[i]] - bobs2
  b2 <- beta.obs.star.bs2[[i]] - bobs.star2
  range.n.bs2[[i]] <- fsst.range.soln(b1, b2, Xi2, p2, d2)

  ## Cone
  beta.res2 <- beta.star.bs2[[i]] - beta.star2 + lam2[1] * beta.r2
  cone.args2 <- fsst.56.args(lpm2, p2, d2, beta.res2, omega2)
  cone.n.bs2[[i]] <- sqrt(n) * do.call(gurobi.qlp, cone.args2)$objval
  ts.bs2[[i]] <- max(cone.n.bs2[[i]], range.n.bs2[[i]])
}

# 9. Compute p-value
pval2 <- mean(unlist(ts.bs2) > ts2)

# 10. Critical values
n99.2 <- ceiling(.99 * reps)
n95.2 <- ceiling(.95 * reps)
n90.2 <- ceiling(.90 * reps)
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
# A list of unit tests - For d < p only
# i: cores
# ---------------- #
tests.fsst.dlp <- function(fsst.out, test.name) {
  # Assign the name
  test.name1 <- sprintf("'beta.obs' as %s:", test.name)

  # 1. p-value
  test_that(sprintf("%s 'd < p': p-value", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(pval2, fsst.out[[i]]$pval$`p-value`[1])
    }
  })

  # 2. CV table
  test_that(sprintf("%s 'd < p': CV table", test.name1), {
    for (i in seq_along(i.cores)) {
      for (l in 1:12) {
        expect_lte(abs(cv2[l] - fsst.out[[i]]$cv.table[l, 3]), 1e-5)
      }
    }
  })

  # 3. Range test statistics
  test_that(sprintf("%s 'd < p': Range test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(range.n2, fsst.out[[i]]$range)
    }
  })

  # 4. Cone test statistics
  test_that(sprintf("%s 'd < p': Cone test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(cone.n2, fsst.out[[i]]$cone)
    }
  })

  # 5. Test statistics
  test_that(sprintf("%s 'd < p': Test statistics", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(ts2, fsst.out[[i]]$test)
    }
  })

  # 6. Solver name
  test_that(sprintf("%s 'd < p': Solver name", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal("gurobi", fsst.out[[i]]$solver.name)
    }
  })

  # 7. Rho parameter
  test_that(sprintf("%s 'd < p': Rho parameter", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(rho2, fsst.out[[i]]$rho)
    }
  })

  # 8. Regularization parameter
  test_that(sprintf("%s 'd < p': Regularization parameter", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(rhobar2, fsst.out[[i]]$rhobar.i)
    }
  })

  # 9. Data-driven lambda
  test_that(sprintf("%s 'd < p': Data-driven lambda", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(lambda.dd2, fsst.out[[i]]$lambda.data)
    }
  })

  # 10. Method of obtaining the beta.var matrix
  test_that(sprintf("%s 'd < p': Method of obtaining beta.var", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(test.name, fsst.out[[i]]$var.method)
    }
  })

  # 11. Omega.i
  test_that(sprintf("%s 'd < p': Omega.i matrix", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(omega2, fsst.out[[i]]$omega.i)
    }
  })

  # 12. Test logical
  test_that(sprintf("%s 'd < p': Test logical", test.name1), {
    for (i in seq_along(i.cores)) {
      expect_equal(1, fsst.out[[i]]$test.logical)
    }
  })
}


# ---------------- #
# Run the tests for d < p
# ---------------- #
# beta.obs is a function
tests.fsst.dlp(fsst.out3, "function")

# beta.obs is a list of bootstrap estimates
tests.fsst.dlp(fsst.out4, "list")
