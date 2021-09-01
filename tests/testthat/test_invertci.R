context("Tests for invertci")
rm(list = ls())

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)
library(future)
library(furrr)

# ---------------- #
# Define the lpmodel object
# ---------------- #
# Define the betaobs function
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
    beta <- c(beta, c(beta_i))
  }
  beta <- as.matrix(beta)
  # Variance
  var <- diag(length(unique(data[,"Y"])))
  return(list(beta = beta,
              var = var))
}

# Declare parameters
N <- dim(sampledata)[1]
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])
reps <- 10
# Compute matrices
yp <- seq(0, 1, 1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)
A_obs_full <- cbind(matrix(rep(0, J1*J1), nrow = J1), diag(1, J1))

# Shape constraints
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
beta_shp <- c(1)

# Define the target beta and the value of tau
tau <- sqrt(log(N)/N)
beta.tgt <- .365
phi <- 2/3
m <- floor(N^phi)
max.iter <- 50
alpha <- c(.05, .1)
tol <- 0.001
init.lb <- c(0, .4)
init.ub <- c(.6, 1)

# lpmodel object
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_full,
                        beta.obs = func_full_info,
                        beta.shp = beta_shp)

# Define arguments for the `subsample` function
farg <- list(data = sampledata,
             lpmodel = lpmodel.full,
             R = reps,
             norm = 2,
             phi = phi,
             replace = FALSE,
             progress = TRUE)

# ---------------- #
# Run invertci procedure
# ---------------- #
invertci.out <- list()
set.seed(1)
invertci.out[[1]] <- invertci(f = subsample,
                              farg = farg,
                              alpha = alpha,
                              init.lb = init.lb,
                              init.ub = init.ub,
                              tol = tol,
                              max.iter = max.iter,
                              pvals = NULL,
                              progress = FALSE)

# ---------------- #
# Run invertci procedure again with the given pvals
# ---------------- #
set.seed(1)
invertci.out[[2]] <- invertci(f = subsample,
                              farg = farg,
                              alpha = alpha,
                              init.lb = init.lb,
                              init.ub = init.ub,
                              tol = tol,
                              max.iter = max.iter,
                              pvals = invertci.out[[1]]$pvals,
                              progress = FALSE)

# ---------------- #
# Construct the result from scratch
# ---------------- #
bisection.ci.test <- function(init.b, category, alpha) {
  a <- init.b[1]
  b <- init.b[2]
  c <- (a + b)/2
  iter <- 1

  while ((b - a) > tol & iter <= max.iter) {
    # Run subsample again
    farg$beta.tgt <- c
    set.seed(1)
    temp <- do.call(subsample, farg)
    pval.temp <- temp$pval[1, 2]
    pval.location <- (pval.temp >= alpha)

    # Decide which is the next interval
    if (category == "lb") {
      if (isTRUE(pval.location)) {
        b <- c
      } else {
        a <- c
      }
    } else if (category == "ub") {
      if (isTRUE(pval.location)) {
        a <- c
      } else {
        b <- c
      }
    }

    # Update iteration number and the test point
    iter <- iter + 1
    c <- (a + b)/2
  }

  return(c)
}

# Get the bounds
lb <- list()
ub <- list()
for (i in seq_along(alpha)) {
  lb[[i]] <- bisection.ci.test(init.lb, "lb", alpha[i])
  ub[[i]] <- bisection.ci.test(init.ub, "ub", alpha[i])
}

# ---------------- #
# Unit tests
# ---------------- #
# 1. Lower bound
test_that("Lower bound", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(lb[[i]], invertci.out[[j]]$ci$lb[i])
    }
  }
})

# 2. Upper bound
test_that("Upper bound", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(ub[[i]], invertci.out[[j]]$ci$ub[i])
    }
  }
})

# 3. Tolerance level
test_that("Tolerance level", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(tol, invertci.out[[j]]$tol)
    }
  }
})

# 4. alpha
test_that("alpha", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(alpha, invertci.out[[j]]$alpha)
    }
  }
})

# 5. max.iter
test_that("max.iter", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(max.iter, invertci.out[[j]]$max.iter)
    }
  }
})

# 6. Reason for termination
test_that("Reason for termination", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      for (k in c("lb", "ub")) {
        expect_equal("Length of interval is below tolerance level",
                     invertci.out[[j]]$termination[[i]][[1]][[k]])
      }
    }
  }
})

# 7. Name of parameter
test_that("Name of parameter", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal("phi", invertci.out[[j]]$para.name)
    }
  }
})

# 8. Values of tuning parameters
test_that("Values of tuning parameters", {
  for (i in seq_along(alpha)) {
    for (j in 1:2) {
      expect_equal(phi, invertci.out[[j]]$para.vals[1, 1])
    }
  }
})

# 9. Same pvals data frame
test_that("Same pvals data frame", {
  expect_equal(invertci.out[[1]]$para.vals, invertci.out[[2]]$para.vals)
})
