context("Tests for error or warning messages")
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
## Full information approach
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
  # Variance
  var <- diag(length(unique(data[,"Y"])))
  return(list(beta = beta,
              var = var))
}

# ---------------- #
# Data preparation and declaration
# ---------------- #
# Declare parameters
N <- nrow(sampledata)
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1
pi <- 1 - mean(sampledata[,"D"])

# Compute matrices required
yp <- seq(0, 1, 1/J)
A_tgt <- matrix(c(yp, yp), nrow = 1)

# Define the observed matrix for each appraoch
A_obs_full <- cbind(matrix(rep(0, J1 * J1), nrow = J1), diag(1, J1))

# ---------------- #
# Shape constraints
# ---------------- #
A_shp_full <- matrix(rep(1, ncol(A_obs_full)), nrow = 1)
beta_shp <- c(1)

# ---------------- #
# Define arguments and produce output
# ---------------- #
# Parameters to test
beta.tgt <- .365
reps <- 100

# Define the lpmodels
lpmodel.full <- lpmodel(A.obs    = A_obs_full,
                        A.tgt    = A_tgt,
                        A.shp    = A_shp_full,
                        beta.obs = func_full_info,
                        beta.shp = beta_shp)

# Define the function that assigns the arguments
get.args <- function() {
  gen.args <- list(data = sampledata,
                   lpmodel = lpmodel.full,
                   beta.tgt = .3)
  dkqs.args <- gen.args
  dkqs.args$tau <- .2
  subs.args <- gen.args
  subs.args$phi <- 2/3
  fsst.args <- gen.args
  fsst.args$lambda <- .2
  minc.args <- gen.args
  minc.args$norm <- 2
  minc.args$beta.tgt <- NULL
  estb.args <- minc.args
  estb.args$kappa <- 1e-5

  return(list(dkqs = dkqs.args,
              subsample = subs.args,
              fsst = fsst.args,
              estbounds = estb.args,
              mincriterion = minc.args))
}

# ---------------- #
# 1. Tuning parameter
# ---------------- #
msg1a <- "The class of the variable '%s' has to be numeric."
msg1b <- "The object '%s' has to be a nonnegative number."
msg1.dkqs <- sprintf(msg1b, "tau")
msg1.subs <- sprintf(msg1a, "phi")
msg1.fsst.lam <- sprintf(msg1a, "lambda")
msg1.fsst.rho <- "The object 'rho' has to be a strictly positive number."
msg1.estb <- sprintf(msg1b, "kappa")
msg1c <- "The variable 'phi' has to be inside the interval (0, 1)."
msg1d <- "The variable 'lambda' has to be inside the interval [0, 1]"

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args1 <- args$fsst
fsst.args2 <- fsst.args1
estb.args <- args$estbounds

# Define function to test the tuning parameter
test.param <- function(param, string) {
  test_that(sprintf("Tuning parameter: %s", string), {
    dkqs.args$tau <- param
    subs.args$phi <- param
    fsst.args1$lambda <- param
    fsst.args2$rho <- param
    estb.args$kappa <- param
    expect_error(do.call(dkqs, dkqs.args), regexp = msg1.dkqs)
    expect_error(do.call(subsample, subs.args), regexp = msg1.subs)
    expect_error(do.call(fsst, fsst.args1), regexp = msg1.fsst.lam)
    expect_error(do.call(fsst, fsst.args2), regexp = msg1.fsst.rho)
    expect_error(do.call(estbounds, estb.args), regexp = msg1.estb)
  })
}

# Run tests
test.param("A", "string")
test.param("1", "string (numeric as string)")
df <- data.frame(matrix(1:6, nrow = 3, ncol = 2))
test.param(df, "data frame")

# NA as the tuning parameter - Should be no problem for FSST lambda
test_that("Tuning parameter as NA", {
  dkqs.args$tau <- NA
  subs.args$phi <- NA
  fsst.args1$lambda <- NA
  fsst.args2$rho <- NA
  estb.args$kappa <- NA
  expect_error(do.call(dkqs, dkqs.args), regexp = msg1.dkqs)
  expect_error(do.call(subsample, subs.args), regexp = msg1.subs)
  expect_error(do.call(fsst, fsst.args1), regexp = NA)
  expect_error(do.call(fsst, fsst.args2), regexp = msg1.fsst.rho)
  expect_error(do.call(estbounds, estb.args), regexp = msg1.estb)
})

# Define function to test the range of the parameter
test.param.range <- function(param) {
  test_that("Wrong range for tuning parameter", {
    dkqs.args$tau <- param
    subs.args$phi <- param
    fsst.args1$lambda <- param
    fsst.args2$rho <- param
    estb.args$kappa <- param
    expect_error(do.call(dkqs, dkqs.args), regexp = msg1.dkqs)
    expect_error(do.call(subsample, subs.args), regexp = msg1c, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args1), regexp = msg1d, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args2), regexp = msg1.fsst.rho)
    expect_error(do.call(estbounds, estb.args), regexp = msg1.estb)
  })
}
test.param.range(-1)

# ---------------- #
# 2. Solver not supported
# ---------------- #
msg2.0 <- paste0("This function with a %s-norm in the estimation is not ",
                 "compatible with '%s'. Please install one of the ",
                 "following packages: 'gurobi' (version 8.1-1 or later); ",
                 "'cplexAPI' (version 1.3.3 or later); 'Rcplex' ",
                 "(version 0.3-3 or later); 'limSolve' (version 1.5.6 or ",
                 "later)")
msg2.1 <- paste0(msg2.0, "; lpSolveAPI (version 5.5.2.0 or later).")
msg2.2 <- paste0(msg2.0, ".")
solver <- "lpinfersolver"
msg2.n1 <- sprintf(msg2.1, 1, solver)
msg2.n2 <- sprintf(msg2.2, 2, solver)
msg2.qc <- paste0("This function with a 2-norm in the estimation procedure ",
                  "and a quadratically constrained quadratic program is only ",
                  "compatible with 'gurobi'. Please install 'gurobi' ",
                  "(version 8.1-1 or later).")

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst
estb.args <- args$estbounds
minc.args <- args$mincriterion

# Run tests
test_that("Solver name not supported", {
  dkqs.args$solver <- solver
  subs.args$solver <- solver
  fsst.args$solver <- solver
  estb.args$solver <- solver
  minc.args$solver <- solver
  expect_error(do.call(dkqs, dkqs.args), regexp = msg2.n2, fixed = TRUE)
  expect_error(do.call(subsample, subs.args), regexp = msg2.n2, fixed = TRUE)
  expect_error(do.call(fsst, fsst.args), regexp = msg2.n2, fixed = TRUE)
  expect_error(do.call(estbounds, estb.args), regexp = msg2.qc, fixed = TRUE)
  expect_error(do.call(mincriterion, minc.args), regexp = msg2.n2, fixed = TRUE)
  estb.args$norm <- 1
  minc.args$norm <- 1
  expect_error(do.call(estbounds, estb.args), regexp = msg2.n1, fixed = TRUE)
  expect_error(do.call(mincriterion, minc.args), regexp = msg2.n1, fixed = TRUE)
})

# ---------------- #
# 3. Incorrect number of bootstraps
# ---------------- #
msg3 <- "The object 'R' has to be a positive integer."

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst

# Define function to test the number of bootstraps
test.r <- function(R, string) {
  test_that(sprintf("Bootstrap: R %s", string), {
    dkqs.args$R <- R
    subs.args$R <- R
    fsst.args$R <- R
    expect_error(do.call(dkqs, dkqs.args), regexp = msg3)
    expect_error(do.call(subsample, subs.args), regexp = msg3)
    expect_error(do.call(fsst, fsst.args), regexp = msg3)
  })
}

# Run tests
test.r(-1, "is negative")
test.r(0, "equals to zero")
test.r("A", "is a string")

# ---------------- #
# 4. Incorrect Rmulti parameter
# ---------------- #
msg4a <- "The class of the variable 'Rmulti' has to be numeric."
msg4b <- "The variable 'Rmulti' has to be inside the interval [1, Inf)."

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst

# Define function to test the 'Rmulti' parameter
test.rm <- function(Rmulti, string, msg) {
  test_that(sprintf("Rmulti %s", string, msg), {
    dkqs.args$Rmulti <- Rmulti
    subs.args$Rmulti <- Rmulti
    fsst.args$Rmulti <- Rmulti
    expect_error(do.call(dkqs, dkqs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
  })
}

# Run tests
test.rm("A", "is a string", msg4a)
test.rm(0, "is less than one", msg4b)

# ---------------- #
# 5. Incorrect norm
# ---------------- #
msg5 <- paste0("Only 1-norm and 2-norm are supported in this function. ",
               "For 1-norm, please use one of the followings: 1, '1', ",
               "'l1', 'one', 'o' or 'taxicab'. For 2-norm, please use ",
               "one of the followings: 2, '2', 'l2', 'two', 't', 'e', ",
               "'euclidean'.")

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst

# Run tests
test_that("Incorrect norm", {
  norm <- 3
  subs.args$norm <- norm
  estb.args$norm <- norm
  minc.args$norm <- norm
  expect_error(do.call(subsample, subs.args), regexp = msg5, fixed = TRUE)
  expect_error(do.call(estbounds, estb.args), regexp = msg5, fixed = TRUE)
  expect_error(do.call(mincriterion, minc.args), regexp = msg5, fixed = TRUE)
})

# ---------------- #
# 6. Incorrect 'beta.tgt'
# ---------------- #
msg6a <- "The class of the variable 'beta.tgt' has to be numeric."
msg6b <- 'argument "beta.tgt" is missing, with no default'
msg6c <- paste0("Computation is skipped because the parameter 'beta.tgt' ",
                "is outside the logical bound.")

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst

# Define function to test 'beta.tgt'
test.btgt.e <- function(beta.tgt, string, msg) {
  test_that(paste0("'beta.tgt': %s", string), {
    dkqs.args$beta.tgt <- beta.tgt
    subs.args$beta.tgt <- beta.tgt
    fsst.args$beta.tgt <- beta.tgt
    expect_error(do.call(dkqs, dkqs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
  })
}

# Run tests
test.btgt.e("1", "non-numeric", msg6a)
test.btgt.e("a", "non-numeric", msg6a)
test.btgt.e(NULL, "missing", msg6b)

# Define function to test the bounds of 'beta.tgt'
test.btgt.w <- function(beta.tgt, msg) {
  test_that("'beta.tgt': outside logical bounds", {
    dkqs.args$beta.tgt <- beta.tgt
    subs.args$beta.tgt <- beta.tgt
    fsst.args$beta.tgt <- beta.tgt
    expect_warning(do.call(dkqs, dkqs.args), regexp = msg, fixed = TRUE)
    expect_warning(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    expect_warning(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
  })
}

# Run tests
test.btgt.w(1e20, msg6c)
test.btgt.w(-1e20, msg6c)

# ---------------- #
# 7. Incorrect 'lpmodel' - Deterministic
# ---------------- #
msg7a <- paste0("The objects '%s' and '%s' in 'lpmodel' need to have ",
               "the same number of rows.")
msg7a.obs <- sprintf(msg7a, "A.obs", "beta.obs")
msg7a.shp <- sprintf(msg7a, "A.shp", "beta.shp")

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst
estb.args <- args$estbounds
minc.args <- args$mincriterion

# Define function to test the 'lpmodel'
test.lpmodel <- function(lpmodel, string, msg, dkqs.msg) {
  test_that(sprintf("lpmodel: %s", string, msg), {
    dkqs.args$lpmodel <- lpmodel
    subs.args$lpmodel <- lpmodel
    fsst.args$lpmodel <- lpmodel
    estb.args$lpmodel <- lpmodel
    minc.args$lpmodel <- lpmodel
    expect_error(do.call(dkqs, dkqs.args), regexp = dkqs.msg, fixed = TRUE)
    expect_error(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(estbounds, estb.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(mincriterion, minc.args), regexp = msg, fixed = TRUE)
  })
}

# A.obs and beta.obs have different number of rows
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$A.obs <- lpmodel.full.temp$A.obs[-1,]

# Run tests
test.lpmodel(lpmodel.full.temp,
             "A.obs and beta.obs have different number of rows",
             msg7a.obs,
             msg7a.obs)

# A.shp and beta.shp have different number of rows
# The shape constraints should not matter for 'dkqs'
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$A.shp <- rbind(lpmodel.full.temp$A.shp,
                                 lpmodel.full.temp$A.shp)
# Run tests
test.lpmodel(lpmodel.full.temp,
             "A.shp and beta.shp have different number of rows",
             msg7a.shp,
             NA)
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$beta.shp <- c(lpmodel.full.temp$beta.shp, 1)
## Run tests
test.lpmodel(lpmodel.full.temp,
             "A.shp and beta.shp have different number of rows",
             msg7a.shp,
             NA)

# Missing components
msg7b <- "The component '%s' is required in the 'lpmodel' object."
## A.obs
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$A.obs <- NULL
## Run tests
test.lpmodel(lpmodel.full.temp,
             "Missing A.obs",
             sprintf(msg7b, "A.obs"),
             sprintf(msg7b, "A.obs"))
## A.shp (Not required in 'dkqs')
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$A.shp <- NULL
## Run tests
test.lpmodel(lpmodel.full.temp,
             "Missing A.shp",
             sprintf(msg7b, "A.shp"),
             NA)
## A.tgt
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$A.tgt <- NULL
## Run tests
test.lpmodel(lpmodel.full.temp,
             "Missing A.tgt",
             sprintf(msg7b, "A.tgt"),
             sprintf(msg7b, "A.tgt"))
## beta.obs
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$beta.obs <- NULL
## Run tests
test.lpmodel(lpmodel.full.temp,
             "Missing beta.obs",
             sprintf(msg7b, "beta.obs"),
             sprintf(msg7b, "beta.obs"))
## beta.shp (Not required in 'dkqs')
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$beta.shp <- NULL
## Run tests
test.lpmodel(lpmodel.full.temp,
             "Missing beta.shp",
             sprintf(msg7b, "beta.shp"),
             NA)

# ---------------- #
# 8. Incorrect 'lpmodel' - Stochastic
# ---------------- #
# 'beta.obs' is a function but has only one output
# Only a problem for subsample (fsst will estimate it by bootstrap and it is
# not required for the rest)
msg8a <- paste0("The output of 'beta.obs' in 'lpmodel' needs to be a list ",
                "of two objects (one vector and one matrix).")
# Define function that does not return variance matrix
func_full_info_novar <- function(data) {
  beta <- func_full_info(data)$beta
  return(beta)
}

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst
estb.args <- args$estbounds
minc.args <- args$mincriterion

lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$beta.obs <- func_full_info_novar

# Define function to test the 'lpmodel' object
test.lpmodel.func <- function(lpmodel, string, msg1, msg2) {
  test_that(sprintf("lpmodel: %s", string), {
    dkqs.args$lpmodel <- lpmodel
    subs.args$lpmodel <- lpmodel
    fsst.args$lpmodel <- lpmodel
    estb.args$lpmodel <- lpmodel
    minc.args$lpmodel <- lpmodel
    expect_error(do.call(dkqs, dkqs.args), regexp = msg1, fixed = TRUE)
    expect_error(do.call(subsample, subs.args), regexp = msg2, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args), regexp = msg1, fixed = TRUE)
    expect_error(do.call(estbounds, estb.args), regexp = msg1, fixed = TRUE)
    expect_error(do.call(mincriterion, minc.args), regexp = msg1, fixed = TRUE)
  })
}

# Run tests
test.lpmodel.func(lpmodel.full.temp, "'beta.obs' without variance", NA, msg8a)

# Define function that returns problematic 'beta.obs' with nonzero probability
func_full_info_prob <- function(data) {
  ret <- func_full_info(data)
  beta <- ret$beta
  ru <- runif(1)
  # Assign one element of 'beta.obs' as NA with some probability
  if (ru < .05) {
    beta[1] <- NA
  }
  var <- ret$var
  return(list(beta = beta,
              var = var))
}
lpmodel.full.temp <- lpmodel.full
lpmodel.full.temp$beta.obs <- func_full_info_prob
dkqs.args$lpmodel <- lpmodel.full.temp
subs.args$lpmodel <- lpmodel.full.temp
fsst.args$lpmodel <- lpmodel.full.temp

# Obtain the outputs
out8 <- list()
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(1)
out8[[1]] <- do.call(dkqs, dkqs.args)
set.seed(1)
out8[[2]] <- do.call(subsample, subs.args)
set.seed(1)
out8[[3]] <- do.call(fsst, fsst.args)

# Expected number of errors and error messages in df.error
error.msg8 <- list()
error.msg8[[1]] <- list("missing value where TRUE/FALSE needed", 11)
error.msg8[[2]] <- list("missing value where TRUE/FALSE needed", 4)
error.msg8[[3]] <- list("model$obj contains NA", 8)

# There should be no error message printed
test.lpmodel.func2 <- function(lpmodel, string, msg) {
  test_that(sprintf("lpmodel: %s", string), {
    dkqs.args$lpmodel <- lpmodel
    subs.args$lpmodel <- lpmodel
    fsst.args$lpmodel <- lpmodel
    set.seed(1)
    expect_error(do.call(dkqs, dkqs.args), regexp = msg, fixed = TRUE)
    set.seed(1)
    expect_error(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    set.seed(1)
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
  })
}

# Run tests
test.lpmodel.func2(lpmodel.full.temp,
                   "'NA' in 'beta.obs' in some bootstraps",
                   NA)

# Check 'df.error'
test_that("Error messages in df.error", {
  for (i in 1:3) {
    expect_equal(sum(out8[[i]]$df.error[, ncol(out8[[i]]$df.error)] ==
                       error.msg8[[i]][[1]]),
                 error.msg8[[i]][[2]])
  }
})

# Test number of failed bootstraps (i.e. number of rows in 'df.error')
test_that("Number of rows in df.error", {
  for (i in 1:3) {
    expect_equal(nrow(out8[[i]]$df.error), error.msg8[[i]][[2]])
  }
})

# Test number of successful bootstraps
test_that("Number of successful bootstraps", {
  for (i in 1:3) {
    expect_equal(out8[[i]]$R.succ, 100)
  }
})

# ---------------- #
# 9. Missing or incorrect 'lpmodel'
# ---------------- #
msg9a <- 'argument "lpmodel" is missing, with no default'
msg9b <- "The object 'lpmodel' has to be an object in the 'lpmodel' class."

# Get arguments
args <- get.args()
dkqs.args <- args$dkqs
subs.args <- args$subsample
fsst.args <- args$fsst
estb.args <- args$estbounds
minc.args <- args$mincriterion

# Define function to test the 'lpmodel' object
test.lpmodel.wrong <- function(lpm, string, msg) {
  test_that(sprintf("%s 'lpmodel'", string), {
    dkqs.args$lpmodel <- lpm
    subs.args$lpmodel <- lpm
    fsst.args$lpmodel <- lpm
    estb.args$lpmodel <- lpm
    minc.args$lpmodel <- lpm
    expect_error(do.call(dkqs, dkqs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(subsample, subs.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(estbounds, estb.args), regexp = msg, fixed = TRUE)
    expect_error(do.call(mincriterion, minc.args), regexp = msg, fixed = TRUE)
  })
}

# Run tests
test.lpmodel.wrong(NULL, "Missing", msg9a)
test.lpmodel.wrong(list(1, 2, 3), "Incorrect", msg9b)

# ---------------- #
# 10. Other test-specific parameters
# ---------------- #
# 10(a): 'weight.matrix' in 'fsst'
msg10a <- "'weight.matrix' has to be one of 'avar', 'diag' and 'identity'."
## Assign arguments
args <- get.args()
fsst.args <- args$fsst
fsst.args$weight.matrix <- "lpinfer"
## Run test
test_that("Incorrect 'weight.matrix' in 'fsst'", {
  expect_error(do.call(fsst, fsst.args), regexp = msg10a, fixed = TRUE)
})

# 10(b): 'n' in 'fsst'
msg10b <- "'n' has to be a positive integer."
## Assign arguments
args <- get.args()
fsst.args <- args$fsst
## Define function to test the 'n' object
test.fsst.n <- function(n, string, msg) {
  test_that(sprintf("n %s", string), {
    fsst.args$n <- n
    expect_error(do.call(fsst, fsst.args), regexp = msg, fixed = TRUE)
  })
}
## Run tests
test.fsst.n(0, "equals to zero", msg10b)
test.fsst.n(-10, "is negative", msg10b)
test.fsst.n("A", "is a string", msg10b)

# 10(c): 'estimate' in 'estbounds'
msg.bool <- "The object '%s' has to be a Boolean expression."
msg10c <- sprintf(msg.bool, "estimate")
## Assign arguments
estb.args <- args$estbounds
estb.args$estimate <- "lpinfer"
## Run test
test_that("Nonboolean 'estimate' in 'estbounds'", {
  expect_error(do.call(estbounds, estb.args), regexp = msg10c, fixed = TRUE)
})

# 10(d): 'replace' in 'subsample'
msg10d <- sprintf(msg.bool, "replace")
## Assign arguments
subs.args <- args$subsample
subs.args$replace <- "lpinfer"
## Run test
test_that("Nonboolean 'estimate' in 'estbounds'", {
  expect_error(do.call(subsample, subs.args), regexp = msg10d, fixed = TRUE)
})
