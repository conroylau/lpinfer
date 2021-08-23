context("Tests for chorussell")

# ---------------- #
# Load relevant packages
# ---------------- #
library(lpinfer)
library(future)
library(furrr)

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
    beta <- c(beta, c(beta_i))
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
  beta <- matrix(c(0, 0), nrow = 2)
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
kappa <- 0
reps <- 100
tol <- 1e-4

# Define arguments for the `chorussell` function
farg <- list(data = sampledata,
             beta.tgt = beta.tgt,
             R = reps,
             kappa = kappa,
             estimate = TRUE,
             tol = tol,
             ci = FALSE,
             solver = "gurobi",
             progress = FALSE)

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


# ---------------- #
# Output 1: beta.obs is a function and output is p-value
# ---------------- #
# List of cores, lpmodel and norm objects to be used
i.cores <- list(1, 8)
j.lpmodel <- list(lpmodel.full, lpmodel.twom)
k.norm <- list(1, 2)

# Generate output
cr.out.pval <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  cr.out.pval[[i]] <- list()
  for (j in seq_along(j.lpmodel)) {
    farg$lpmodel <- j.lpmodel[[j]]
    cr.out.pval[[i]][[j]] <- list()
    for (k in seq_along(k.norm)) {
      set.seed(1)
      farg$norm <- k.norm[[k]]
      cr.out.pval[[i]][[j]][[k]] <- do.call(chorussell, farg)
    }
  }
}

# ---------------- #
# Output 2: beta.obs is a function and output is CI
# ---------------- #
farg$ci <- TRUE
# Generate output
cr.out.ci <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  cr.out.ci[[i]] <- list()
  for (j in seq_along(j.lpmodel)) {
    farg$lpmodel <- j.lpmodel[[j]]
    cr.out.ci[[i]][[j]] <- list()
    for (k in seq_along(k.norm)) {
      set.seed(1)
      farg$norm <- k.norm[[k]]
      cr.out.ci[[i]][[j]][[k]] <- do.call(chorussell, farg)
    }
  }
}

# ---------------- #
# Output 3: beta.obs is a list that contains the sample and bootstrap estimates
# and output is p-value
# ---------------- #
# Function to draw bootstrap data
draw.bs.data <- function(x, f, data) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  bobs.bs <- f(data.bs)$beta
  return(bobs.bs)
}

# Draw bootstrap data for the full information and two moments method
set.seed(1)
bobs.bs.full.list <- furrr::future_map(1:reps,
                                       .f = draw.bs.data,
                                       f = func_full_info,
                                       data = sampledata,
                                       .options =
                                         furrr::furrr_options(seed = TRUE))
set.seed(1)
bobs.bs.twom.list <- furrr::future_map(1:reps,
                                       .f = draw.bs.data,
                                       f = func_two_moment,
                                       data = sampledata,
                                       .options =
                                         furrr::furrr_options(seed = TRUE))

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
              R = reps,
              kappa = kappa,
              estimate = TRUE,
              n = nrow(sampledata),
              solver = "gurobi",
              tol = tol,
              ci = FALSE,
              progress = FALSE)

# Compute the chorussell output again
cr.out.pval2 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  cr.out.pval2[[i]] <- list()
  for (j in seq_along(j.lpmodel2)) {
    farg2$lpmodel <- j.lpmodel2[[j]]
    cr.out.pval2[[i]][[j]] <- list()
    for (k in seq_along(k.norm)) {
      set.seed(1)
      farg2$norm <- k.norm[[k]]
      cr.out.pval2[[i]][[j]][[k]] <- do.call(chorussell, farg2)
    }
  }
}

# ---------------- #
# Output 4: beta.obs is a list that contains the sample and bootstrap estimates
# and output is a list
# ---------------- #
farg2$ci <- TRUE
cr.out.ci2 <- list()
for (i in seq_along(i.cores)) {
  plan(multisession, workers = i.cores[[i]])
  cr.out.ci2[[i]] <- list()
  for (j in seq_along(j.lpmodel2)) {
    farg2$lpmodel <- j.lpmodel2[[j]]
    cr.out.ci2[[i]][[j]] <- list()
    for (k in seq_along(k.norm)) {
      set.seed(1)
      farg2$norm <- k.norm[[k]]
      cr.out.ci2[[i]][[j]][[k]] <- do.call(chorussell, farg2)
    }
  }
}

# ---------------- #
# Obtain the results without using the chorussell function
# ---------------- #
# 1. Obtain the sample and bootstrap bounds
## Function to get the bootstrap bounds
cr.bs.fn <- function(x, data, lpmodel, kappa, norm, estimate) {
  data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
  temp <- estbounds(data.bs, lpmodel, kappa, norm, estimate)
  return(list(ub = temp$ub,
              lb = temp$lb))
}

## Run
ub <- list()
lb <- list()
ub.bs <- list()
lb.bs <- list()
delta <- list()
for (j in seq_along(j.lpmodel)) {
  ub[[j]] <- list()
  lb[[j]] <- list()
  ub.bs[[j]] <- list()
  lb.bs[[j]] <- list()
  delta[[j]] <- list()
  for (k in seq_along(k.norm)) {
    # Sample bounds
    temp <- lpinfer::estbounds(sampledata, j.lpmodel[[j]], kappa = kappa,
                               norm = k.norm[[k]])
    ub[[j]][[k]] <- temp$ub
    lb[[j]][[k]] <- temp$lb
    delta[[j]][[k]] <- temp$ub - temp$lb
    # Bootstrap bounds
    set.seed(1)
    temp <- furrr::future_map(1:reps,
                             .f = cr.bs.fn,
                             lpmodel = j.lpmodel[[j]],
                             data = sampledata,
                             kappa = kappa,
                             norm = k.norm[[k]],
                             estimate = TRUE,
                             .options = furrr::furrr_options(seed = TRUE))
    ub.bs[[j]][[k]] <- unlist(sapply(temp, "[", "ub"), use.names = FALSE)
    lb.bs[[j]][[k]] <- unlist(sapply(temp, "[", "lb"), use.names = FALSE)
  }
}

# 2. Compute the list of candidate bounds
n <- nrow(sampledata)
ub.can1 <- list()
ub.can2 <- list()
lb.can1 <- list()
lb.can2 <- list()
for (j in seq_along(j.lpmodel)) {
  ub.can1[[j]] <- list()
  ub.can2[[j]] <- list()
  lb.can1[[j]] <- list()
  lb.can2[[j]] <- list()
  for (k in seq_along(k.norm)) {
    ub.can1[[j]][[k]] <- sqrt(n) * (ub.bs[[j]][[k]] - ub[[j]][[k]])
    ub.can2[[j]][[k]] <- sqrt(n) * (ub.bs[[j]][[k]] - ub[[j]][[k]] +
                                      delta[[j]][[k]])
    lb.can1[[j]][[k]] <- sqrt(n) * (lb.bs[[j]][[k]] - lb[[j]][[k]])
    lb.can2[[j]][[k]] <- sqrt(n) * (lb.bs[[j]][[k]] - lb[[j]][[k]] -
                                      delta[[j]][[k]])
    }
}

# 3. Create solver for the optimization problem
cr.optim <- function(lb.can1, lb.can2, ub.can1, ub.can2, lb, ub, alpha) {
  ## Combine the list of candidates
  lb.can <- c(lb.can1, lb.can2)
  ub.can <- -c(ub.can1, ub.can2)

  ## Data frame to store the results
  df <- data.frame(matrix(vector(), nrow = 0, ncol = 3))
  colnames(df) <- c("lb", "ub", "sum")

  for (i in seq_along(lb.can)) {
    for (j in seq_along(ub.can)) {
      ind1 <- (mean((lb.can1 <= lb.can[i]) * (-ub.can[j] <= ub.can2))
               >= (1 - alpha))
      ind2 <- (mean((lb.can2 <= lb.can[i]) * (-ub.can[j] <= ub.can1))
               >= (1 - alpha))
      if (isTRUE(ind1) & isTRUE(ind2)) {
        k <- nrow(df) + 1
        df[k, 1] <- lb.can[i]
        df[k, 2] <- ub.can[j]
        df[k, 3] <- lb.can[i] + ub.can[j]
      }
    }
  }
  bd <- df[df[,3] == min(df[,3]), 1:2]
  return(list(lb = lb - bd[1,1]/sqrt(n),
              ub = ub + bd[1,2]/sqrt(n)))
}

# 4. Construct confidence interval
c.lb <- list()
c.ub <- list()
for (j in seq_along(j.lpmodel)) {
  c.lb[[j]] <- list()
  c.ub[[j]] <- list()
  for (k in seq_along(k.norm)) {
    temp <- cr.optim(lb.can1[[j]][[k]], lb.can2[[j]][[k]], ub.can1[[j]][[k]],
                     ub.can2[[j]][[k]], lb[[j]][[k]], ub[[j]][[k]], alpha = .05)
    c.lb[[j]][[k]] <- temp$lb
    c.ub[[j]][[k]] <- temp$ub
  }
}

# 5. Construct p-value by bisection method
pval <- list()
for (j in seq_along(j.lpmodel)) {
  pval[[j]] <- list()
  for (k in seq_along(k.norm)) {
    # Check end points
    a <- 0
    a.lp <- cr.optim(lb.can1[[j]][[k]], lb.can2[[j]][[k]], ub.can1[[j]][[k]],
                     ub.can2[[j]][[k]], lb[[j]][[k]], ub[[j]][[k]], alpha = a)
    a.inout <- (beta.tgt <= a.lp$ub & beta.tgt >= a.lp$lb)
    if (isFALSE(a.inout)) {
      pval[[j]][[k]] <- a
      next()
    }
    b <- 1
    b.lp <- cr.optim(lb.can1[[j]][[k]], lb.can2[[j]][[k]], ub.can1[[j]][[k]],
                     ub.can2[[j]][[k]], lb[[j]][[k]], ub[[j]][[k]], alpha = b)
    b.inout <- (beta.tgt <= b.lp$ub & beta.tgt >= b.lp$lb)
    if (isTRUE(b.inout)) {
      pval[[j]][[k]] <- b
      next()
    }
    while (abs(b - a) > tol) {
      c <- (a + b)/2
      c.lp <- cr.optim(lb.can1[[j]][[k]], lb.can2[[j]][[k]], ub.can1[[j]][[k]],
                       ub.can2[[j]][[k]], lb[[j]][[k]], ub[[j]][[k]], alpha = c)
      c.inout <- (beta.tgt <= c.lp$ub & beta.tgt >= c.lp$lb)
      if (isFALSE(c.inout)) {
        b <- c
      } else {
        a <- c
      }
    }
    pval[[j]][[k]] <- c
  }
}

# ---------------- #
# A list of unit tests for p-value
# ---------------- #
test.cr.pval <- function(cr.out, test.name, test.type) {
  # Assign the name
  test.name1 <- sprintf("'beta.obs' as %s and test for %s:",
                        test.name,
                        test.type)

  # 1. Full information approach p-values
  test_that(sprintf("%s Full information approach", test.name1), {
    for (i in seq_along(i.cores)) {
      j <- 1
      for (k in seq_along(k.norm)) {
        expect_equal(pval[[j]][[k]], cr.out[[i]][[j]][[k]]$pval[1, 2])
      }
    }
  })

  # 2. Two moments approach p-values
  test_that(sprintf("%s Two moments approach", test.name1), {
    for (i in seq_along(i.cores)) {
      j <- 2
      for (k in seq_along(k.norm)) {
        expect_equal(pval[[j]][[k]], cr.out[[i]][[j]][[k]]$pval[1, 2])
      }
    }
  })

  # 3. ci indicator
  test_that(sprintf("%s ci indicator", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(FALSE, cr.out[[i]][[j]][[k]]$ci)
        }
      }
    }
  })

  # 4. Test other common objects
  test.cr.common(cr.out, test.name, test.type)
}

# ---------------- #
# A list of unit tests for confidence interval
# ---------------- #
test.cr.ci <- function(cr.out, test.name, test.type) {
  # Assign the name
  test.name1 <- sprintf("'beta.obs' as %s and test for %s:",
                        test.name,
                        test.type)

  # 1. Full information approach p-values
  test_that(sprintf("%s Full information approach", test.name1), {
    for (i in seq_along(i.cores)) {
      j <- 1
      for (k in seq_along(k.norm)) {
        expect_equal(c.lb[[j]][[k]], cr.out[[i]][[j]][[k]]$ci.df[1, 3])
        expect_equal(c.ub[[j]][[k]], cr.out[[i]][[j]][[k]]$ci.df[1, 4])
      }
    }
  })

  # 2. Two moments approach p-values
  test_that(sprintf("%s Two moments approach", test.name1), {
    for (i in seq_along(i.cores)) {
      j <- 2
      for (k in seq_along(k.norm)) {
        expect_equal(c.lb[[j]][[k]], cr.out[[i]][[j]][[k]]$ci.df[1, 3])
        expect_equal(c.ub[[j]][[k]], cr.out[[i]][[j]][[k]]$ci.df[1, 4])
      }
    }
  })

  # 3. ci indicator
  test_that(sprintf("%s ci indicator", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(TRUE, cr.out[[i]][[j]][[k]]$ci)
        }
      }
    }
  })

  # 4. Significance level
  test_that(sprintf("%s Norm used", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(0.05, cr.out[[i]][[j]][[k]]$alpha)
        }
      }
    }
  })

  # 5. Test other common objects
  test.cr.common(cr.out, test.name, test.type)
}

# ---------------- #
# A list of unit tests that are common for confidence interval and p-value
# ---------------- #
test.cr.common <- function(cr.out, test.name, test.type) {
  # Assign the name
  test.name1 <- sprintf("'beta.obs' as %s and test for %s:",
                        test.name,
                        test.type)

  # 1. Lower bound
  test_that(sprintf("%s Lower bound", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(lb[[j]][[k]], cr.out[[i]][[j]][[k]]$lb[[1]])
        }
      }
    }
  })

  # 2. Upper bound
  test_that(sprintf("%s Upper bound", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(ub[[j]][[k]], cr.out[[i]][[j]][[k]]$ub[[1]])
        }
      }
    }
  })

  # 3. Test logical
  test_that(sprintf("%s Test logical", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(1, cr.out[[i]][[j]][[k]]$test.logical)
        }
      }
    }
  })

  # 4. df.error
  test_that(sprintf("%s df.error", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(NULL, cr.out[[i]][[j]][[k]]$df.error)
        }
      }
    }
  })

  # 5. Norm used
  test_that(sprintf("%s Norm used", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(k.norm[[k]], cr.out[[i]][[j]][[k]]$norm)
        }
      }
    }
  })

  # 6. Number of successful bootstrap replications
  test_that(sprintf("%s Number of successful bootstrap replications",
                    test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(reps, cr.out[[i]][[j]][[k]]$R.succ)
        }
      }
    }
  })

  # 7. kappa parameter
  test_that(sprintf("%s kappa parameter", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal(kappa, cr.out[[i]][[j]][[k]]$kappa)
        }
      }
    }
  })

  # 8. Solver name
  test_that(sprintf("%s Solver name", test.name1), {
    for (i in seq_along(i.cores)) {
      for (j in seq_along(j.lpmodel)) {
        for (k in seq_along(k.norm)) {
          expect_equal("gurobi", cr.out[[i]][[j]][[k]]$solver)
        }
      }
    }
  })
}

# ---------------- #
# Run the tests for p-value
# ---------------- #
# beta.obs is a function
test.cr.pval(cr.out.pval, "function", "p-value")

# beta.obs is a list of bootstrap estimates
test.cr.pval(cr.out.pval2, "list", "p-value")

# ---------------- #
# Run the tests for confidence interval
# ---------------- #
# beta.obs is a function
test.cr.ci(cr.out.ci, "function", "confidence interval")

# beta.obs is a list of bootstrap estimates
test.cr.ci(cr.out.ci2, "list", "confidence interval")

# ---------------- #
# Make sure the results obtained from the brute force approach is the
# same as the refinement approach
# ---------------- #
cr.out.brute.pval <- list()
cr.out.brute.ci <- list()

# p-value
farg$ci <- FALSE
farg$remove.const <- FALSE
for (j in seq_along(j.lpmodel)) {
  farg$lpmodel <- j.lpmodel[[j]]
  cr.out.brute.pval[[j]] <- list()
  for (k in seq_along(k.norm)) {
    set.seed(1)
    farg$norm <- k.norm[[k]]
    cr.out.brute.pval[[j]][[k]] <- do.call(chorussell, farg)
  }
}

# Confidence intervals
farg$ci <- TRUE
farg$remove.const <- FALSE
for (j in seq_along(j.lpmodel)) {
  farg$lpmodel <- j.lpmodel[[j]]
  cr.out.brute.ci[[j]] <- list()
  for (k in seq_along(k.norm)) {
    set.seed(1)
    farg$norm <- k.norm[[k]]
    cr.out.brute.ci[[j]][[k]] <- do.call(chorussell, farg)
  }
}

# Compare the confidence intervals
test_that("Same answers in brute force and refinement approaches", {
  for (i in seq_along(i.cores)) {
    for (j in seq_along(j.lpmodel)) {
      for (k in seq_along(k.norm)) {
        # p-value
        expect_equal(cr.out.brute.pval[[j]][[k]]$pval[1, 2],
                     cr.out.pval[[i]][[j]][[k]]$pval[1, 2])
        expect_equal(cr.out.brute.pval[[j]][[k]]$pval[1, 2],
                     cr.out.pval2[[i]][[j]][[k]]$pval[1, 2])

        # Lower-bound of the confidence intervals
        expect_equal(cr.out.brute.ci[[j]][[k]]$ci.df[1, 3],
                     cr.out.ci[[i]][[j]][[k]]$ci.df[1, 3])
        expect_equal(cr.out.brute.ci[[j]][[k]]$ci.df[1, 3],
                     cr.out.ci2[[i]][[j]][[k]]$ci.df[1, 3])
        
        # Upper-bound of the confidence intervals
        expect_equal(cr.out.brute.ci[[j]][[k]]$ci.df[1, 4],
                     cr.out.ci[[i]][[j]][[k]]$ci.df[1, 4])
        expect_equal(cr.out.brute.ci[[j]][[k]]$ci.df[1, 4],
                     cr.out.ci2[[i]][[j]][[k]]$ci.df[1, 4])
      }
    }
  }
})

# ---------------- #
# Make sure the results from `chorussell` and `invertci` are the same
# ---------------- #
set.seed(1)
cr.out.direct <- do.call(chorussell, farg)
set.seed(1)
cr.out.invertci <- invertci(f = chorussell,
                            farg = farg,
                            alpha = .05,
                            max.iter = 100,
                            tol = 1e-4)

# Same confidence interval
expect_equal(cr.out.direct$ci.df[1,], cr.out.invertci$ci[1,])
