#' Conducts inference using the subsampling procedure
#'
#' @description This function conducts inference and returns the
#'   \eqn{p}-value using the subsampling procedure.
#'
#' @param lpmodel The \code{lpmodel} object used in the test. The following
#'   components are required in the \code{lpmodel} for the subsampling test:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{A.shp}}
#'      \item{\code{beta.obs}}
#'      \item{\code{beta.shp}}
#'    }
#'   matrix of the estimator \eqn{\hat{\bm{\beta}}_{\mathrm{obs}}}.
#' @param phi Tuning parameter for the subsampling test. The size of each
#'   subsample is \eqn{n^\phi} where \eqn{\phi \in [0,1]}.
#' @param replace Boolean variable to indicate whether the function samples
#'   the data with or without replacement.
#' @inheritParams dkqs
#' @inheritParams estbounds
#'
#' @details There are three possible combinations for the parameters
#' \code{phi} and \code{replace}:
#' \itemize{
#'   \item{If \code{replace} is set as \code{FALSE}, it refers to the
#'     subsampling procedure. In this case, \code{phi} has to be inside the
#'     interval \eqn{(0, 1)}.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} is set as 1,
#'     then it refers to the bootstrap procedure.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} is inside the
#'     interval \eqn{(0, 1)}, then it refers to the \eqn{m} out of \eqn{n}
#'     bootstrap procedure, where \eqn{m} is the size of the subsample and
#'     \eqn{n} is the total number of observations.}
#' }
#'
#' @return Returns a list of output calculated from the function:
#'   \item{pval}{\eqn{p}-value.}
#'   \item{T.n}{Test statistic \eqn{T_n}.}
#'   \item{T.bs}{List of bootstrap estimates of the test statistics
#'     from the subsampling procedure.}
#'   \item{solver}{Solver used.}
#'   \item{cv.table}{Table of critical values.}
#'   \item{call}{The function that has been called.}
#'   \item{phi}{The \eqn{\phi} parameter used.}
#'   \item{norm}{Norm used.}
#'   \item{subsample.size}{Size of subsample}
#'   \item{test.logical}{Indicator variable for whether the computation has
#'     been conducted. If \code{test.logical} is 1, it refers to the case
#'     where \code{beta.tgt} is inside the logical bound. If
#'     \code{test.logical} is 0, it refers to the case where
#'     \code{beta.tgt} is outside the logical bound.}
#'   \item{logical.lb}{Logical lower bound.}
#'   \item{logical.ub}{Logical upper bound.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @example ./inst/example/subsample_example.R
#'
#' @export
#'
subsample <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                      norm = 2, phi = 2/3, n = NULL, replace = FALSE,
                      solver = NULL, progress = TRUE) {
  # ---------------- #
  # Step 1: Obtain call, check and update the dependencies
  # ---------------- #
  # Extract the current RNG state
  rngstate <- .Random.seed

  # Obtain the call information
  call <- match.call()

  # Check the arguments
  subsample.return <- subsample.check(data, lpmodel, beta.tgt, R, Rmulti,
                                      solver, norm, phi, n, replace, progress)

  # Update the arguments
  data <- subsample.return$data
  lpmodel <- subsample.return$lpmodel
  solver <- subsample.return$solver
  solver.name <- subsample.return$solver.name
  norm <- subsample.return$norm
  test.logical <- subsample.return$test.logical
  logical.lb <- subsample.return$logical.lb
  logical.ub <- subsample.return$logical.ub

  # Compute size of each subsample
  if (!is.null(data)) {
    n <- nrow(data)
  }
  m <- floor(n^(phi))

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  # List of RNG states and the corresponding number of iterations
  seed.list <- list()
  seed.list[[1]] <- list(seed = rngstate,
                         iter = R)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds
  if (test.logical == 1) {
    # Evaluate the beta.obs and variance from lpmodel
    beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)
    beta.obs.hat <- beta.obs.return$beta.obs
    omega.hat0 <- beta.obs.return$omega

    # Initialize parameters
    R.succ <- -1
    i1 <- -1
    error.id <- NULL
    new.error.bs <- 0
    beta.obs.bs <- list()
    beta.obs.full <- list()

    # Initialize a table to contain the error messages
    df.error <- data.frame(matrix(vector(), ncol = 3))
    colnames(df.error) <- c("Iteration", "phi", "Error message")

    # Start for-loop
    while (((R.succ < R) & (i1 < maxR)) | new.error.bs != 0) {
      # Assign the indices
      if (R.succ == -1) {
        # -1 corresponds to the initial bootstrap replications
        i0 <- 1
        i1 <- R
        iseq <- 1:maxR
        eval.count <- 0
      } else {
        # Update the sequence of indices
        i0 <- min(maxR, i1 + 1)
        i1 <- min(maxR, i0 + (R - R.succ) - 1)
        if (class(lpmodel$beta.obs) == "list") {
          i0 <- i1 + 1
          i1 <- i0 + (R - R.succ)
        }
        iseq <- i0:i1
        eval.count <- eval.count + 1
        seed.list[[length(seed.list)]]$iter <- i1 - i0 + 1
      }

      # ---------------- #
      # Step 3: Bootstrap variance matrix if it is not provided
      # ---------------- #
      if (is.null(omega.hat0)) {
        ## Case 1: Variance matrix is not provided by the user
        assign(x = ".Random.seed", value = rngstate, envir = .GlobalEnv)
        # Write a new function just for computing the beta bs
        beta.obs.return <- beta.bs(data, lpmodel, seed.list, df.error)

        # Store the new RNG state
        seed.list[[length(seed.list) + 1]] <- list(seed = .Random.seed)

        # Obtain the list of beta.obs
        beta.obs.bs.new <- beta.obs.return$beta.obs.bs
        beta.obs.full <- c(beta.obs.full, beta.obs.bs.new)

        # Consolidate the error messages (if any)
        df.error <- beta.obs.return$df.error
        error.id.new <- beta.obs.return$error.id
        error.id <- c(error.id, error.id.new)
        R.succ <- beta.obs.return$R.succ

        # Next if there is some problematic draws
        new.error.bs <- beta.obs.return$R.eval - R.succ
        if (new.error.bs != 0) {
          next
        }

        # Remove the problematic draws in computing variance matrix
        if (!is.null(error.id)) {
          beta.obs.bs <- beta.obs.full[-c(error.id)]
        } else {
          beta.obs.bs <- beta.obs.full
        }
        beta.obs.list <- c(list(beta.obs.hat), beta.obs.bs)
        omega.hat <- sigma.summation(n, beta.obs.list, progress, eval.count)
      } else {
        ## Case 2: Variance matrix is provided by the user
        omega.hat <- omega.hat0
      }

      # ---------------- #
      # Step 3: Solve for T.n
      # ---------------- #
      # Solve the main problem with the full sample
      Treturn <- subsample.prob(data, lpmodel, beta.tgt, norm, solver, n,
                                beta.obs.hat, omega.hat)

      # ---------------- #
      # Step 4: Subsampling procedure
      # ---------------- #
      # Solve the subsampling problem
      Tsub.return <- subsample.bs(data, i1, lpmodel, beta.tgt, norm, solver,
                                  replace, progress, m, n, omega.hat, df.error,
                                  eval.count, error.id, seed.list)

      T.sub <- Tsub.return$T.sub
      df.error <- Tsub.return$df.error
      error.id <- Tsub.return$error.id
      R.eval <- Tsub.return$R.eval
      R.succ <- Tsub.return$R.succ

      # Next if there is some problematic draws
      if (Tsub.return$new.error != 0) {
        next
      }
    }

    if (R.succ != 0) {
      # ---------------- #
      # Step 5: Compute the p-value (using the p_eval function in dkqs)
      # ---------------- #
      # # Obtain the result
      pval_return <- pval(T.sub, Treturn$objval)

      # Create a data frame of p-values
      pval.df <- data.frame(matrix(vector(), nrow = 1, ncol = 2))
      colnames(pval.df) <- c("phi", "p-value")
      pval.df[1, 1] <- phi
      pval.df[1, 2] <- pval_return$p

      # ---------------- #
      # Step 6: Generate a table of critical values
      # ---------------- #
      cv.table <- construct.cv.table(phi, "phi", Treturn$objval,
                                     T.sub)

    } else {
      pval <- NA
      cv.table <- NA
    }

    # ---------------- #
    # Step 7: Assign the return list
    # ---------------- #
    output <- list(pval = pval.df,
                   T.n = as.numeric(Treturn$objval),
                   T.bs = T.sub,
                   cv.table = cv.table,
                   df.error = df.error,
                   R.succ = R.succ)
  } else {
    ### Case 2: test.logical == 0. Set the p-value as 0 directly because
    ### beta.tgt is outside the logical bounds
    output <- list(pval = 0)

    # Print warning message
    infeasible.betatgt.warning()
  }
  # Assign the common objects in the output list
  output <- append(output,
                   list(solver = solver.name,
                        call = call,
                        phi = phi,
                        norm = norm,
                        subsample.size = m,
                        test.logical = test.logical,
                        logical.lb = logical.lb,
                        logical.ub = logical.ub))

  # Assign class
  attr(output, "class") <- "subsample"

  return(output)
}

#' Formulates and solves the subsampling problem
#'
#' @description This function formulates and solves the linear or quadratic
#'   program in the subsampling procedure. If the user chooses a 1-norm, this
#'   function solves a linear program. If the user chooses a 2-norm, this
#'   function solves a quadratic program.
#'
#' @inheritParams subsample
#' @inheritParams dkqs.qlp
#' @param omega.hat Estimator of the asymptotic variance.
#'
#' @return Returns the following list of outputs:
#'   \item{status}{Status of the optimization problem.}
#'   \item{x}{Optimal point.}
#'   \item{objval}{Optimal objective value.}
#'   \item{larg}{List of arguments passed to the optimizer.}
#'   \item{beta}{The beta vector \eqn{\widehat{\bm{\beta}}_{\mathrm{obs}}}
#'     used in the optimization problem that is obtained from the
#'     \code{beta.obs} component of the \code{lpmodel} object.}
#'   \item{omega}{The Omega matrix \eqn{\widehat{\bm{\Omega}}_n} used in the
#'     optimization problem that is obtained from the \code{beta.obs}
#'     component of the \code{lpmodel} object.}
#'
#' @export
#'
subsample.prob <- function(data, lpmodel, beta.tgt, norm, solver, n,
                           beta.obs.hat, omega.hat) {
  # ---------------- #
  # Step 1: Determine whether each argument is a function or a list
  # ---------------- #
  # Always evaluate the first object in lpmodel - the 'lpmodel' object in the
  # bootstrap replications are explicitly passed
  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)
  k <- length(beta.obs.hat)

  # ---------------- #
  # Step 2: Define the inverse omega matrix
  # ---------------- #
  # Obtain the inverse of the diagonal entries
  diag.omega <- diag(omega.hat)
  g <- 1/diag.omega
  # Replace the entries by 0 for those that are equal to zero in 1/Omega
  g[diag.omega == 0] <- 0
  G <- diag(g)
  # Create the new A and b matrices
  GA <- G %*% A.obs.hat
  Gb <- G %*% beta.obs.hat

  # ---------------- #
  # Step 3: Form the objective function and constraints
  # ---------------- #
  # Model sense
  modelsense.new <- "min"

  # Set the objective function and constraints
  if (norm == 1) {
    ### L1-norm
    # Objective function - cost matrix
    c.new <- c(rep(0, ncol(A.obs.hat)), rep(1, k), rep(-1, k))

    # Constraints
    A.zero.shp <- matrix(rep(0, k*nrow(A.shp.hat)), nrow = nrow(A.shp.hat))
    A.zero.tgt <- matrix(rep(0, k*nrow(A.tgt.hat)), nrow = nrow(A.tgt.hat))
    A1.shp <- cbind(A.shp.hat, A.zero.shp, A.zero.shp)
    A1.tgt <- cbind(A.tgt.hat, A.zero.tgt, A.zero.tgt)
    A1.obs <- cbind(GA, -diag(k), diag(k))
    A.new <- rbind(A1.shp, A1.tgt, A1.obs)

    # RHS vector
    rhs.new <- c(beta.shp.hat, beta.tgt, Gb)

    # Lower bounds
    lb.new <- rep(0, length(c.new))

    # Sense
    sense.new <- rep("=", nrow(A.new))

    # Set the list to pass to the solver
    l_arg <- list(Af = NULL,
                  bf = c.new,
                  nf = sqrt(n),
                  A = A.new,
                  rhs = rhs.new,
                  sense = sense.new,
                  modelsense = modelsense.new,
                  lb = lb.new)
  } else if (norm == 2) {
    ### L2-norm
    # Constraints
    A.new <- rbind(A.shp.hat, A.tgt.hat)

    # RHS vector
    rhs.new <- c(beta.shp.hat, beta.tgt)

    # Lower bounds
    lb.new <- rep(0, ncol(A.shp.hat))

    # Sense
    sense.new <- rep("=", nrow(A.new))

    # Set the list to pass to the solver
    l_arg <- list(Af = GA,
                  bf = Gb,
                  nf = n,
                  A = A.new,
                  rhs = rhs.new,
                  sense = sense.new,
                  modelsense = modelsense.new,
                  lb = lb.new)
  }

  # ---------------- #
  # Step 4: Solve the model and return the results
  # ---------------- #
  # Solve the model
  ans <- do.call(solver, l_arg)

  invisible(list(status = ans$status,
                 x = ans$x,
                 objval = ans$objval,
                 larg = l_arg,
                 beta = beta.obs.hat,
                 omega = omega.hat))
}

#' Bootstrap procedure for the subsampling test
#'
#' @description This function carries out the bootstrap procedure of the
#'   subsampling test. This function supports parallel programming via the
#'   \code{future.apply} package.
#'
#' @import future.apply progressr
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @inheritParams subsample.prob
#' @inheritParams dkqs.bs
#' @param m Size of each subsample.
#' @param error.id List of ID that corresponds to problematic bootstrap draws.
#' @param seed.list List of RNG states and the corresponding iterations.
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T.sub}{Bootstrap test statistics from the subsampling procedure.}
#'   \item{beta.sub}{Bootstrap estimators for the \code{beta} component.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.eval}{Number of bootstrap replications that have been conducted.}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
subsample.bs <- function(data, i1, lpmodel, beta.tgt, norm, solver,
                         replace, progress, m, n, omega.hat, df.error,
                         eval.count, error.id, seed.list) {
  # ---------------- #
  # Step 1: Initialize the quantities and compute the lists
  # ---------------- #
  R.succ <- 0
  R.eval <- 0
  T.sub <- list()
  beta.sub.temp <- list()
  error.list <- list()

  # Check if there is any list objects in 'lpmodel'
  any.list <- lpmodel.anylist(lpmodel)

  # If there is some list objects, set maxR as the max length of the list
  if (isTRUE(any.list$list)) {
    maxR <- length(any.list$consol)
  }

  # ---------------- #
  # Step 2: Subsampling
  # ---------------- #
  iter.list <- unlist(sapply(seed.list, "[", "iter"))
  for (i in seq_along(iter.list)) {
    # Set seed
    assign(x = ".Random.seed", value = seed.list[[i]]$seed, envir = .GlobalEnv)

    # Set the indices
    bs.temp <- bs.assign(seed.list[[i]]$iter, R.eval, R.succ,
                         seed.list[[i]]$iter, any.list, lpmodel, data, m,
                         replace)
    i0 <- bs.temp$i0
    i1 <- bs.temp$i1
    bs.list <- bs.temp$bs.list

    # Set the default for progress bar
    progressr::handlers("progress")

    # Obtain results from the bootstrap replications
    progressr::with_progress({
      if (isTRUE(progress)) {
        pbar <- progressr::progressor(along = i0:i1)
      } else {
        pbar <- NULL
      }

      # Obtain results from future_lapply
      subsample.return <- future.apply::future_lapply(bs.list,
                                                      FUN = subsample.bs.fn,
                                                      future.seed = TRUE,
                                                      data = data,
                                                      lpmodel = lpmodel,
                                                      beta.tgt = beta.tgt,
                                                      norm = norm,
                                                      m = m,
                                                      solver = solver,
                                                      replace = replace,
                                                      n = n,
                                                      pbar = pbar,
                                                      progress = progress,
                                                      eval.count = eval.count,
                                                      n.bs = i1 - i0 + 1,
                                                      omega.hat = omega.hat)
    })

    # Extract the test satistics and list of errors
    T.sub <- c(T.sub, sapply(subsample.return, "[", "Ts"))
    error.list <- c(error.list, sapply(subsample.return, "[", "msg"))
  }

  # ---------------- #
  # Step 3: Consolidate the error messages
  # ---------------- #
  # Count the previous number of errors
  if (is.null(df.error)) {
    error0 <- 0
  } else {
    error0 <- nrow(df.error)
  }

  ns <- length(iter.list)
  if (length(unlist(error.list)) != 0) {
    # New error messages
    new.ind <- (sum(iter.list[1:(ns - 1)]) + 1):(sum(iter.list))
    df.error1 <- data.frame(id = NA,
                            lambda = NA,
                            message = unlist(error.list[new.ind]))
    df.error1 <- error.id.match(error.list[new.ind], df.error1)

    # Merge with the previous error messages
    df.error <- rbind(df.error, df.error1)

    # New errors
    new.error <- nrow(df.error) - error0

    # Consolidate the error ids and get the list of nonproblematic test
    # statistics
    error.id1 <- df.error$id
    error.id <- unique(c(error.id, error.id1))
    T.sub <- unlist(T.sub[-c(error.id)], use.names = FALSE)
  } else {
    df.error <- NULL
    T.sub <- unlist(T.sub, use.names = FALSE)
    error.id <- error.id
    new.error <- 0
  }

  # Compute the number of evaluations and successful evaluations
  R.succ <- length(T.sub)
  R.eval <- i1

  return(list(T.sub = T.sub,
              df.error = df.error,
              error.id = error.id,
              R.eval = R.eval,
              R.succ = R.succ,
              new.error = new.error))
}

#' Carries out one bootstrap replication for the subsampling test
#'
#' @description This function carries out the one bootstrap replication of the
#'   subsampling test. This function is used in the \code{subsample.bs}
#'   function via the \code{future_lapply} command.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @inheritParams subsample.prob
#' @inheritParams subsample.bs
#' @inheritParams dkqs.bs
#' @inheritParams dkqs.bs.fn
#' @param x This is either the list of indices that represent the bootstrap
#'   replications, or the list of bootstrap components of the \code{lpmodel}
#'   object passed from the user.
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{Ts}{Bootstrap test statistic.}
#'   \item{beta}{Bootstrap estimate of \code{beta.obs}.}
#'   \item{msg}{Error message (if applicable).}
#'
#' @export
#'
subsample.bs.fn <- function(x, data, lpmodel, beta.tgt, norm, m, solver,
                            replace, n, pbar, eval.count, n.bs, progress,
                            omega.hat) {
  # ---------------- #
  # Step 1: Print progress bar
  # ---------------- #
  # If replace = FALSE, it refers to subsampling. Otherwise, it is referring
  # to bootstrap or m out of n bootstrap
  if (isTRUE(replace)) {
    procedure <- "bootstrap"
  } else {
    procedure <- "subsample"
  }

  # Print progress bar
  if (isTRUE(progress)) {
    if (eval.count == 0) {
      pbar(sprintf("(Computing %s %s estimates)", n.bs, procedure))
    } else {
      pbar(sprintf("(Computing %s extra %s estimates)", n.bs, procedure))
    }
  }

  # ---------------- #
  # Step 2: Initialize the parameters
  # ---------------- #
  # Replace lpmodel by x if x is a list
  if (is.list(x)) {
    lpm <- lpmodel.update(lpmodel, x)
  } else {
    lpm <- lpmodel
  }

  # ---------------- #
  # Step 3: Conduct one bootstrap/subsample replication
  # ---------------- #
  # Draw data
  if (!is.null(data)) {
    data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
    rownames(data.bs) <- 1:nrow(data.bs)
  } else {
    data.bs <- NULL
  }

  # Bootstrap estimator
  result <- tryCatch({
    # The 'lpm' object does not contain the asymptotic variance estimator
    beta.obs.hat <- lpmodel.beta.eval(data.bs, lpm$beta.obs, 1)$beta.obs
    sub.return <- subsample.prob(data.bs, lpm, beta.tgt, norm, solver, m,
                                 beta.obs.hat, omega.hat)
    sub.return
  }, warning = function(w) {
    return(list(status = "warning",
                msg = w))
  }, error = function(e) {
    return(list(status = "error",
                msg = e))
  })

  if (!(result$status %in% c("warning", "error"))) {
    Ts <- result$objval
    beta <- result$beta
    msg <- NULL
  } else {
    Ts <- NULL
    beta <- NULL
    msg <- result$msg$message
  }

  return(list(Ts = Ts,
              beta = beta,
              msg = msg))
}

#' Print results from \code{subsample}
#'
#' @description This function prints the \eqn{p}-values from \code{subsample}.
#'
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.subsample <- function(x, ...) {
  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    cat(sprintf("p-value: %s\n", round(x$pval[1, 2], digits = 5)))
  } else {
    # Case 2: 'beta.tgt' is outside the logical bound
    infeasible.pval.msg()
  }
}

#' Summary of results from \code{subsample}
#'
#' @description This function prints a summary of the results obtained from
#'   \code{subsample}.
#'
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.subsample <- function(x, ...) {
  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-values
    print(x)

    # Print test statistic, solver used, norm used, and the number of
    # successful bootstrap replications
    cat(sprintf("Test statistic: %s\n", round(x$T.n, digits = 5)))
    cat(sprintf("Solver used: %s\n", x$solver))
    cat(sprintf("Norm used: %s\n", x$norm))
    cat(sprintf("Phi used: %s\n", round(x$phi, digits = 5)))
    cat(sprintf("Size of each subsample: %s\n", x$subsample.size))
    cat(sprintf("Number of successful bootstrap replications: %s\n", x$R.succ))

    # Number of failed bootstrap replications
    if (!is.null(x$df.error)) {
      nerr <- nrow(x$df.error)
      errstring <- "Number of failed bootstrap"
      if (nerr == 1) {
        cat(sprintf(paste(errstring, "replication: %s\n"), nerr))
      } else {
        cat(sprintf(paste(errstring, "replications: %s\n"), nerr))
      }
    }
  } else if (x$test.logical == 0) {
    # Case 2: 'beta.tgt' is outside the logical bound
    infeasible.pval.msg()
    cat(sprintf("\nSolver used: %s\n", x$solver))
  }
}

#' Checks and updates the input in \code{subsample}
#'
#' @description This function checks and updates the input of the user. If
#'    there is any invalid input, this function will be terminated and
#'    generates appropriate error messages. This function is mainly a wrapper
#'    of the selected functions from the \code{checks} files to conduct the
#'    checks and updates.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#'
#' @return Returns the updated parameters and objects back to the function
#' \code{subsample}. The following information are updated:
#'    \itemize{
#'       \item{\code{data}}
#'       \item{\code{lpmodel}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{norm}}
#'       \item{\code{test.logical}}
#'       \item{\code{logical.lb}}
#'       \item{\code{logical.ub}}
#'    }
#'
#' @export
#'
subsample.check <- function(data, lpmodel, beta.tgt, R, Rmulti, solver, norm,
                            phi, n, replace, progress) {
  # ---------------- #
  # Step 1: Conduct the checks
  # ---------------- #
  # Check data. If data is NULL, check if n is a positive integer
  if (!is.null(data)) {
    data <- check.dataframe(data)
  } else {
    check.samplesize(n, "n")
  }

  # Check lpmodel
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = c("matrix", "function_mat", "list"),
                           A.obs.cat = c("matrix", "function_mat", "list"),
                           A.shp.cat = c("matrix", "function_mat", "list"),
                           beta.obs.cat = c("list", "function_obs_var_bs"),
                           beta.shp.cat = c("matrix", "function_mat", "list"),
                           R = R)

  # Check solver
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # Check phi
  if (isTRUE(replace)) {
    check.numrange(phi, "phi", "open", 0, "closed", 1)
  } else if (isFALSE(replace)) {
    check.numrange(phi, "phi", "open", 0, "open", 1)
  }

  # Check Rmulti
  check.numrange(Rmulti, "Rmulti", "closed", 1, "open", Inf)

  # Check the number of bootstrap replications
  check.positiveinteger(R, "R")

  # Check norm
  norm <- check.norm(norm, "norm")

  # Check Boolean
  check.boolean(replace, "replace")
  check.boolean(progress, "progress")

  # Check whether beta.tgt is within the logical bounds
  check.numeric(beta.tgt, "beta.tgt")
  test.return <- check.betatgt(data, lpmodel, beta.tgt, solver)
  test.logical <- test.return$inout
  logical.lb <- test.return$lb
  logical.ub <- test.return$ub

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              norm = norm,
              test.logical = test.logical,
              logical.lb = logical.lb,
              logical.ub = logical.ub))
}

#' Evaluate the bootstrap betas
#'
#' @description This function computes the bootstrap estimates of
#'   \eqn{\beta_{\rm obs}}.
#'
#' @inheritParams dkqs
#' @inheritParams subsample
#' @inheritParams subsample.bs
#'
#' @return Returns the following objects:
#'   \item{beta.obs.bs}{Bootstrap estimates of \eqn{\beta_{\rm obs}}.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{error.id}{List of problematic IDs.}
#'   \item{R.eval}{Number of bootstrap replications that have been conducted.}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
beta.bs <- function(data, lpmodel, seed.list, df.error) {
  # Initialize
  beta.list <- list()
  msg.list <- list()

  # Assign the RNG state
  assign(x = ".Random.seed",
         value = seed.list[[length(seed.list)]]$seed,
         envir = .GlobalEnv)

  ind <- unlist(sapply(seed.list, "[", "iter"), use.names = FALSE)

  # Compute the bootstrap betas
  for (i in 1:ind[length(seed.list)]) {
    data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])

    # Compute the bootstrap beta
    result <- tryCatch({
      beta.obs.return <- lpmodel.beta.eval(data.bs, lpmodel$beta.obs, 1)
      list(beta.obs.return = beta.obs.return)
    }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
    }, error = function(e) {
      return(list(status = "error",
                  msg = e))
    })

    # Assign either the beta or the error message
    if (is.null(result$status)) {
      beta.list[[i]] <- result$beta.obs.return$beta.obs
      msg.list[[i]] <- NULL
    } else {
      beta.list[[i]] <- NULL
      msg.list[[i]] <- result$msg$message
    }
  }

  # Other parameters
  R.eval <- ind[length(seed.list)]
  R.succ <- R.eval - length(unlist(msg.list))

  # ---------------- #
  # Step 3: Consolidate the error messages
  # ---------------- #
  if (R.eval != R.succ) {
    # Create data.frame for error messages
    df.error1 <- data.frame(id = NA,
                            lambda = NA,
                            message = unlist(msg.list))

    # Match the id of the error messages
    df.error1 <- error.id.match(msg.list, df.error1)
    error.id <- df.error1$id + sum(ind) - ind[length(seed.list)]
    if (!is.null(df.error)) {
      # Case 1: There are errors in this function and the earlier parts
      df.error <- rbind(df.error, df.error1)
      rownames(df.error) <- 1:nrow(df.error)
    } else {
      # Case 2: There are no errors in the earlier steps and there are
      # errors in this step
      df.error <- df.error1
    }
    new.error <- R.eval - R.succ
  } else {
    # Case 3: There are no new errors in this procedure
    df.error <- df.error
    error.id <- NULL
  }

  return(list(beta.obs.bs = beta.list,
              df.error = df.error,
              error.id = error.id,
              R.eval = R.eval,
              R.succ = R.succ))
}
