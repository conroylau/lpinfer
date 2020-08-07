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
#'   \item{T.sub}{List of bootstrap estimates of the test statistics
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
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
subsample <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                      norm = 2, phi = 2/3, replace = FALSE, solver = NULL,
                      progress = TRUE) {
  # ---------------- #
  # Step 1: Obtain call, check and update the dependencies
  # ---------------- #
  # Obtain the call information
  call <- match.call()

  # Check the arguments
  subsample.return <- subsample.check(data, lpmodel, beta.tgt, R, Rmulti,
                                      solver, norm, phi, replace, progress)

  # Update the arguments and seed
  data <- subsample.return$data
  lpmodel <- subsample.return$lpmodel
  solver <- subsample.return$solver
  solver.name <- subsample.return$solver.name
  norm <- subsample.return$norm
  test.logical <- subsample.return$test.logical
  seed <- subsample.return$seed

  # Compute size of each subsample
  n = nrow(data)
  m = floor(n^(phi))

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds
  if (test.logical == 1) {
    # ---------------- #
    # Step 2: Solve for T.n
    # ---------------- #
    ## Solve the main problem with the full sample
    Treturn <- subsample.prob(data, lpmodel, beta.tgt, norm, solver)

    # ---------------- #
    # Step 3: Subsampling procedure
    # ---------------- #
    T_subsample <- subsample.bs(data, R, maxR, lpmodel, beta.tgt, norm, solver,
                                replace, progress, m, seed)

    # ---------------- #
    # Step 4: Compute the p-value (using the p_eval function in dkqs)
    # ---------------- #
    pval_return <- pval(T_subsample$T.sub, Treturn$objval)
    pval <- pval_return$p

    # ---------------- #
    # Step 5: Generate a table of critical values
    # ---------------- #
    cv.table <- construct.cv.table(phi, "phi", Treturn$objval,
                                   T_subsample$T.sub)

    # ---------------- #
    # Step 6: Assign the return list
    # ---------------- #
    output <- list(pval = as.numeric(pval),
                   T.n = as.numeric(Treturn$objval),
                   T.sub = T_subsample$T.sub,
                   solver = solver.name,
                   cv.table = cv.table,
                   call = call,
                   phi = phi,
                   norm = norm,
                   subsample.size = m,
                   test.logical = test.logical,
                   df.error = T_subsample$df.error,
                   R.succ = T_subsample$R.succ)
  } else {
    ### Case 2: test.logical == 0. Set the p-value as 0 directly because
    ### beta.tgt is outside the logical bounds
    output <- list(pval = 0,
                   solver = solver.name,
                   call = call,
                   phi = phi,
                   norm = norm,
                   subsample.size = m,
                   test.logical = test.logical)

    # Print warning message
    infeasible.betatgt.warning()
  }

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
subsample.prob <- function(data, lpmodel, beta.tgt, norm, solver) {
  # ---------------- #
  # Step 1: Determine whether each argument is a function or a list
  # ---------------- #
  # beta.obs
  beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)
  beta.obs.hat <- beta.obs.return$beta.obs
  omega.hat <- beta.obs.return$omega

  # Always evaluate the first object in lpmodel - the 'lpmodel' object in the
  # bootstrap replications are explicitly passed
  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)
  k <- length(beta.obs.hat)

  # Count the number of rows
  n <- nrow(data)

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
#' @import future.apply
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @inheritParams subsample.prob
#' @inheritParams dkqs.bs
#' @param m Size of each subsample.
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
subsample.bs <- function(data, R, maxR, lpmodel, beta.tgt, norm, solver,
                         replace, progress, m, seed) {
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
    maxR <- length(any.list$consol[any.list$name[1]])
  }

  # ---------------- #
  # Step 2: Bootstrap replicatoins
  # ---------------- #
  while ((R.succ < R) & (R.eval != maxR)) {
    # Evaluate the list of indices to be passed to 'future_lapply'
    bs.temp <- bs.assign(R, R.eval, R.succ, maxR, any.list)
    i0 <- bs.temp$i0
    i1 <- bs.temp$i1
    bs.list <- bs.temp$bs.list

    # Obtain results from the bootstrap replications
    subsample.return <- future.apply::future_lapply(bs.list,
                                                    FUN = subsample.bs.fn,
                                                    future.seed = seed,
                                                    data = data,
                                                    lpmodel = lpmodel,
                                                    beta.tgt = beta.tgt,
                                                    norm = norm,
                                                    m = m,
                                                    solver = solver,
                                                    replace = replace)

    # Update the list and parameters
    post.return <- post.bs(subsample.return, i0, i1, R.eval, T.sub,
                           beta.sub.temp, error.list)
    T.sub <- post.return$T.list
    beta.sub.temp <- post.return$beta.list
    error.list <- post.return$error.list
    R.succ <- post.return$R.succ
    R.eval <- post.return$R.eval

    # Update seed
    seed <- .Random.seed
  }

  # ---------------- #
  # Step 3: Retrieve the bootstrap betas
  # ---------------- #
  if (R.succ != 0) {
    beta.sub <- data.frame(matrix(unlist(beta.sub.temp),
                                  nrow = R.succ,
                                  byrow = TRUE))
  } else {
    beta.sub <- NULL
  }

  # ---------------- #
  # Step 4: Consolidate the error messages
  # ---------------- #
  if (R.eval != R.succ) {
    # Create data.frame for error messages
    df.error <- data.frame(id = NA, message = unlist(error.list))

    # Remove 'NULL' from the list before passing to future_lapply
    df.error.nonnull <- error.list[-which(sapply(error.list, is.null))]

    # Match the id of the error messages
    df.error$id <- unlist(future.apply::future_lapply(df.error.nonnull,
                                                      FUN = match,
                                                      error.list))
  } else {
    df.error <- NULL
  }

  return(list(T.sub = T.sub,
              beta.sub = beta.sub,
              df.error = df.error,
              R.eval = R.eval,
              R.succ = R.succ))
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
                            replace) {
  # Replace lpmodel by x if x is a list
  if (is.list(x)) {
    lpm <- x
  } else {
    lpm <- lpmodel
  }

  # Draw data
  data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
  rownames(data.bs) <- 1:nrow(data.bs)

  # Bootstrap estimator
  result <- tryCatch({
    sub.return <- subsample.prob(data.bs, lpm, beta.tgt, norm, solver)
    sub.return
  }, warning = function(w) {
    return(list(status = "warning",
                msg = w))
  }, error = function(e) {
    return(list(status = "error",
                msg = e))
  })

  if (is.null(result$status)) {
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
    cat(sprintf("p-value: %s\n", round(x$pval, digits = 5)))
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
#'       \item{\code{seed}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{norm}}
#'       \item{\code{test.logical}}
#'    }
#'
#' @export
#'
subsample.check <- function(data, lpmodel, beta.tgt, R, Rmulti, solver, norm,
                            phi, replace, progress) {
  # ---------------- #
  # Step 1: Conduct the checks
  # ---------------- #
  # Check the data frame
  data <- check.dataframe(data)

  # Check lpmodel
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = c("matrix", "function_mat", "list"),
                           A.obs.cat = c("matrix", "function_mat", "list"),
                           A.shp.cat = c("matrix", "function_mat", "list"),
                           beta.obs.cat = c("list", "function_obs_var"),
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
  test.logical <- check.betatgt(data, lpmodel, beta.tgt, solver)

  # Obtain the seed
  seed <- lpinfer.seed()

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              norm = norm,
              test.logical = test.logical,
              seed = seed))
}
