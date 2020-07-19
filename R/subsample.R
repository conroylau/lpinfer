#' Computes the \eqn{p}-value of the subsampling procedure
#'
#' @description This function conducts inference and returns the
#'   \eqn{p}-value using the subsampling procedure.
#'
#' @import foreach doMC parallel
#'
#' @param lpmodel A list of objects that are used in inference of linear
#'    programming problems. The list of objects required in the \code{dkqs}
#'    procedure are:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{A.shp}}
#'      \item{\code{beta.obs}}
#'      \item{\code{beta.shp}}
#'    }
#'   matrix of the estimator \eqn{\hat{\beta}_{\mathrm{obs}}}.
#' @param phi Power for the sample. \eqn{n^\phi} represents the size
#'   of each subsample.
#' @param replace Boolean variable to indicate whether the function samples
#'   the data with or without replacement.
#' @inheritParams dkqs
#' @inheritParams estbounds
#'
#' @details There are three possible combinations for the parameters
#' \code{phi} and \code{replace}:
#' \itemize{
#'   \item{If \code{replace} is set as \code{FALSE}, it refers to the
#'     subsampling procedure. In this case, \code{phi} has to be in the
#'     interval \eqn{(0, 1)}.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} is set as 1,
#'     then it refers to the bootstrap procedure.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} is in the
#'     interval \eqn{(0, 1)}, then it refers to the \eqn{m} out of \eqn{n}
#'     bootstrap procedure, where \eqn{m} is the size of the subsample and
#'     \eqn{n} is the total number of observations.}
#' }
#'
#' @return Returns a list of output calculated from the function:
#'   \item{pval}{\eqn{p}-value.}
#'   \item{decision}{Decision of the test.}
#'   \item{T.n}{Test statistic \eqn{T_n}.}
#'   \item{T.sub}{The list of test statistics from the subsampling procedure.}
#'   \item{solver}{Solver used in solving the linear and quadratic programs.}
#'   \item{cores}{Number of cores used.}
#'   \item{cv.table}{Table of critical values.}
#'   \item{call}{The function that has been called.}
#'   \item{phi}{The \eqn{\phi} parameter used.}
#'   \item{norm}{Norm used.}
#'   \item{subsample.size}{Size of subsample}
#'   \item{test.logical}{Indicator variable for whether the computation has
#'     been conducted. If '\code{test.logical}' is 1, it refers to the case
#'     where '\code{beta.tgt}' is inside the logical bound. If
#'     '\code{test.logical}' is 0, it refers to the case where '
#'     \code{beta.tgt}' is outside the logical bound.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
subsample <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                      norm = 2, phi = 2/3, replace = FALSE, solver = NULL,
                      cores = 1, progress = TRUE){

  # ---------------- #
  # Step 1: Obtain call, check and update the dependencies
  # ---------------- #
  # Obtain the call information
  call <- match.call()

  # Check the arguments
  subsample.return <- subsample.check(data, lpmodel, beta.tgt, R, Rmulti,
                                      solver, cores, norm, phi, replace,
                                      progress)

  # Update the arguments
  data <- subsample.return$data
  lpmodel <- subsample.return$lpmodel
  solver <- subsample.return$solver
  solver.name <- subsample.return$solver.name
  cores <- subsample.return$cores
  norm <- subsample.return$norm
  test.logical <- subsample.return$test.logical

  # Compute size of each subsample
  n = nrow(data)
  m = floor(n^(phi))

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds
  if (test.logical == 1) {
    # = = = = = =
    # Step 2: Solve for T.n
    # = = = = = =
    ## Solve the main problem with the full sample
    Treturn <- subsample.prob(data, lpmodel, beta.tgt, norm, solver, 1)

    # ---------------- #
    # Step 3: Subsampling procedure
    # ---------------- #
    if (cores == 1){
      # One core
      T_subsample <- subsample.onecore(data, R, maxR, lpmodel, beta.tgt, norm,
                                       solver, replace, progress, m)
    } else {
      # Many cores
      T_subsample <- subsample.manycores(data, R, maxR, lpmodel, beta.tgt,
                                         norm, solver, cores, replace,
                                         progress, m)
    }

    # ---------------- #
    # Step 4: Compute the p-value (using the p_eval function in dkqs)
    # ---------------- #
    pval_return <- pval(T_subsample$T.sub, Treturn$objval)
    pval <- pval_return$p
    decision <- pval_return$decision

    # ---------------- #
    # Step 5: Generate a table of critical values
    # ---------------- #
    cv.table <- construct.cv.table(phi, "phi", Treturn$objval,
                                   T_subsample$T.sub)

    # ---------------- #
    # Step 6: Close the progress bar that is used in the subsampling procedure
    # ---------------- #
    if (progress == TRUE){
      close(T_subsample$pb)
      cat("                                            \n\b\r")
    }

    # ---------------- #
    # Step 7: Assign the return list
    # ---------------- #
    output <- list(pval = as.numeric(pval),
                   decision = decision,
                   T.n = as.numeric(Treturn$objval),
                   T.sub = T_subsample$T.sub,
                   solver = solver.name,
                   cores = cores,
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
                   cores = cores,
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

#' Formulate and solve the subsampling problem
#'
#' @description Based on the sample data given by the data frame \code{data},
#'   this function formulates and solves linear/quadratic program in the
#'   subsampling procedure.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @param i Index that represents whether the current problem is the
#'   bootstrap problem or the first-step problem.
#'
#' @return Returns a list of output that are obtained from the optimizer:
#'   \item{Status}{Status of the optimization problem.}
#'   \item{x}{Optimal point calculated from the optimizer.}
#'   \item{objval}{Optimal value calculated from the optimizer.}
#'   \item{larg}{List of arguments passed to the optimizer.}
#'   \item{beta}{The beta vector \eqn{\widehat{\bm{\beta}}_{\mathrm{obs}}}
#'     used in the optimization problem that is calculated from the
#'     \code{func_obs} function with data \code{data}.}
#'   \item{omega}{The Omega matrix \eqn{\widehat{\bm{\Omega}}_n} used in the
#'     optimization problem that is calculated from the \code{func_obs}
#'     function with data \code{data}.}
#'
#' @export
#'
subsample.prob <- function(data, lpmodel, beta.tgt, norm, solver, i){
  # ---------------- #
  # Step 1: Determine whether each argument is a function or a list
  # ---------------- #
  # beta.obs
  beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, i)
  beta.obs.hat <- beta.obs.return$beta.obs
  omega.hat <- beta.obs.return$omega

  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, i)
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, i)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, i)
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, i)
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
  if (norm == 1){
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
                  nf = 1,
                  A = A.new,
                  rhs = rhs.new,
                  sense = sense.new,
                  modelsense = modelsense.new,
                  lb = lb.new)
  } else if (norm == 2){
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
                  nf = sqrt(n),
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

#' Subsampling procedure without parallel programming
#'
#' @description This function computes the list of test statistics that are
#'   obtained from the subsampling procedure without using parallel
#'   programming.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @inheritParams subsample.prob
#' @inheritParams beta.bs
#' @param m Size of each subsample.
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T.sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
subsample.onecore <- function(data, R, maxR, lpmodel, beta.tgt, norm, solver,
                              replace, progress, m){
  # ---------------- #
  # Step 1: Initialize the vectors and the progress bar
  # ---------------- #
  # Initialize the vectors
  T.sub <- NULL
  beta.sub <- NULL

  # Initialize the progress bar
  if (progress == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = maxR, style = 3, width = 20)
    cat("\r")
  } else {
    pb <- NULL
  }

  # Initialize a table to contain the error messages
  df.error <- data.frame(matrix(vector(), ncol = 2))
  colnames(df.error) <- c("Iteration", "Error message")

  # ---------------- #
  # Step 2: Conduct the subsampling procedure
  # ---------------- #
  for (i in 1:maxR){
    # (2.1) Re-sample the data
    data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
    rownames(data.bs) <- 1:nrow(data.bs)

    # (2.2) Compute the bootstrap estimates
    result <- tryCatch(
      expr = {
        sub.return <- subsample.prob(data.bs, lpmodel, beta.tgt, norm, solver,
                                     i+1)
      },
      error = function(e) {
        return(list(status = "ERROR",
                    msg = e))
      },
      finally = {
        sub.return
      }
    )

    # (2.3) Store the results or error message depending on the status
    if (result$status %in% c("ERROR")) {
      df.error[nrow(df.error) + 1, 1] <- i
      df.error[nrow(df.error), 2] <- result$msg$message
    } else {
      T.sub <- c(T.sub, result$objval)
      beta.sub <- cbind(beta.sub, result$beta)
    }

    # (2.4) Update progress bar
    if (progress == TRUE){
      if ((i == maxR) | (length(T.sub) == R)){
        utils::setTxtProgressBar(pb, maxR)
        cat("\r\b")
      } else {
        utils::setTxtProgressBar(pb, i)
        cat("\r\r")
      }
    }

    # (2.5) Break the loop if R successful replications are made
    if (length(T.sub) == R) {
      if (progress == TRUE) {
        utils::setTxtProgressBar(pb, maxR)
        cat("\r\b")
      }
      break()
    }
  }

  # (2.6) Number of successful bootstrap replications
  R.succ <- length(T.sub)

  # (2.7) Set df.error as NULL if no failed bootstrap replications
  if (nrow(df.error) == 0) {
    df.error <- NULL
  }

  # ---------------- #
  # Step 3: Return the results
  # ---------------- #
  return(list(T.sub = T.sub,
              pb = pb,
              df.error = df.error,
              R.succ = R.succ))
}

#' Subsampling procedure with parallel programming
#'
#' @description This function computes the list of test statistics that are
#'   obtained from the subsampling procedure using parallel programming.
#'
#' @inheritParams dkqs
#' @inheritParams estbounds
#' @inheritParams subsample
#' @inheritParams subsample.onecore
#' @inheritParams beta.bs
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T.sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
subsample.manycores <- function(data, R, maxR, lpmodel, beta.tgt, norm, solver,
                                cores, replace, progress, m){
  # ---------------- #
  # Step 1: Initialize the parallel programming package
  # ---------------- #
  options(warn=-1)
  # Assign dorng
  `%dorng%` <- doRNG::`%dorng%`

  # Register core
  doMC::registerDoMC(cores)

  # ---------------- #
  # Step 2: Initialize the vectors and the progress bar
  # ---------------- #
  # Initialize the vectors
  T.sub <- NULL
  beta_sub <- NULL
  # Initialize the progress bar
  if (progress == TRUE) {
    # Initialize the counter
    cl <- PtProcess::makeSOCKcluster(8)
    doSNOW::registerDoSNOW(cl)

    # Set the counter and progress bar
    pb <- utils::txtProgressBar(max = maxR, style = 3, width = 20)

    cat("\r")
    progress <- function(n){
      utils::setTxtProgressBar(pb, n)
      cat("\r\r")
    }
    opts <- list(progress = progress)
  } else {
    pb <- NULL
    opts <- NULL
  }

  # Comb function for using parallel programming
  comb <- function(x, ...) {
    lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...),
                                                      function(y) y[[i]])))
  }

  # Initialize the lpmodel.bs and beta.obs.var object
  lpmodel.bs <- list()

  # ---------------- #
  # Step 3: Subsampling procedure
  # ---------------- #
  # Initialize the data frames
  k <- 0
  df.error1 <- data.frame(matrix(vector(), ncol = 2))
  colnames(df.error1) <- c("Iteration", "Error message")
  error.21 <- NULL
  error.22 <- NULL

  # Loop until the number of bootstrap replications match R or if maxR has been
  # reached
  while (k != R) {
    # Denote the starting index and ending index
    i0 <- k + 1
    i1 <- min(i0 + (R - k) - 1, maxR)

    # Use a normal for-loop to construct the list of lpmodel objects
    for (i in i0:i1) {
      # Construct bootstrap data
      data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
      rownames(data.bs) <- 1:nrow(data.bs)

      # Assign the lpmodel objects
      beta.obs.result <- tryCatch(
        expr <- {
          lpmodel.bs[[i]] <- list()
          lpmodel.bs[[i]]$A.obs <- lpmodel.eval(data.bs, lpmodel$A.obs, i + 1)
          lpmodel.bs[[i]]$A.shp <- lpmodel.eval(data.bs, lpmodel$A.shp, i + 1)
          lpmodel.bs[[i]]$A.tgt <- lpmodel.eval(data.bs, lpmodel$A.tgt, i + 1)
          lpmodel.bs[[i]]$beta.shp <- lpmodel.eval(data.bs, lpmodel$beta.shp,
                                                   i + 1)
          beta.obs.return <- lpmodel.beta.eval(data.bs, lpmodel$beta.obs, i + 1)
          lpmodel.bs[[i]]$beta.obs <- beta.obs.return
          beta.obs.ls <- list(status = "NOERROR",
                              lpmodel.bs = lpmodel.bs)
        },
        error = function(e) {
          return(list(status = "ERROR",
                      msg = e))
        },
        finally = {
          beta.obs.ls
        }
      )
      
      # Record error (if any)
      if (beta.obs.result$status == "ERROR") {
        df.error1[nrow(df.error1) + 1, 1] <- i
        df.error1[nrow(df.error1), 2] <- beta.obs.result$msg$message
      }
    }

    # Subsampling procedure
    listans = foreach::foreach(i = i0:i1, .multicombine = TRUE,
                               .combine = "comb", .options.snow = opts,
                               .packages = "lpinfer") %dorng%
      {
        # Only consider the subsample problem if there is no error in forming
        # the beta.obs parts
        if (!(i %in% df.error1[,1])) {
          ## (3.1) Compute the bootstrap estimates
          result <- tryCatch(
            expr = {
              sub.return <- subsample.prob(data.bs, lpmodel.bs[[i]], beta.tgt,
                                           norm, solver, 1)
            },
            error = function(e) {
              return(list(status = "ERROR",
                          msg = e))
            },
            finally = {
              sub.return
            }
          )

          ## (3.2) Store the results or error message depending on the status
          if (result$status %in% c("ERROR")) {
            ind <- i
            ind.msg <- result$msg$message
            T.sub <- NULL
            beta.sub <- NULL
          } else {
            ind <- NULL
            ind.msg <- NULL
            T.sub <- data.frame(result$objval)
            beta.sub <- data.frame(c(result$beta))
          }
        } else {
          ind <- NULL
          ind.msg <- NULL
          T.sub <- NULL
          beta.sub <- NULL
        }
        list(T.sub, beta.sub, ind, ind.msg)
      }

    ## (3.3) Extract the results
    T.sub.temp <- as.vector(unlist(listans[[1]]))
    beta.sub.temp <- data.frame(matrix(unlist(listans[[2]]),
                                       ncol = R,
                                       byrow = FALSE))

    ## (3.4) Combine with the previous results
    if (i0 == 1) {
      T.sub <- T.sub.temp
      beta.sub <- beta.sub.temp
    } else {
      T.sub <- c(T.sub, T.sub.temp)
      beta.sub <- cbind(beta.sub, beta.sub.temp)
    }
    k <- length(T.sub)

    ## (3.5) Consolidate the list of error messages
    if (length(unlist(listans[[3]])) != 0) {
      error.21.temp <- unlist(listans[[3]])
      error.22.temp <- unlist(listans[[4]])
      error.21 <- c(error.21, error.21.temp)
      error.22 <- c(error.22, error.22.temp)
    }

    ## (3.6) Break the while-loop if it reached maxR
    if (i1 == maxR) {
      break()
    }
  }
  # Close the progress bar
  cat("\r\b")

  # (3.6) Number of successful bootstrap replications
  R.succ <- length(T.sub)

  # ---------------- #
  # Step 4: Combine the error messages
  # ---------------- #
  error.length <- length(error.21)
  df.error2 <- data.frame(matrix(vector(), nrow = error.length, ncol = 2))
  colnames(df.error2) <- c("Iteration", "Error message")
  if (error.length != 0) {
    df.error2[,1] <- error.21
    df.error2[,2] <- error.22
  }

  df.error <- rbind(df.error1, df.error2)

  # ---------------- #
  # Step 5: Return the results
  # ---------------- #
  return(list(T.sub = T.sub,
              pb = pb,
              df.error = df.error,
              R.succ = R.succ))
}

#' Print results from \code{subsample}
#'
#' @description This function uses the print method on the return list of the
#'    function \code{subsample}.
#'
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.subsample <- function(x, ...){
  cat("\r\r")
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
#' @description This function uses the summary method on the return list of
#'    the function \code{subsample}.
#'
#' @param x Object returned from \code{subsample}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.subsample <- function(x, ...){
  cat("\r\r")
  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-values
    print(x)

    # Print test statistic, solver used, norm used, number of
    # cores used and the number of bootstrap replications
    cat(sprintf("Test statistic: %s\n", round(x$T.n, digits = 5)))
    cat(sprintf("Solver used: %s\n", x$solver))
    cat(sprintf("Norm used: %s\n", x$norm))
    cat(sprintf("Phi used: %s\n", round(x$phi, digits = 5)))
    cat(sprintf("Size of each subsample: %s\n", x$subsample.size))
    coremsg <- "Number of %s used: %s\n"
    ncore <- x$cores
    if (ncore == 1) {
      cat(sprintf(coremsg, "core", ncore))
    } else{
      cat(sprintf(coremsg, "cores", ncore))
    }
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

#' Checks and updates the input from \code{subsample}
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
#'       \item{\code{cores}}
#'       \item{\code{norm}}
#'       \item{\code{test.logical}}
#'    }
#'
#' @export
#'
subsample.check <- function(data, lpmodel, beta.tgt, R, Rmulti, solver, cores,
                            norm, phi, replace, progress){

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

  # Check other numbers
  check.positiveinteger(R, "R")
  cores <- check.cores(cores)

  # Check norm
  norm <- check.norm(norm, "norm")

  # Check Boolean
  check.boolean(replace, "replace")
  check.boolean(progress, "progress")

  # Check whether beta.tgt is within the logical bounds
  test.logical <- check.betatgt(data, lpmodel, beta.tgt, solver)

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              cores = cores,
              norm = norm,
              test.logical = test.logical))
}
