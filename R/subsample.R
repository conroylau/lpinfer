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
#' @param alpha Significance level.
#' @param replace Boolean variable to indicate whether the function samples
#'   the data with or without replacement.
#' @inheritParams dkqs
#' @inheritParams estbounds
#' 
#' @details There are three possible cases for the parameter \code{replace}:
#' \itemize{
#'   \item{If \code{replace} is set as \code{FALSE}, it refers to the 
#'     subsampling procedure.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} is set as 1,
#'     then it refers to the bootstrap procedure.}
#'   \item{If \code{replace} is set as \code{TRUE} and \code{phi} set between
#'     0 and 1, then it refers to the \eqn{m} out of \eqn{n} bootstrap 
#'     procedure, where \eqn{m} is the size of the subsample and \eqn{n} is
#'     the total number of observations.}   
#' }
#'
#' @return Returns a list of output calculated from the function:
#'   \item{pval}{\eqn{p}-value.}
#'   \item{decision}{Decision of the test.}
#'   \item{alpha}{Significance level.}
#'   \item{T.n}{Test statistic \eqn{T_n}.}
#'   \item{T.sub}{The list of test statistics from the subsampling procedure.}
#'   \item{solver}{Solver used in solving the linear and quadratic programs.}
#'   \item{cores}{Number of cores used.}
#'   \item{call}{The function that has been called.}
#'   \item{phi}{The \eqn{\phi} parameter used.}
#'   \item{norm}{Norm used.}
#'   \item{subsample.size}{Size of subsample}
#'
#' @export
#'
subsample <- function(data, lpmodel, beta.tgt, R = 100, solver = NULL,
                      cores = 8, norm = 2, phi = 2/3, alpha = .05,
                      replace = FALSE, progress = FALSE){

  # ---------------- #
  # Step 1: Obtain call, check and update the dependencies
  # ---------------- #
  # Obtain the call information
  call <- match.call()

  # Check the arguments
  subsample.return <- subsample.check(data, lpmodel, beta.tgt, R, solver, 
                                      cores, norm, phi, alpha, replace,
                                      progress)

  # Update the arguments
  data <- subsample.return$data
  lpmodel <- subsample.return$lpmodel
  solver <- subsample.return$solver
  solver.name <- subsample.return$solver.name
  cores <- subsample.return$cores
  norm <- subsample.return$norm

  # = = = = = =
  # Step 2: Solve for T.n
  # = = = = = =
  ## Solve the main problem with the full sample
  Treturn <- subsample.prob(data, lpmodel, beta.tgt, norm, solver, 1)

  # ---------------- #
  # Step 3: Subsampling procedure
  # ---------------- #
  n = nrow(data)
  m = floor(n^(phi))

  if (cores == 1){
    # One core
    T_subsample <- subsample.onecore(data, R, lpmodel, beta.tgt, norm,
                                     solver, replace, progress, m)

  } else {
    # Many cores
    T_subsample <- subsample.manycores(data, R, lpmodel, beta.tgt, norm,
                                       solver, cores, replace, progress, m)
  }

  # ---------------- #
  # Step 4: Compute the p-value (using the p_eval function in dkqs)
  # ---------------- #
  pval_return <- pval(T_subsample$T.sub, Treturn$objval, alpha)
  pval <- pval_return$p
  decision <- pval_return$decision

  # ---------------- #
  # Step 5: Close the progress bar that is used in the subsampling procedure
  # ---------------- #
  if (progress == TRUE){
    close(T_subsample$pb)
    cat("                                            \n\b\r")
  }

  # ---------------- #
  # Step 6: Assign the return list
  # ---------------- #
  output <- list(pval = as.numeric(pval),
                 decision = decision,
                 alpha = alpha,
                 T.n = as.numeric(Treturn$objval),
                 T.sub = T_subsample$T.sub,
                 solver = solver.name,
                 cores = cores,
                 call = call,
                 phi = phi,
                 norm = norm,
                 subsample.size = m)

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

  # Return the results
  invisible(list(x = ans$x,
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
#' @param m Size of each subsample.
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T.sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#'
#' @export
#'
subsample.onecore <- function(data, R, lpmodel, beta.tgt, norm, solver,
                              replace, progress, m){
  # ---------------- #
  # Step 1: Initialize the vectors and the progress bar
  # ---------------- #
  # Initialize the vectors
  T.sub <- NULL
  beta.sub <- NULL
  # Initialize the progress bar
  if (progress == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = R, style = 3, width = 20)
    cat("\r")
  } else {
    pb <- NULL
  }

  # ---------------- #
  # Step 2: Conduct the subsampling procedure
  # ---------------- #
  for (i in 1:R){
    # (2.1) Re-sample the data
    data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
    rownames(data.bs) <- 1:nrow(data.bs)

    # (2.2) Compute the bootstrap estimates
    # Compute the value of beta_bs_star using the function func_obs
    sub.return <- subsample.prob(data.bs, lpmodel, beta.tgt, norm, solver, i+1)
    T.sub <- c(T.sub, sub.return$objval)
    beta.sub <- cbind(beta.sub, sub.return$beta)

    # (2.3) Update progress bar
    if (progress == TRUE){
      if (i != R){
        utils::setTxtProgressBar(pb, i)
        cat("\r\r")
      } else {
        utils::setTxtProgressBar(pb, i)
        cat("\r\b")
      }
    }
  }

  # ---------------- #
  # Step 3: Return the results
  # ---------------- #
  return(list(T.sub = T.sub,
              pb = pb))
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
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{T.sub}{List of test statistic from the subsampling procedure.}
#'   \item{pb}{Progress bar object.}
#'
#' @export
#'
subsample.manycores <- function(data, R, lpmodel, beta.tgt, norm, solver,
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
    pb <- utils::txtProgressBar(max=R, style=3, width = 20)

    cat("\r")
    progress <- function(n){
      utils::setTxtProgressBar(pb, n)
      if (n < R) {
        cat("\r\r")
      } else {
        cat("\r\b")
      }
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

  # Use a normal for-loop to construct the list of lpmodel objects
  for (i in 1:R) {
    # Construct bootstrap data
    data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
    rownames(data.bs) <- 1:nrow(data.bs)

    # Assign the lpmodel objects
    lpmodel.bs[[i]] <- list()
    lpmodel.bs[[i]]$A.obs <- lpmodel.eval(data.bs, lpmodel$A.obs, i + 1)
    lpmodel.bs[[i]]$A.shp <- lpmodel.eval(data.bs, lpmodel$A.shp, i + 1)
    lpmodel.bs[[i]]$A.tgt <- lpmodel.eval(data.bs, lpmodel$A.tgt, i + 1)
    lpmodel.bs[[i]]$beta.shp <- lpmodel.eval(data.bs, lpmodel$beta.shp, i + 1)
    beta.obs.return <- lpmodel.beta.eval(data.bs, lpmodel$beta.obs, i + 1)
    lpmodel.bs[[i]]$beta.obs <- beta.obs.return
  }

  # ---------------- #
  # Step 3: Conduct the subsampling procedure
  # ---------------- #
  listans = foreach::foreach(i = 1:R, .multicombine = TRUE,
                             .combine = "comb", .options.snow = opts,
                             .packages = "lpinfer") %dorng% {
   ## (3.1) Compute the bootstrap estimates
   # Compute the value of beta_bs_star using the function func_obs
   sub.return <- subsample.prob(data.bs, lpmodel.bs[[i]], beta.tgt, norm,
                                solver, 1)
   T.sub <- data.frame(sub.return$objval)
   beta.sub <- data.frame(c(sub.return$beta))

   ## (3.2) Combine the results
   list(T.sub, beta.sub)
  }

  # ---------------- #
  # Step 4: Retrieve results from output
  # ---------------- #
  T.sub = as.vector(unlist(listans[[1]]))
  beta.sub = data.frame(matrix(unlist(listans[[2]]),
                               ncol = R,
                               byrow = FALSE))

  # ---------------- #
  # Step 5: Return the results
  # ---------------- #
  return(list(T.sub = T.sub,
              pb = pb))
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
  cat(sprintf("p-value: %s\n", round(x$pval, digits = 5)))
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
  # Print the p-values
  print(x)

  # Print test statistic, solver used, norm used and number of
  # cores used
  cat(sprintf("Test statistic: %s\n", round(x$T.n, digits = 5)))
  cat(sprintf("Solver used: %s\n", x$solver))
  cat(sprintf("Norm used: %s\n", x$norm))
  cat(sprintf("Phi used: %s\n", round(x$phi, digits = 5)))
  cat(sprintf("Size of each subsample: %s\n", x$subsample.size))
  cat(sprintf("Number of cores used: %s\n", x$cores))
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
#'    }
#'
#' @export
#'
subsample.check <- function(data, lpmodel, beta.tgt, R, solver, cores,
                            norm, phi, alpha, replace, progress){

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

  # Check numerics
  check.numrange(phi, "phi", "closed", 0, "closed", 1)
  check.numeric(phi, "phi")
  check.positiveinteger(R, "R")
  cores <- check.cores(cores)

  # Check norm
  norm <- check.norm(norm, "norm")

  # Check Boolean
  check.boolean(replace, "replace")
  check.boolean(progress, "progress")

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              cores = cores,
              norm = norm))
}
