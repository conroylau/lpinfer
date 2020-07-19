#' Computes the \eqn{p}-value of the \code{dkqs} procedure
#'
#' @description This module conducts inference in quadratic programs using
#'    the procedure suggested by Torgovitsky (2019) that incorporates
#'    the cone-tightening procedure proposed by Deb, Kitamura, Quah
#'    and Stoye (2018).
#'
#' @import foreach doMC parallel
#'
#' @param data The data being used in the inference.
#' @param lpmodel A list of objects that are used in inference of linear
#'    programming problems. The list of objects required in the \code{dkqs}
#'    procedure are:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{beta.obs}}
#'    }
#' @param beta.tgt Value of beta to be tested.
#' @param R Number of bootstraps chosen by the users.
#' @param Rmulti Multiplier for the number of bootstrap replications. The
#'   product of `\code{Rmulti}' and `\code{R}' refers to the maximum
#'   number of bootstrap replications.
#' @param tau The value of tau chosen by the user.
#' @param solver The name of the linear and quadratic programming solver that
#'    are used to obtain the solution to linear and quadratic programs.
#'    The solvers supported by this module are `\code{cplexAPI}',
#'    `\code{gurobi}', `\code{limSolve}' and `\code{Rcplex}'.
#' @param cores Number of cores to be used in the parallelized for-loop.
#'    Parallelized for-loop is used if \code{cores} is set to be an integer
#'    greater than or equal to 2.
#' @param progress The boolean variable for whether the result messages
#'    should be displayed in the inference procedure. If it is set as
#'    \code{TRUE}, the messages are displayed throughout the procedure.
#'    Otherwise, the messages will not be displayed.
#'
#' @return Returns a list of output calculated from the function:
#'   \item{pval}{A table of \eqn{p}-values for each \eqn{\tau}.}
#'   \item{tau.feasible}{The list of \eqn{\tau}s that are feasible.}
#'   \item{tau.ineffective}{The list of \eqn{\tau}s that are infeasible.}
#'   \item{tau.max}{Maximum value of the feasible tau for the problem.}
#'   \item{T.n}{Test statistic \eqn{T.n}.}
#'   \item{T.bs}{The list of bootstrap test statistics
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}} for each \eqn{\tau}.}
#'   \item{beta.bs.bar}{The list of \eqn{\tau}-tightened re-centered bootstrap
#'      estimators \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'   \item{lb0}{Logical lower bound of the problem for each \eqn{\tau}.}
#'   \item{ub0}{Logical upper bound of the problem for each \eqn{\tau}.}
#'   \item{solver}{Solver used in solving the linear and quadratic programs.}
#'   \item{cores}{Number of cores used.}
#'   \item{cv.table}{Table of critical values.}
#'   \item{call}{The function that has been called.}
#'   \item{test.logical}{Indicator variable for whether the computation has
#'     been conducted. If '\code{test.logical}' is 1, it refers to the case
#'     where '\code{beta.tgt}' is inside the logical bound. If
#'     '\code{test.logical}' is 0, it refers to the case where '
#'     \code{beta.tgt}' is outside the logical bound.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @details If the value of the test statistic \eqn{T.n} is zero, the
#'    bootstrap procedure will be skipped.
#'
#' @export
#'
dkqs <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                 tau = NULL, solver = NULL, cores = 1, progress = TRUE){

  # ---------------- #
  # Step 1: Update call, check and update the arguments
  # ---------------- #
  # Obtain call information
  call <- match.call()

  # Check the arguments
  dkqs.return <- dkqs.check(data, lpmodel, beta.tgt, R, Rmulti, tau, solver,
                            cores, progress)

  # Update the arguments
  data <- dkqs.return$data
  lpmodel <- dkqs.return$lpmodel
  solver <- dkqs.return$solver
  solver.name <- dkqs.return$solver.name
  cores <- dkqs.return$cores
  test.logical <- dkqs.return$test.logical

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds
  if (test.logical == 1) {
    # ---------------- #
    # Step 2: Initialization
    # ---------------- #
    n <- nrow(data)
    # Choose whether the data frame has "Y" or "y" as the column name for the
    # outcomes
    if ("y" %in% colnames(data)) {
      coly <- "y"
    } else if  ("Y" %in% colnames(data)) {
      coly <- "Y"
    } else {
      stop(paste0("'data' needs to have a column called 'y' or 'Y' to ",
                  "represent the outcomes."))
    }
    J <- length(unique(data[,coly])) - 1

    # Initialize beta.obs
    beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)
    beta.obs.hat <- beta.obs.return[[1]]

    # ---------------- #
    # Step 3: Choose the value of tau
    # ---------------- #
    tau.return <- dkqs.qlp(lpmodel, beta.tgt, beta.obs.hat, 1, "tau",
                           n, solver)
    if (is.null(tau)) {
      tau.feasible <- as.numeric(tau.return$objval)
      tau.infeasible <- NULL
    } else {
      tau.feasible <- c()
      tau.infeasible <- c()
      for (i in 1:length(tau)) {
        if (tau[i] > tau.return$objval) {
          tau.feasible <- c(tau.feasible, tau.return$objval)
          tau.infeasible <- c(tau.infeasible, tau[i])
        } else if (tau[i] <= tau.return$objval) {
          tau.feasible <- c(tau.feasible, tau[i])
        } else {
          # Error message when the problem is infeasible.
          stop("The problem is infeasible. Choose other values of tau.")
        }
      }

      # Remove the duplicates
      tau.feasible <- sort(unique(tau.feasible))
      tau.infeasible <- sort(unique(tau.infeasible))
    }

    n.tau <- length(tau.feasible)

    # ---------------- #
    # Step 4: Compute T.n, x.star and s.star
    # ---------------- #
    # Initialize the data.frames and lists to contain the optimal value of x,
    # bootstrap test statistics and p-values for each tau
    pval.df <- data.frame(matrix(vector(), nrow = n.tau, ncol = 2))
    pval.df[,1] <- tau.feasible
    lb.df <- pval.df
    ub.df <- pval.df
    colnames(pval.df) <- c("tau", "p-value")
    colnames(lb.df) <- c("tau", "lb")
    colnames(ub.df) <- c("tau", "ub")

    x.star.list <- list()
    s.star.list <- list()
    T.bs.list <- list()
    beta.bs.bar.list <- list()

    # Compute T.n (Here, the value of tau does not affect the problem)
    T.n <- dkqs.qlp(lpmodel, beta.tgt, beta.obs.hat, tau.feasible[1], "test",
                    n, solver)$objval

    if (T.n == 0) {
      cat(paste0("Bootstrap is skipped because the ",
                 "value of the test statistic is zero.\n"))
      T.n <- 0
      T.bs.list[[i]] <- NULL
      beta.bs.bar.list[[i]] <- NULL
    } else {
      # Compute s.star, x.star, bootstrap test statistics and p-values for each
      # tau
      for (i in 1:n.tau) {
        # ---------------- #
        # Step 5: Compute x.star and s.star for each tau
        # ---------------- #
        x.return <- dkqs.qlp(lpmodel, beta.tgt, beta.obs.hat, tau.feasible[i],
                             "cone", n, solver)
        x.star.list[[i]] <- x.return$x
        s.star.list[[i]] <- lpmodel$A.obs %*% x.star.list[[i]]

        # ---------------- #
        # Step 6: Obtain logical bounds for the function invertci
        # ---------------- #
        lb.df[i,2] <- x.return$lb0$objval
        ub.df[i,2] <- x.return$ub0$objval
      }

      # ---------------- #
      # Step 7: Compute the bootstrap beta and estimates
      # ---------------- #
      if (cores == 1) {
        # No parallelization
        T.bs.return <- beta.bs(data, lpmodel, beta.tgt, R, maxR, J,
                               s.star.list, tau.feasible, solver, progress)
      } else {
        # Parallelization
        T.bs.return <- beta.bs.parallel(data, lpmodel, beta.tgt, R, maxR, J,
                                        s.star.list, tau.feasible,
                                        solver, progress, cores)
      }

      # ---------------- #
      # Step 8: Compute p-value
      # ---------------- #
      for (i in 1:n.tau) {
        pval.df[i,2] <- pval(T.bs.return$T.bs[[i]], T.n)$p
      }

      # ---------------- #
      # Step 9: Generate a table of critical values
      # ---------------- #
      cv.table <- construct.cv.table(tau.feasible, "tau", rep(T.n, n.tau),
                                     T.bs.return$T.bs)
      # ---------------- #
      # Step 10: Close the progress bar that is used in the bootstrap procedure
      # ---------------- #
      if ((progress == TRUE) & (i == n.tau)){
        close(T.bs.return$pb)
        # Remove progress bar
        cat("                               \n\b\r")
      }

    }

    # ---------------- #
    # Step 11: Assign the return list and return output
    # ---------------- #
    output <- list(pval = pval.df,
                   tau.feasible = tau.feasible,
                   tau.infeasible = tau.infeasible,
                   tau.max = tau.return$objval,
                   T.n = T.n,
                   T.bs = T.bs.return$T.bs,
                   beta.bs.bar = T.bs.return$beta.bs.bar.list,
                   lb0 = lb.df,
                   ub0 = ub.df,
                   solver = solver.name,
                   cores = cores,
                   cv.table = cv.table,
                   call = call,
                   test.logical = test.logical,
                   df.error = T.bs.return$df.error,
                   R.succ = T.bs.return$R.succ)
  } else {
    ### Case 2: test.logical == 0. Set the p-value as 0 directly because
    ### beta.tgt is outside the logical bounds
    output <- list(pval = 0,
                   solver = solver.name,
                   cores = cores,
                   call = call,
                   test.logical = test.logical)

    # Print warning message
    infeasible.betatgt.warning()

  }
  attr(output, "class") <- "dkqs"

  # Return output
  return(output)
}

#' Formulating and solving linear and quadratic programs
#'
#' @description This function formulates the matrices and vectors, and
#'    solves the quadratic programs (4) or linear programs (5) and (6) in
#'    Torgovitsky (2019).
#'
#' @param tau The value of tau to be used in the linear program.
#' @param problem The problem that the function will be solved.
#' @param solver Name of the solver that solves the linear and quadratic
#'    programs.
#' @param beta.obs.hat The value of \eqn{\hat{\beta}_{\mathrm{obs}}} based on
#'    that is either the value of the observed.
#' @param n The number of observations in the dataframe.
#' @inheritParams dkqs
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @details The argument \code{problem} must be one of the followings:
#' \itemize{
#'   \item{\code{test} --- this computes the quadratic program for the test
#'     statistic, i.e. quadratic program (4) in Torgovitsky (2019)}
#'   \item{\code{cone} --- this computes the quadratic program for the
#'     bootstrap test statistics, i.e. quadratic program (5) in
#'     Torgovitsky (2019)}
#'   \item{\code{tau} --- this computes the value of tau based on the
#'     procedure suggested by Kamat (2018), i.e. linear program (6) in
#'     Torgovitsky (2019)}
#' }
#'
#' @export
#'
dkqs.qlp <- function(lpmodel, beta.tgt, beta.obs.hat, tau, problem, n,
                     solver){
  # ---------------- #
  # Step 1: Obtain the dimension of the
  # ---------------- #
  A.tgt.dim <- dim(lpmodel$A.tgt)
  if (is.null(A.tgt.dim)) {
    A.tgt.nr <- 1
    A.tgt.nc <- length(lpmodel$A.tgt)
  } else {
    A.tgt.nr <- A.tgt.dim[1]
    A.tgt.nc <- A.tgt.dim[2]
  }

  # ---------------- #
  # Step 2: Formulation of constraints
  # ---------------- #

  ones <- matrix(rep(1, A.tgt.nc), nrow = 1)
  lb <- matrix(rep(0, A.tgt.nc), nrow = 1)

  # Theta parameters
  theta.down <- do.call(solver, list(Af  = NULL,
                                     bf  = lpmodel$A.tgt,
                                     nf  = n,
                                     A   = ones,
                                     rhs = c(1),
                                     sense = "=",
                                     modelsense = "min",
                                     lb = lb))
  theta.up <- do.call(solver, list(Af  = NULL,
                                   bf  = lpmodel$A.tgt,
                                   nf  = n,
                                   A   = ones,
                                   rhs = c(1),
                                   sense = "=",
                                   modelsense ="max",
                                   lb = lb))

  # Obtain required set of indices
  x.ind <- 1:A.tgt.nc
  ind.up <- which(lpmodel$A.tgt %in% theta.up$objval)
  ind.down <- which(lpmodel$A.tgt %in% theta.down$objval)
  ind.0 <- x.ind[-c(ind.up, ind.down)]

  # Updated lb for certain indices
  rhs.up <- (beta.tgt - theta.down$objval) * tau / length(c(ind.0, ind.up))
  rhs.down <- (theta.up$objval - beta.tgt) * tau / length(c(ind.0, ind.down))
  rhs.0 <- (1 - rhs.up * length(ind.up) -
              rhs.down * length(ind.down)) * tau / length(ind.0)

  # RHS vector and sense for the linear or quadratic program
  lp.rhs <- c(beta.tgt, 1)
  lp.sense <- c("=", "=")

  # ---------------- #
  # Step 3: Solve the QP
  # - If problem == "test", this function solves LP (4)
  # - If problem == "cone", this function solves LP (5)
  # - If problem == "tau", this function solves LP (6)
  # ---------------- #
  if (problem == "test"){
    ans <- do.call(solver, list(Af  = lpmodel$A.obs,
                                bf  = beta.obs.hat,
                                nf  = n,
                                A   = rbind(lpmodel$A.tgt, ones),
                                rhs = lp.rhs,
                                sense = lp.sense,
                                modelsense ="min",
                                lb = lb))
  } else if (problem == "cone"){
    # Update lb
    lb.new <- lb
    lb.new[ind.up] <- rhs.up
    lb.new[ind.down] <- rhs.down
    lb.new[ind.0] <- rhs.0

    # Find the solution using the solver
    ans <- do.call(solver, list(Af  = lpmodel$A.obs,
                                bf  = beta.obs.hat,
                                nf  = n,
                                A   = rbind(lpmodel$A.tgt, ones),
                                rhs = lp.rhs,
                                sense = lp.sense,
                                modelsense = "min",
                                lb = lb.new))
  } else if (problem == "tau") {
    lp.lhs.tau <- rbind(lpmodel$A.tgt, ones)
    lp.lhs.tau <- cbind(matrix(c(0,0), nrow = 2), lp.lhs.tau)

    # Add one unit because of the additional position for tau
    len.tau <- A.tgt.nc + 1
    lb.tau <- rep(0, len.tau)
    lp.rhs.tau <- lp.rhs
    lp.sense.tau <- lp.sense
    # Inequality constraints for ind.up
    if (length(ind.up) != 0) {
      for (i in 1:length(ind.up)) {
        new.const <- tau_constraints(len.tau, rhs.up, -1, ind.up[i] + 1, 0,
                                     "<=", lp.lhs.tau, lp.rhs.tau,
                                     lp.sense.tau)
        lp.lhs.tau <- new.const$lp.lhs.tau
        lp.rhs.tau <- new.const$lp.rhs.tau
        lp.sense.tau <- new.const$lp.sense.tau
      }
    }
    # Inequality constraints for ind.down
    if (length(ind.down) != 0) {
      for (i in 1:length(ind.down)) {
        new.const <- tau_constraints(len.tau, rhs.down, -1, ind.down[i] + 1, 0,
                                     "<=", lp.lhs.tau, lp.rhs.tau,
                                     lp.sense.tau)
        lp.lhs.tau <- new.const$lp.lhs.tau
        lp.rhs.tau <- new.const$lp.rhs.tau
        lp.sense.tau <- new.const$lp.sense.tau
      }
    }
    if (length(ind.0) != 0) {
      # Inequality constraints for ind.0
      for (i in 1:length(ind.0)) {
        new.const <- tau_constraints(len.tau, rhs.0, -1, ind.0[i] + 1, 0, "<=",
                                     lp.lhs.tau, lp.rhs.tau, lp.sense.tau)
        lp.lhs.tau <- new.const$lp.lhs.tau
        lp.rhs.tau <- new.const$lp.rhs.tau
        lp.sense.tau <- new.const$lp.sense.tau
      }
    }
    ans <- do.call(solver, list(Af  = NULL,
                                bf  = c(1, rep(0, A.tgt.nc)),
                                nf  = n,
                                A   = lp.lhs.tau,
                                rhs = lp.rhs.tau,
                                sense = lp.sense.tau,
                                modelsense = "max",
                                lb = lb.tau))
  }

  # ---------------- #
  # Step 4: Append the logical bounds to the results
  # ---------------- #
  ans$lb0 <- theta.down
  ans$ub0 <- theta.up
  return(ans)
}

#' Computes the bootstrap test statistics
#'
#' @description This function computes the bootstrap test statistics.
#'
#' @import modelr dplyr
#'
#' @param J The number of distinct nonzero values in vector \eqn{\bm{y}}.
#' @param s.star.list The list of values of
#'    \eqn{\hat{s}^\ast \equiv A_{\mathrm{obs}}\hat{\bm{x}}_n^\ast}
#'    in the cone-tightening procedure for each tau.
#' @param tau.list The list of feasible taus.
#' @param maxR Maximum number of bootstrap replications in case error occured.
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the list of estimates from bootstrap:
#'   \item{T.bs}{A list of bootstrap test statistics
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'  \item{beta.bs.bar.list}{A list of \eqn{\tau_n}-tightened recentered
#'     bootstrap estimates \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'  \item{pb}{Progress bar object.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
beta.bs <- function(data, lpmodel, beta.tgt, R, maxR, J, s.star.list, tau.list,
                    solver, progress){
  # ---------------- #
  # Step 1: Initialize the vectors and progress counters
  # ---------------- #
  T.bs <- NULL
  beta.bs.bar.list <- NULL

  # Initialize the progress bar
  if (progress == TRUE) {
    pb <- utils::txtProgressBar(initial = 0, max = R, style = 3, width = 20)
    cat("\r")
  } else {
    pb <- NULL
  }

  # Count the total number of taus
  n.tau <- length(tau.list)

  # Initialize the lists to store the test statistics
  T.bs <- list()
  beta.bs.bar.list <- list()

  # Initialize a table to contain the error messages
  df.error <- data.frame(matrix(vector(), ncol = 3))
  colnames(df.error) <- c("Iteration", "tau", "Error message")

  # ---------------- #
  # Step 2: Bootstrap procedure
  # ---------------- #
  lpmodel.bs <- lpmodel

  # Use the for-loop to construct the bootstrap test statistic
  # There are two error-handling parts here. If there is an error from
  # constructing beta.obs from the bootstrapped data, go to the next bootstrap
  # replication.
  # If not, construct the bootstrap estimators for each tau. If there is an
  # error in the procedure for one of the taus. Go to the next bootstrap
  # replication.
  for (i in 1:maxR) {
    data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
    rownames(data.bs) <- 1:nrow(data.bs)

    ## (2.1) Evaluate beta.obs from bootstrap data
    beta.obs.result <- tryCatch(
      expr = {
        # Compute the bootstrap test statistic
        if (class(lpmodel$beta.obs) == "function"){
          beta.obs.bs <- lpmodel.beta.eval(data.bs, lpmodel$beta.obs, 1)[[1]]
          beta.obs <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
        } else if (class(lpmodel$beta.obs) == "list") {
          beta.obs.bs <- lpmodel.beta.eval(data, lpmodel$beta.obs, i + 1)[[1]]
          beta.obs <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
        }
        beta.obs.ls <- list(status = "NOERROR",
                            beta.obs = beta.obs,
                            beta.obs.bs = beta.obs.bs)
      },
      error = function(e) {
        return(list(status = "ERROR",
                    msg = e))
      },
      finally = {
        beta.obs.ls
      }
    )

    ## (2.2) Check if there is any error in forming beta.obs
    if (beta.obs.result$status == "NOERROR") {
      # If there is no error, start to construct the estimators for each tau
      for (j in 1:n.tau) {
        # Compute beta.bs.bar and test statistic
        beta.bs.bar <- beta.obs.bs - beta.obs + s.star.list[[j]]

        # Check if there is any error in the program
        bstau.result <- tryCatch(
          expr = {
            T.bs.j <- dkqs.qlp(lpmodel.bs, beta.tgt, beta.bs.bar, tau.list[j],
                               "cone", nrow(data), solver)$objval
            T.bs.ls <- list(status = "NOERROR",
                            T.bs.j = T.bs.j)
          },
          error = function(e) {
            return(list(status = "ERROR",
                        msg = e))
          },
          finally = {
            T.bs.ls
          }
        )

        # Append the results if no error
        if (bstau.result$status == "NOERROR") {
          if (i == 1) {
            T.bs[[j]] <- T.bs.j
            beta.bs.bar.list[[j]] <- beta.bs.bar
          } else {
            T.bs[[j]] <- c(T.bs[[j]], T.bs.j)
            beta.bs.bar.list[[j]] <- cbind(beta.bs.bar.list[[j]], beta.bs.bar)
          }
        } else {
          # Break the loop if there is an error in one of the taus
          break()
        }
      }
    }

    # (2.3) Save the index, tau, and the error message (if any)
    if (beta.obs.result$status == "ERROR" | bstau.result$status == "ERROR") {
      df.error[nrow(df.error) + 1, 1] <- i
      if (beta.obs.result$status != "ERROR") {
        df.error[nrow(df.error) + 1, 2] <- tau.list[j]
        df.error[nrow(df.error), 3] <- bstau.result$msg$message
      } else {

        df.error[nrow(df.error) + 1, 2] <- NA
        df.error[nrow(df.error), 3] <- beta.obs.result$msg$message
      }
    }

    # (2.4) Update progress bar
    if (progress == TRUE) {
      if (i < maxR) {
        utils::setTxtProgressBar(pb, i)
        cat("\r\r")
      } else if ((i == maxR)) {
        utils::setTxtProgressBar(pb, i)
        cat("\r\b")
      }
    }

    # (2.5) Break the loop if R successful replications are made
    if (length(T.bs[[1]]) == R) {
      if (progress == TRUE) {
        utils::setTxtProgressBar(pb, maxR)
        cat("\r\b")
      }
      break()
    }
  }

  # (2.6) Compute the number of successful bootstrap replications
  R.succ <- length(T.bs[[1]])

  # ---------------- #
  # Step 3: Return results
  # ---------------- #
  return(list(T.bs = T.bs,
              beta.bs.bar.list = beta.bs.bar.list,
              pb = pb,
              R.succ = R.succ,
              df.error = df.error))
}

#' Computes the bootstrap test statistics with parallelization
#'
#' @description This function computes the bootstrap test statistics with
#'    parallelized for-loops.
#'
#' @import doMC foreach doRNG
#'
#' @inheritParams beta.bs
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the list of estimates from bootstrap:
#'   \item{T.bs}{A list of bootstrap test statistics
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'  \item{beta.bs.bar.set}{A list of \eqn{\tau_n}-tightened recentered
#'     bootstrap estimates \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'  \item{pb}{Progress bar object.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
beta.bs.parallel <- function(data, lpmodel, beta.tgt, R, maxR, J, s.star.list,
                             tau.list, solver, progress, cores){
  # ---------------- #
  # Step 1: Register the number of cores and extract information
  # ---------------- #
  options(warn=-1)

  # Register core
  doMC::registerDoMC(cores)

  # Compute dimension
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.nr <- length(lpmodel$beta.obs(data))
  } else if (class(lpmodel$beta.obs) == "list"){
    beta.obs.nr <- lpmodel$beta.obs[[1]]
  }

  # Count the total number of taus
  n.tau <- length(tau.list)

  # ---------------- #
  # Step 2: Initialize progress bar, comb function and assign doRNG
  # ---------------- #
  if (progress == TRUE){
    # Initialize the counter
    cl <- PtProcess::makeSOCKcluster(8)
    doSNOW::registerDoSNOW(cl)

    # Set the counter and progress bar
    pb <- utils::txtProgressBar(initial = 0,
                                max = maxR, style = 3, width = 20)
    cat("\r")
    progress <- function(n){
      utils::setTxtProgressBar(pb, n)
      if (n < maxR){
        cat("\r\r")
      } else if (n == maxR) {
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

  # Assign doRnG
  `%dorng%` <- doRNG::`%dorng%`

  # ---------------- #
  # Step 3: Bootstrap procedure
  # ---------------- #
  # Initialize the parameters and data frame
  k <- 0
  df.error1 <- data.frame(matrix(vector(), ncol = 3))
  colnames(df.error1) <- c("Iteration", "tau", "Error message")
  error.21 <- NULL
  error.22 <- NULL
  error.23 <- NULL

  # Set the estimator of beta.obs
  beta.bs.list <- list()
  beta.bs.list[[1]] <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]

  # Initialize list for final results
  T.bs <- list()
  beta.bs.bar.list <- list()

  # Loop until the number of bootstrap replications match R or if maxR has been
  # reached
  while (k != R) {
    # Denote the starting index and ending index
    i0 <- k + 1
    i1 <- min(i0 + (R - k) - 1, maxR)

    # Use a normal for-loop to construct the list of lpmodel objects
    for (i in i0:i1) {

      # Compute the bootstrapped beta.obs objects
      beta.bs.result <- tryCatch(
        expr = {
          if (class(lpmodel$beta.obs) == "function") {
            data.bs <- as.data.frame(data[sample(1:nrow(data),
                                                 replace = TRUE),])
            rownames(data.bs) <- 1:nrow(data.bs)
            beta.obs.return <- lpmodel.beta.eval(data.bs, lpmodel$beta.obs, 1)
            beta.bs.list[[i + 1]] <- beta.obs.return[[1]]
          } else if (class(lpmodel$beta.obs) == "list") {
            beta.bs.list[[i + 1]] <- lpmodel.beta.eval(data,
                                                       lpmodel$beta.obs,
                                                       i + 1)[[1]]
          }
          beta.bs.ls <- list(status = "NOERROR",
                             lpmodel.bs = beta.bs.list[[i + 1]])
        },
        error = function(e) {
          return(list(status = "ERROR",
                      msg = e))
        },
        finally = {
          beta.bs.ls
        }
      )

      # Record error (if any)
      if (beta.bs.result$status == "ERROR") {
        df.error1[nrow(df.error1) + 1, 1] <- i
        df.error1[nrow(df.error1), 2] <- NA
        df.error1[nrow(df.error1), 3] <- beta.bs.result$msg$message
      }
    }

    # Bootstrap procedure of DKQS
    listans <- foreach::foreach(i = i0:i1, .multicombine = TRUE,
                                .combine = "comb", .options.snow = opts,
                                .packages = c("lpinfer", "doRNG")) %dorng%
      {
        # Initialize the parameters that contain the results for each tau
        T.bs.list <- NULL
        beta.bs.bar.list <- list()

        # Only consider the subsample problem if there is no error in forming
        # the beta.obs parts
        if (!(i %in% df.error1[,1])) {
          ## (3.1) Compute the bootstrap estimates
          for (j in 1:n.tau) {
            result <- tryCatch(
              expr = {
                beta.bs.bar <- beta.bs.list[[i + 1]] - beta.bs.list[[1]] +
                  s.star.list[[j]]

                T.bs <- dkqs.qlp(lpmodel, beta.tgt,  beta.bs.bar, tau.list[j],
                                  "cone", nrow(data), solver)$objval

                T.bs.list <- c(T.bs.list, T.bs)
                beta.bs.bar.list[[j]] <- beta.bs.bar
                beta.bs.ls <- list(status = "NOERROR")
              },
              error = function(e) {
                return(list(status = "ERROR",
                            msg = e))
              },
              finally = {
              }
            )

            # Break the loop if there is an error in one of the taus
            if (result$status %in% c("ERROR")) {
              break()
            }
        }

        ## (3.2) Store the results or error message depending on the status
        if (result$status %in% c("ERROR")) {
          ind <- i
          ind.tau <- j
          ind.msg <- result$msg$message
          T.bs.list <- NULL
          beta.bs.bar.list <- NULL
        } else {
          ind <- NULL
          ind.tau <- NULL
          ind.msg <- NULL
          T.bs.list <- T.bs.list
          beta.bs.bar.list <- beta.bs.bar.list
        }
      } else {
        ind <- NULL
        ind.tau <- NULL
        ind.msg <- NULL
        T.bs.list <- NULL
        beta.bs.bar.list <- NULL
      }

      list(T.bs.list, beta.bs.bar.list, ind, ind.tau, ind.msg)
    }

    ## (3.3) Consolidate the results
    T.bs.temp <- list()
    beta.bs.bar.list.temp <- list()
    for (j in 1:n.tau) {
      # Consolidate the test statistics
      T.bs.temp[[j]] <- unlist(listans[[1]])[seq(j, n.tau*R, n.tau)]

      # Consolidate the betas
      beta.bs.bar.list.temp[[j]] <- listans[[2]][[j]]
      for (i in 2:R) {
        beta.bs.bar.list.temp[[j]] <- cbind(beta.bs.bar.list.temp[[j]],
                                            listans[[2]][[n.tau + i - 1]][[j]])
      }

      # Consolidate the test statistics
      if (i0 == 1) {
        T.bs[[j]] <- T.bs.temp[[j]]
        beta.bs.bar.list[[j]] <- beta.bs.bar.list.temp[[j]]
      } else {
        T.bs[[j]] <- c(T.bs[[j]], T.bs.temp[[j]])
        beta.bs.bar.list[[j]] <- cbind(beta.bs.bar.list[[j]],
                                       beta.bs.bar.list.temp[[j]])
      }
    }


    k <- length(T.bs[[1]])

    ## (3.4) Consolidate the list of error messages
    if (length(unlist(listans[[3]])) != 0) {
      error.21.temp <- unlist(listans[[3]])
      error.22.temp <- unlist(listans[[4]])
      error.23.temp <- unlist(listans[[5]])
      error.21 <- c(error.21, error.21.temp)
      error.22 <- c(error.22, error.22.temp)
      error.23 <- c(error.23, error.23.temp)
    }

    ## (3.5) Break the while-loop if it reached maxR
    if (i1 == maxR) {
      break()
    }

  }
  # Close the progress bar
  cat("\r\b")

  # Number of successful bootstrap replications
  R.succ <- length(T.bs[[1]])

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
  # Step 5: Return results
  # ---------------- #
  return(list(T.bs = T.bs,
              beta.bs.bar.list = beta.bs.bar.list,
              pb = pb,
              df.error = df.error,
              R.succ = R.succ))
}

#' Auxiliary function to calculate the p-value
#'
#' @description This function computes the \eqn{p}-value of the test based on
#'    the bootstrap estimates.
#'
#' @param T_bs The test statistics obtained from bootstrap.
#' @param T.n The test statistics obtained from quadratic program (5).
#' @param alpha The significance level.
#'
#' @return Returns the \eqn{p}-value:
#'   \item{p}{\eqn{p}-value.}
#'   \item{decision}{Decision to reject or not.}
#'
#' @export
#'
pval <- function(T.bs, T.n, alpha = .05) {
  # Compute p-value
  p <- mean(T.bs > T.n)

  # Update decision
  if (p > alpha){
    decision <- 1
  } else {
    decision <- 0
  }
  return(list(p = p,
              decision = decision))
}

#' Auxiliary function to create the constraints for the linear program for tau
#'
#' @description This function generates the constraints in the linear program
#'   for computing the value of tau based on the procedure suggested by Kamat
#'   (2018), i.e. linear program (6) of Torgovitsky (2019).
#'
#' @param length.tau The number of variables in the constraint.
#' @param coeff.tau The coefficient in front of tau in the constraint.
#' @param coeff.x The coefficient in front of \eqn{x_i} in the constraint.
#' @param ind.x The index of \eqn{x_i}, i.e. the value of \eqn{i}.
#' @param rhs The RHS vector of the new constraint.
#' @param sense The equality or inequality symbol to be used in the new
#'    constraint.
#' @param lp.lhs.tau The constraint matrix to be updated.
#' @param lp.rhs.tau The RHS vector of the linear constraints to be updated.
#' @param lp.sense.tau The sense vector fo the linear constraints to be
#'    updated.
#'
#' @return Returns the list of matrices that corresponds to the updated
#'   constraints:
#'   \item{lp.lhs.tau}{Upated constraint matrix.}
#'   \item{lp.rhs.tau}{Update RHS vector.}
#'   \item{lp.sense.tau}{Update sense for the constraints.}
#'
#' @export
#'
tau_constraints <- function(length.tau, coeff.tau, coeff.x, ind.x, rhs, sense,
                            lp.lhs.tau, lp.rhs.tau, lp.sense.tau){
  # ---------------- #
  # Step 1: Create and append the new constraints
  # ---------------- #
  temp <- rep(0, length.tau)
  temp[1] <- coeff.tau
  temp[ind.x] <- coeff.x

  # Update the lhs, rhs and sense of the constraints
  lp.lhs.tau <- rbind(lp.lhs.tau, c(temp))
  lp.rhs.tau <- c(lp.rhs.tau, 0)
  lp.sense.tau <- c(lp.sense.tau, sense)

  # ---------------- #
  # Step 2: Return the list of updated constraints for the tau-problem
  # ---------------- #
  return(list(lp.lhs.tau = lp.lhs.tau,
              lp.rhs.tau = lp.rhs.tau,
              lp.sense.tau = lp.sense.tau))
}

#' Checks and updates the input
#'
#' @description This function checks and updates the input of the user. If
#'    there is any invalid input, this function will be terminated and
#'    generates appropriate error messages.
#'
#' @inheritParams dkqs
#'
#' @return Returns the list of updated parameters as follows:
#'   \item{data}{Upated data in class \code{data.frame}}
#'   \item{lpmodel}{A list of linear programming objects in this
#'      `\code{lpinfer}` package.}
#'   \item{solver}{Updated name of solver in lower case.}
#'   \item{cores}{Updated number of cores to be used in the parallelization
#'      of the for-loop in the bootstrap procedure.}
#'
#' @export
#'
dkqs.check <- function(data, lpmodel, beta.tgt, R, Rmulti, tau, solver, cores,
                       progress){
  # ---------------- #
  # Step 1: Check the arguments
  # ---------------- #
  # Check data
  data <- check.dataframe(data)

  # Check lpmodel
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = "matrix",
                           A.obs.cat = "matrix",
                           A.shp.cat = "not_used",
                           beta.obs.cat = c("function_mat",
                                            "list",
                                            "function_obs_var"),
                           beta.shp.cat = "not_used",
                           R = R)

  # Check solver
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # Check Rmulti
  check.numrange(Rmulti, "Rmulti", "closed", 1, "open", Inf)

  # Check numerics
  check.numeric(beta.tgt, "beta.tgt")
  check.positiveinteger(R, "R")
  for (i in 1:length(tau)) {
    check.numeric(tau[i], "tau")
  }
  cores <- check.cores(cores)

  # Check Boolean
  check.boolean(progress, "progress")

  # Check whether beta.tgt is within the logical bounds
  # Create temporary lpmodel to incorporate the shape constraints in the DKQS
  # procedure
  lpmodel.temp <- lpmodel
  lpmodel.temp$A.shp <- rep(1, ncol(lpmodel$A.obs))
  lpmodel.temp$beta.shp <- 1
  test.logical <- check.betatgt(data, lpmodel.temp, beta.tgt, solver)

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              cores = cores,
              test.logical = test.logical))
}

#' Print results from \code{dkqs}
#'
#' @description This function uses the print method on the return list of the
#'    function \code{dkqs}.
#'
#' @param x Object returned from \code{dkqs}.
#' @param ... Additional arguments.
#'
#' @details The following information are printed:
#'  \itemize{
#'     \item{Test statistic}
#'     \item{\eqn{p}-value}
#'     \item{\eqn{\tau}}
#'     \item{Solver used}
#'     \item{Number of cores used}
#'  }
#'
#' @return Print the basic set of results from \code{dkqs}.
#'
#' @export
#'
print.dkqs <- function(x, ...){
  cat("\r\r")

  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-values
    df.pval <- x$pval

    # Merge the table with the infeasible taus
    if (!is.null(x$tau.infeasible)) {
      df.infeasible <- cbind(x$tau.infeasible, "infeasible")
      colnames(df.infeasible) <- colnames(df.pval)
      df.pval <- rbind(df.pval, df.infeasible)
    }

    # Print the p-values
    if (nrow(df.pval) == 1) {
      cat(sprintf(" p-value: %s\n", df.pval[1,2]))
    } else {
      print(df.pval, row.names = FALSE)
    }
  } else {
    # Case 2: 'beta.tgt' is outside the logical bound
    infeasible.pval.msg()
  }
}

#' Summary of results from \code{dkqs}
#'
#' @description This function uses the summary method on the return list of the
#'    function \code{dkqs}. This is a wrapper of the \code{print} command.
#'
#' @param x Object returned from \code{dkqs}.
#' @param ... Additional arguments.
#'
#' @return Print the summary of the basic set of results from \code{dkqs}.
#'
#' @export
#'
summary.dkqs <- function(x, ...){

  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-values
    print(x)

    # Print maximum feasible tau, test statistic, solver used, number of
    # cores used and the number of bootstrap replications
    cat(sprintf(" Maximum feasible tau: %s\n", round(x$tau.max, digits = 5)))
    cat(sprintf(" Test statistic: %s\n", round(x$T.n, digits = 5)))
    cat(sprintf(" Solver used: %s\n", x$solver))
    cat(sprintf(" Number of cores used: %s\n", x$cores))
    cat(sprintf(" Number of successful bootstrap replications: %s\n",
                x$R.succ))

    # Number of failed bootstrap replications
    if (!is.null(x$df.error) & nrow(x$df.error) != 0) {
      nerr <- nrow(x$df.error)
      errstring <- " Number of failed bootstrap"
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
