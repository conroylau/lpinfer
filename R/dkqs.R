#' Conducts inference using the DKQS procedure
#'
#' @description This module conducts inference using the cone-tightening
#'   procedure proposed by Deb, Kitamura, Quah and Stoye (2018).
#'
#' @param data Data used in the tests.
#' @param lpmodel The \code{lpmodel} object used in the test. The following
#'   components are required in the \code{lpmodel} for the DKQS test:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{beta.obs}}
#'    }
#' @param beta.tgt Value of beta to be tested.
#' @param R Number of bootstrap replications.
#' @param Rmulti Multiplier for the number of bootstrap replications. The
#'   product of \code{Rmulti} and \code{R} refers to the maximum
#'   number of bootstrap replications to be conducted if there are errors.
#' @param tau Value of the tuning parameter \eqn{\tau} in the DKQS
#'   procedure.
#' @param solver Name of the linear and quadratic programming solver that
#'    are used to obtain the solution to linear and quadratic programs.
#'    The solvers supported by this package are \code{cplexAPI}, \code{gurobi},
#'    \code{limSolve} and \code{Rcplex}.
#' @param progress The boolean variable for whether the progress bars should
#'    be displayed. If it is set as \code{TRUE}, the progress bars will be
#'    displayed while the code is running.
#' @param n Sample size (only required if \code{data} is omitted in the input).
#' @return Returns the following list of outputs:
#'   \item{pval}{A table of \eqn{p}-values for each \eqn{\tau}.}
#'   \item{tau.feasible}{The list of \eqn{\tau} that are feasible.}
#'   \item{tau.infeasible}{The list of \eqn{\tau} that are infeasible.}
#'   \item{tau.max}{Maximum value of the feasible \eqn{\tau} for the problem.}
#'   \item{T.n}{Test statistic \eqn{T_n}.}
#'   \item{T.bs}{The list of bootstrap test statistics
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}} for each \eqn{\tau}.}
#'   \item{beta.bs.bar}{The list of \eqn{\tau}-tightened re-centered bootstrap
#'      estimators \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'   \item{lb0}{Logical lower bound of the problem for each \eqn{\tau}.}
#'   \item{ub0}{Logical upper bound of the problem for each \eqn{\tau}.}
#'   \item{solver}{Solver used.}
#'   \item{cv.table}{Table of critical values.}
#'   \item{call}{The function that has been called.}
#'   \item{test.logical}{Indicator variable for whether the computation has
#'     been conducted. If \code{test.logical} is 1, it refers to the case
#'     where \code{beta.tgt} is inside the logical bounds. If
#'     \code{test.logical} is 0, it refers to the case where
#'     \code{beta.tgt} is outside the logical bounds.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @details If the value of the test statistic \eqn{T_n} is zero, the
#'    bootstrap procedure will be skipped and the \eqn{p}-value is zero.
#'
#' @export
#'
dkqs <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                 tau = NULL, n = NULL, solver = NULL, progress = TRUE) {
  # ---------------- #
  # Step 1: Update call, check and update the arguments
  # ---------------- #
  # Obtain call information
  call <- match.call()

  # Check the arguments
  dkqs.return <- dkqs.check(data, lpmodel, beta.tgt, R, Rmulti, tau, n, solver,
                            progress)

  # Update the arguments
  data <- dkqs.return$data
  lpmodel <- dkqs.return$lpmodel
  solver <- dkqs.return$solver
  solver.name <- dkqs.return$solver.name
  test.logical <- dkqs.return$test.logical

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds
  if (test.logical == 1) {
    # ---------------- #
    # Step 2: Initialization
    # ---------------- #
    if (!is.null(data)) {
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
    }

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
    pval.df[, 1] <- tau.feasible
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
      # Compute s.star, x.star, their bootstrap estimates and p-values for each
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
        lb.df[i, 2] <- x.return$lb0$objval
        ub.df[i, 2] <- x.return$ub0$objval
      }

      # ---------------- #
      # Step 7: Compute the bootstrap beta and estimates
      # ---------------- #
      T.bs.return <- dkqs.bs(data, lpmodel, beta.tgt, R, maxR, s.star.list,
                             tau.feasible, solver, progress, n)
      R.succ <- T.bs.return$R.succ

      if (R.succ != 0) {
        # ---------------- #
        # Step 8: Compute p-value
        # ---------------- #
        for (i in 1:n.tau) {
          pval.df[i, 2] <- pval(T.bs.return$T.bs[[i]], T.n)$p
        }

        # ---------------- #
        # Step 9: Generate a table of critical values
        # ---------------- #
        cv.table <- construct.cv.table(tau.feasible, "tau", rep(T.n, n.tau),
                                       T.bs.return$T.bs)
      } else {
        pval.df[,2] <- NA
        cv.table <- NULL
      }
    }

    # ---------------- #
    # Step 10: Assign the return list and return output
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
                   call = call,
                   test.logical = test.logical)

    # Print warning message
    infeasible.betatgt.warning()
  }
  attr(output, "class") <- "dkqs"

  # Return output
  return(output)
}

#' Formulates and solves the linear and quadratic programs in the \code{dkqs}
#' procedure
#'
#' @description This function formulates the matrices and vectors used in
#'    the \code{dkqs} test and solves the quadratic programs (4) or linear
#'    programs (5) and (6) in Torgovitsky (2019).
#'
#' @param tau The value of the tuning parameter tau to be used in the
#'    \code{dkqs} procedure.
#' @param problem The problem that the function will be solved.
#' @param beta.obs.hat The value of sample \eqn{\hat{\bm{\beta}}_{\mathrm{obs}}}
#'    from the \code{lpmodel} object.
#' @param n Number of observations for the data.
#' @inheritParams dkqs
#'
#' @return Returns the optimal point and optimal value.
#'  \item{objval}{Optimal objective value.}
#'  \item{x}{Optimal point.}
#'
#' @details The argument \code{problem} must be one of the followings:
#' \itemize{
#'   \item{\code{test} --- this computes the solution to the quadratic program
#'     that solves the test statistic, i.e. quadratic program (4) in
#'     Torgovitsky (2019)}
#'   \item{\code{cone} --- this computes the solution to the quadratic program
#'     for the bootstrap test statistics, i.e. quadratic program (5) in
#'     Torgovitsky (2019)}
#'   \item{\code{tau} --- this computes the value of tau based on the
#'     procedure suggested by Kamat (2018), i.e. linear program (6) in
#'     Torgovitsky (2019)}
#' }
#'
#' @export
#'
dkqs.qlp <- function(lpmodel, beta.tgt, beta.obs.hat, tau, problem, n,
                     solver) {
  # ---------------- #
  # Step 1: Obtain the dimension of the A matrix
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
                                     nf  = 1,
                                     A   = ones,
                                     rhs = c(1),
                                     sense = "=",
                                     modelsense = "min",
                                     lb = lb))
  theta.up <- do.call(solver, list(Af  = NULL,
                                   bf  = lpmodel$A.tgt,
                                   nf  = 1,
                                   A   = ones,
                                   rhs = c(1),
                                   sense = "=",
                                   modelsense = "max",
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
  if (problem == "test") {
    ans <- do.call(solver, list(Af  = lpmodel$A.obs,
                                bf  = beta.obs.hat,
                                nf  = n,
                                A   = rbind(lpmodel$A.tgt, ones),
                                rhs = lp.rhs,
                                sense = lp.sense,
                                modelsense = "min",
                                lb = lb))
  } else if (problem == "cone") {
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
    lp.lhs.tau <- cbind(matrix(c(0, 0), nrow = 2), lp.lhs.tau)

    # Add one unit because of the additional position for tau
    len.tau <- A.tgt.nc + 1
    lb.tau <- rep(0, len.tau)
    lp.rhs.tau <- lp.rhs
    lp.sense.tau <- lp.sense
    # Inequality constraints for ind.up
    if (length(ind.up) != 0) {
      for (i in 1:length(ind.up)) {
        new.const <- tau.constraints(len.tau, rhs.up, -1, ind.up[i] + 1, 0,
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
        new.const <- tau.constraints(len.tau, rhs.down, -1, ind.down[i] + 1, 0,
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
        new.const <- tau.constraints(len.tau, rhs.0, -1, ind.0[i] + 1, 0, "<=",
                                     lp.lhs.tau, lp.rhs.tau, lp.sense.tau)
        lp.lhs.tau <- new.const$lp.lhs.tau
        lp.rhs.tau <- new.const$lp.rhs.tau
        lp.sense.tau <- new.const$lp.sense.tau
      }
    }
    ans <- do.call(solver, list(Af  = NULL,
                                bf  = c(1, rep(0, A.tgt.nc)),
                                nf  = 1,
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

#' Bootstrap procedure for the DKQS test
#'
#' @description This function carries out the bootstrap procedure of the
#'   DKQS test. This function supports parallel programming via the
#'   \code{future.apply} package.
#'
#' @import future.apply progressr
#'
#' @param s.star.list The list of values of
#'    \eqn{\hat{\bm{s}}^\star \equiv \bm{A}_{\mathrm{obs}}\hat{\bm{x}}_n^\star}
#'    in the cone-tightening procedure for each \eqn{\tau}.
#' @param tau.list The list of feasible parameters \eqn{\tau}.
#' @param maxR Maximum number of bootstrap replications to be considered in
#'    case there are some errors.
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the list of estimates from bootstrap:
#'   \item{T.bs}{A list of bootstrap test statistics
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'  \item{beta.bs.bar.list}{A list of \eqn{\tau_n}-tightened recentered
#'     bootstrap estimates \eqn{\bar{\bm{\beta}}^\star_{\mathrm{obs},n,b}}.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'
#' @export
#'
dkqs.bs <- function(data, lpmodel, beta.tgt, R, maxR, s.star.list, tau.list,
                    solver, progress, n) {
  # ---------------- #
  # Step 1: Initialize and assigning the lists
  # ---------------- #
  R.succ <- 0
  R.eval <- 0
  T.list <- list()
  beta.bs.bar.temp <- list()
  error.list <- list()
  tau.error <- list()

  # Check if there is any list objects in 'lpmodel'
  any.list <- lpmodel.anylist(lpmodel)

  # If there is some list objects, set maxR as the max length of the list
  if (isTRUE(any.list$list)) {
    maxR <- length(any.list$consol)
  }

  # Obtain the sample beta.obs with full data
  beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]

  # ---------------- #
  # Step 2: Bootstrap replicatoins
  # ---------------- #
  eval.count <- 0
  while ((R.succ < R) & (R.eval != maxR)) {

    # Evaluate the list of indices to be passed to 'future_lapply'
    bs.temp <- bs.assign(R, R.eval, R.succ, maxR, any.list, lpmodel, data,
                         n, TRUE)
    i0 <- bs.temp$i0
    i1 <- bs.temp$i1
    bs.list <- bs.temp$bs.list
    
    # Obtain results from the bootstrap replications
    progressr::with_progress({
      if (isTRUE(progress)) {
        pbar <- progressr::progressor(along = i0:i1)
      } else {
        pbar <- NULL
      }
      dkqs.return <- future.apply::future_lapply(bs.list,
                                                 FUN = dkqs.bs.fn,
                                                 future.seed = TRUE,
                                                 data = data,
                                                 lpmodel = lpmodel,
                                                 beta.obs = beta.obs.hat,
                                                 beta.tgt = beta.tgt,
                                                 s.star.list = s.star.list,
                                                 tau.list = tau.list,
                                                 n = n,
                                                 solver = solver,
                                                 pbar = pbar,
                                                 progress = progress,
                                                 eval.count = eval.count,
                                                 n.bs = i1 - i0 + 1)
      eval.count <- eval.count + 1
    })

    # Update the list and parameters
    post.return <- post.bs(dkqs.return, i0, i1, R.eval, T.list,
                           beta.bs.bar.temp, error.list, tau.list, tau.error)
    T.list <- post.return$T.list
    beta.bs.bar.temp <- post.return$beta.list
    tau.error <- post.return$error.param
    error.list <- post.return$error.list
    R.succ <- post.return$R.succ
    R.eval <- post.return$R.eval
  }

  # ---------------- #
  # Step 3: Consolidate the test statistics
  # ---------------- #
  T.bs <- list()
  for (i in seq_along(unlist(tau.list))) {
    # Make sure that the remainder is 0 if i == length(tau.list)
    if (i == length(tau.list)) {
      k <- 0
    } else {
      k <- i
    }
    T.bs[[i]] <- T.list[seq_along(T.list) %% length(tau.list) == k]
  }

  # ---------------- #
  # Step 4: Consolidate the error messages
  # ---------------- #
  if (R.eval != R.succ) {
    # Create data.frame for error messages
    df.error <- data.frame(id = NA,
                           tau = unlist(tau.error),
                           message = unlist(error.list))

    # Match the id of the error messages
    df.error <- error.id.match(error.list, df.error)
  } else {
    df.error <- NULL
  }

  return(list(T.bs = T.bs,
              beta.bs.bar.list = beta.bs.bar.temp,
              df.error = df.error,
              R.succ = R.succ))
}

#' Carries out one bootstrap replication for the DKQS procedure
#'
#' @description This function carries out the one bootstrap replication of the
#'   DKQS procedure This function is used in the \code{dkqs.bs} function via
#'   the \code{future_lapply} command.
#'
#' @inheritParams dkqs
#' @inheritParams dkqs.bs
#' @inheritParams dkqs.qlp
#' @param x This is either the list of indices that represent the bootstrap
#'   replications, or the list of bootstrap components of the \code{lpmodel}
#'   object passed from the user.
#' @param pbar Progress bar object.
#' @param eval.count Count for the number of times the \code{future_lapply}
#'   function has been called. If this object is zero, it means that the
#'   \code{future_lapply} function is being called for the first time in this
#'   subprocedure. Otherwise, it means that the \code{future_lapply} function
#'   has been called for more than once. This situation typically refers to the
#'   situations where there are some errors in the first time of the
#'   replications.
#' @param n.bs Total number of replications to be conducted in this procedure.
#'
#' @return Returns a list of output that are obtained from the DKQS
#'   procedure:
#'   \item{Ts}{Bootstrap test statistic.}
#'   \item{beta}{Bootstrap estimator.}
#'   \item{param}{List of problematic parameters in the DKQS test.}
#'   \item{msg}{Error message (if applicable).}
#'
#' @export
#'
dkqs.bs.fn <- function(x, data, lpmodel, beta.obs.hat, beta.tgt, s.star.list,
                       tau.list, solver, n, pbar, eval.count, n.bs, progress) {
  # ---------------- #
  # Step 1: Print progress bar
  # ---------------- #
  if (isTRUE(progress)) {
    if (eval.count == 0) {
      pbar(sprintf("(Computing %s bootstrap estimates)", n.bs))
    } else {
      pbar(sprintf("(Computing %s extra bootstrap estimates)", n.bs))
    }
  }

  # ---------------- #
  # Step 2: Initialize the parameters
  # ---------------- #
  # Initialization
  T.bs.list <- list()
  beta.bs.bar.list <- list()
  tau.error <- NULL
  msg <- NULL

  # Replace lpmodel by x if x is a list
  if (is.list(x)) {
    lpm <- lpmodel.update(lpmodel, x)
  } else {
    lpm <- lpmodel
  }

  # ---------------- #
  # Step 3: Conduct one bootstrap replication
  # ---------------- #
  # Draw data
  if (!is.null(data)) {
    data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
    rownames(data.bs) <- 1:nrow(data.bs)
  }

  # Bootstrap estimator
  beta.result <- tryCatch({
    beta.obs.bs <- lpmodel.beta.eval(data.bs, lpm$beta.obs, 1)[[1]]
    list(beta.obs.bs = beta.obs.bs)
  }, warning = function(w) {
    return(list(status = "warning",
                msg = w))
  }, error = function(e) {
    return(list(status = "error",
                msg = e))
  })

  if (is.null(beta.result$status)) {
    # Loop through each tau if there is no error in getting 'beta.obs'
    for (i in seq_along(unlist(tau.list))) {
      # Compute tau-tightened recentered bootstrap estimate
      beta.bs.bar <- beta.result$beta.obs.bs - beta.obs.hat + s.star.list[[i]]

      # Bootstrap test statistic
      T.result <- tryCatch({
        T.bs.i <- dkqs.qlp(lpm, beta.tgt, beta.bs.bar, tau.list[i], "cone",
                           n, solver)$objval
        list(T.bs.i = T.bs.i)
      }, warning = function(w) {
        return(list(status = "warning",
                    msg = w))
      }, error = function(e) {
        return(list(status = "error",
                    msg = e))
      })

      # Assign return list
      if (is.null(T.result$status)) {
        # Case 1: No error
        T.bs.list[[i]] <- T.result$T.bs.i
        beta.bs.bar.list[[i]] <- beta.bs.bar
        tau.error <- NULL
        msg <- NULL
      } else {
        # Case 2: No error in getting the bootstrap 'beta.obs' but there is an
        # error in solving the optimization problem with the bootstrap info
        T.bs.list <- NULL
        beta.bs.bar.list <- NULL
        tau.error <- tau.list[i]
        msg <- T.result$msg$message
        break()
      }
    }
  } else {
    # Case 3: There is an error in getting the bootstrap 'beta.obs'
    T.bs.list <- NULL
    beta.bs.bar.list <- NULL
    tau.error <- "NA"
    msg <- beta.result$msg$message
  }

  return(list(Ts = T.bs.list,
              beta = beta.bs.bar.list,
              param = tau.error,
              msg = msg))
}

#' Calculates the \eqn{p}-value
#'
#' @description This function computes the \eqn{p}-value of the test based on
#'    the bootstrap estimates.
#'
#' @param T_bs Bootstrap estimates of the test statistic.
#' @param T.n Sample test statistics.
#'
#' @return Returns the \eqn{p}-value:
#'   \item{p}{\eqn{p}-value.}
#'
#' @export
#'
pval <- function(T.bs, T.n) {
  # Compute p-value
  p <- mean(T.bs > T.n)

  return(list(p = p))
}

#' Creates the constraints for the linear program of \code{tau} in the
#' DKQS procedure
#'
#' @description This function generates the constraints in the linear program
#'   for computing the value of tau based on the procedure suggested by Kamat
#'   (2018) for the DKQS procedure, i.e. linear program (6) of
#'   Torgovitsky (2019).
#'
#' @param length.tau The number of variables in the constraint.
#' @param coeff.tau The coefficient in front of tau in the constraint.
#' @param coeff.x The coefficient in front of \eqn{x_i} in the constraint.
#' @param ind.x The index of \eqn{x_i}, i.e. the value of \eqn{i}.
#' @param rhs The rhs vector of the new constraint.
#' @param sense The equality or inequality symbol to be used in the new
#'    constraint.
#' @param lp.lhs.tau The constraint matrix to be updated.
#' @param lp.rhs.tau The rhs vector of the linear constraints to be updated.
#' @param lp.sense.tau The sense vector fo the linear constraints to be
#'    updated.
#'
#' @return Returns the list of matrices that corresponds to the updated
#'   constraints:
#'   \item{lp.lhs.tau}{Updated constraint matrix.}
#'   \item{lp.rhs.tau}{Updated rhs vector.}
#'   \item{lp.sense.tau}{Update sense for the constraints.}
#'
#' @export
#'
tau.constraints <- function(length.tau, coeff.tau, coeff.x, ind.x, rhs, sense,
                            lp.lhs.tau, lp.rhs.tau, lp.sense.tau) {
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

#' Checks and updates the input in \code{dkqs}
#'
#' @description This function checks and updates the input of the user. If
#'    there is any invalid input, this function will be terminated and
#'    generates appropriate error messages.
#'
#' @inheritParams dkqs
#'
#' @return Returns the updated parameters back to the function
#' \code{dkqs}. The following information are updated:
#'    \itemize{
#'       \item{\code{data}}
#'       \item{\code{lpmodel}}
#'       \item{\code{solver}}``
#'       \item{\code{solver.name}}
#'       \item{\code{test.logical}}
#'    }
#'
#' @export
#'
dkqs.check <- function(data, lpmodel, beta.tgt, R, Rmulti, tau, n, solver,
                       progress) {
  # ---------------- #
  # Step 1: Check the arguments
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
    check.nonnegative(tau[i], "tau")
  }

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
              test.logical = test.logical))
}

#' Print results from \code{dkqs}
#'
#' @description This function prints the \eqn{p}-values from \code{dkqs}.
#'
#' @param x Object returned from \code{dkqs}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.dkqs <- function(x, ...) {
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
      cat(sprintf(" p-value: %s\n", df.pval[1, 2]))
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
#' @description This function prints a summary of the results obtained from
#'   \code{dkqs}.
#'
#' @param x Objects returned from \code{dkqs}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.dkqs <- function(x, ...) {
  if (x$test.logical == 1) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-values
    print(x)

    # Print maximum feasible tau, test statistic, solver used used and the
    # number of bootstrap replications
    cat(sprintf(" Maximum feasible tau: %s\n", round(x$tau.max, digits = 5)))
    cat(sprintf(" Test statistic: %s\n", round(x$T.n, digits = 5)))
    cat(sprintf(" Solver used: %s\n", x$solver))
    cat(sprintf(" Number of successful bootstrap replications: %s\n",
                x$R.succ))

    # Number of failed bootstrap replications
    if (!is.null(x$df.error)) {
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
