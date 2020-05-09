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
#'   \item{pval}{\eqn{p}-value.}
#'   \item{tau}{The value of tau used \eqn{\tau^\ast} in the linear and 
#'      quadratic programs.}
#'   \item{T.n}{Test statistic \eqn{T.n}.}
#'   \item{T.bs}{The list of bootstrap test statistics 
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'   \item{beta.bs.bar}{The list of \eqn{\tau}-tightened re-centered bootstrap 
#'      estimators \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'   \item{lb0}{Logical lower bound of the problem.}
#'   \item{ub0}{Logical upper bound of the problem.}
#'   \item{solver}{Solver used in solving the linear and quadratic programs.}
#'   \item{cores}{Number of cores used.}
#'   \item{call}{The function that has been called.}
#' 
#' @details If the value of the test statistic \eqn{T.n} is zero, the
#'    bootstrap procedure will be skipped.
#' 
#' @export
#' 
dkqs <- function(data, lpmodel, beta.tgt, R = 100, tau = .5, solver = NULL, 
                 cores = 1, progress = FALSE){
  
  # ---------------- #
  # Step 1: Update call, check and update the arguments
  # ---------------- #
  # Obtain call information
  call <- match.call()
  
  # Check the arguments
  dkqs.return <- dkqs.check(data, lpmodel, beta.tgt, R, tau, solver,
                               cores, progress)
  
  # Update the arguments
  data <- dkqs.return$data
  lpmodel <- dkqs.return$lpmodel
  solver <- dkqs.return$solver
  solver.name <- dkqs.return$solver.name
  cores <- dkqs.return$cores
  
  # ---------------- #
  # Step 2: Initialization
  # ---------------- #
  n <- nrow(data)
  J <- length(unique(data$Y)) - 1
  
  # Initialize beta.obs
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.hat <- lpmodel$beta.obs(data)
  } else if (class(lpmodel$beta.obs) == "list"){
    beta.obs.hat <- lpmodel$beta.obs[[1]]
  }
  
  # ---------------- #
  # Step 3: Choose the value of tau
  # ---------------- #
  tau.return <- dkqs.qlp(data, lpmodel, beta.tgt, beta.obs.hat, tau, "tau", 
                         n, solver)
  if (tau > tau.return$objval){
    tau <- tau.return$objval
  } else if (tau <= tau.return$objval){
    tau <- tau
  } else {
    # Error message when the problem is infeasible.
    stop("The problem is infeasible. Choose other values of tau.")
  }
  
  # ---------------- #
  # Step 4: Compute T.n, x.star and s.star
  # ---------------- #
  # Compute T.n
  T.n <- dkqs.qlp(data, lpmodel, beta.tgt, beta.obs.hat, tau, "test", 
                  n, solver)$objval
  
  # The problem is infeasible if T.n is NULL
  if (is.null(T.n) == TRUE){
    stop("The problem is infeasible. Choose other values of tau.")
  }
  
  # Compute s.star
  x.return <- dkqs.qlp(data, lpmodel, beta.tgt, beta.obs.hat, tau, "cone", 
                       n, solver)
  x.star <- x.return$x
  s.star <- lpmodel$A.obs %*% x.star
  
  # ---------------- #
  # Step 5: Compute the bootstrap beta and estimates
  # ---------------- #
  # if (T.n != 0){
  if (TRUE){
    if (cores == 1){
      # No parallelization
      T.bs.return <- beta.bs(data, lpmodel, beta.tgt, R, J, s.star,
                             tau, solver, progress)
    } else {
      # Parallelization
      T.bs.return <- beta.bs.parallel(data, lpmodel, beta.tgt, R, J, s.star,  
                                      tau, solver, progress)
    }
    # Retrive answer
    T.bs <- T.bs.return$T.bs
    beta.bs.bar <- T.bs.return$beta.bs.bar.list
  } else {
    if (progress == TRUE){
      cat(paste("Bootstrap is skipped because the value of the test ",
                "statistic is zero.\n"))
      T.n <- 0
      T.bs.return <- NULL
      beta.bs.bar <- NULL
    }
  }
  # ---------------- #
  # Step 6: Compute p-value
  # ---------------- #
  p <- pval(T.bs, T.n)$p
  
  # ---------------- #
  # Step 7: Obtain logical bounds for the function invertci
  # ---------------- #
  lb0 <- x.return$lb0
  ub0 <- x.return$ub0
  
  # ---------------- #
  # Step 8: Close the progress bar that is used in the bootstrap procedure
  # ---------------- #
  if (progress == TRUE){
    close(T.bs.return$pb) 
  }
  
  # ---------------- #
  # Step 9: Assign the return list and return output
  # ---------------- #
  output <- list(pval = as.numeric(p), 
                 tau = as.numeric(tau), 
                 T.n = as.numeric(T.n), 
                 T.bs = T.bs,
                 beta.bs.bar = beta.bs.bar,
                 lb0 = lb0$objval,
                 ub0 = ub0$objval,
                 solver = solver.name,
                 cores = cores,
                 call = call)
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
dkqs.qlp <- function(data, lpmodel, beta.tgt, beta.obs.hat, tau, problem, n, 
                     solver){
  # ---------------- #
  # Step 1: Obtain the dimension of the 
  # ---------------- #
  A.tgt.nr <- nrow(lpmodel$A.tgt)
  A.tgt.nc <- ncol(lpmodel$A.tgt)
  
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
    for (i in 1:length(ind.up)){
      new.const <- tau_constraints(len.tau, rhs.up, -1, ind.up[i]+1, 0, "<=",
                                   lp.lhs.tau, lp.rhs.tau, lp.sense.tau)
      lp.lhs.tau <- new.const$lp.lhs.tau
      lp.rhs.tau <- new.const$lp.rhs.tau
      lp.sense.tau <- new.const$lp.sense.tau
    }
    # Inequality constraints for ind.down
    for (i in 1:length(ind.down)){
      new.const <- tau_constraints(len.tau, rhs.down, -1, ind.down[i]+1, 0,
                                   "<=", lp.lhs.tau, lp.rhs.tau, lp.sense.tau)
      lp.lhs.tau <- new.const$lp.lhs.tau
      lp.rhs.tau <- new.const$lp.rhs.tau
      lp.sense.tau <- new.const$lp.sense.tau
    }
    # Inequality constraints for ind.0
    for (i in 1:length(ind.0)){
      new.const <- tau_constraints(len.tau, rhs.0, -1, ind.0[i]+1, 0, "<=",
                                   lp.lhs.tau, lp.rhs.tau, lp.sense.tau)
      lp.lhs.tau <- new.const$lp.lhs.tau
      lp.rhs.tau <- new.const$lp.rhs.tau
      lp.sense.tau <- new.const$lp.sense.tau
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
#' @param s_star The value of 
#'    \eqn{\hat{s}^\ast \equiv A_{\mathrm{obs}}\hat{\bm{x}}_n^\ast} 
#'    in the cone-tightening procedure.
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the list of estimates from bootstrap:
#'   \item{T.bs}{A list of bootstrap test statistics 
#'      \eqn{\{\overline{T}_{n,b}(\tau_n)\}^B_{b=1}}.}
#'  \item{beta.bs.bar.list}{A list of \eqn{\tau_n}-tightened recentered 
#'     bootstrap estimates \eqn{\bar{\beta}^\ast_{\mathrm{obs},n,b}}.}
#'  \item{pb}{Progress bar object.}
#'
#' @export
#' 
beta.bs <- function(data, lpmodel, beta.tgt, R, J, s.star, tau, solver,
                    progress){
  # ---------------- #
  # Step 1: Initialize the vectors and progress counters
  # ---------------- #
  T.bs <- NULL
  beta.bs.bar.list <- NULL
  
  # Initialize the progress bar
  if (progress == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = R, style = 3, width = 20)
    cat("\r") 
  } else {
    pb <- NULL
  }
  
  # ---------------- #
  # Step 2: Bootstrap procedure
  # ---------------- #
  lpmodel.bs <- lpmodel
  
  # Use the for-loop to construct the bootstrap test statistic
  for (i in 1:R){
    data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
    rownames(data.bs) <- 1:nrow(data.bs)
    
    # Compute the bootstrap test statistic
    if (class(lpmodel$beta.obs) == "function"){
      beta.obs.bs <- lpmodel$beta.obs(data.bs)
      beta.obs <- lpmodel$beta.obs(data)
    } else if (class(lpmodel$beta.obs) == "list") {
      beta.obs.bs <- lpmodel$beta.obs[[i+1]]
      beta.obs <- lpmodel$beta.obs[[1]]
    }
    
    # Compute beta.bs.bar and test statistic
    beta.bs.bar <- beta.obs.bs - beta.obs + s.star
    T.bs.i <- dkqs.qlp(data.bs, lpmodel.bs, beta.tgt, beta.bs.bar, tau, "cone", 
                       nrow(data), solver)$objval
    
    # Append results
    T.bs <- c(T.bs, T.bs.i)
    beta.bs.bar.list <- cbind(beta.bs.bar.list, beta.bs.bar)
    
    # Update progress bar
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
  # Step 3: Return results
  # ---------------- #
  return(list(T.bs = T.bs,
              beta.bs.bar.list = beta.bs.bar.list,
              pb = pb))
}

#' Computes the bootstrap test statistics with parallelization
#' 
#' @description This function computes the bootstrap test statistics with 
#'    parallelized for-loops.
#' 
#' @import doMC foreach
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
#' 
#' @export
#' 
beta.bs.parallel <- function(data, lpmodel, beta.tgt, R, J, s.star, tau,
                             solver, progress, cores){
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
  
  # ---------------- #
  # Step 2: Initialize progress bar, comb function and assign doRNG
  # ---------------- #
  if (progress == TRUE){
    # Initialize the counter
    cl <- PtProcess::makeSOCKcluster(8)
    doSNOW::registerDoSNOW(cl)
    
    # Set the counter and progress bar
    pb <- utils::txtProgressBar(max = R, style = 3, width = 20) 
    cat("\r")
    progress <- function(n){
      utils::setTxtProgressBar(pb, n) 
      if (n < R){
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
  
  # Assign doRnG
  `%dorng%` <- doRNG::`%dorng%`
  
  # ---------------- #
  # Step 3: Bootstrap procedure
  # ---------------- #
  listans <- foreach(i = 1:R, .multicombine = TRUE, .combine = "comb", 
                     .options.snow = opts, .packages = "lpinfer") %dorng% {
    lpmodel.bs <- lpmodel
   
    # Re-sample the data 
    data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
    rownames(data.bs) <- 1:nrow(data.bs)
   
    # Compute the bootstrap test statistic
    if (class(lpmodel$beta.obs) == "function"){
      beta.obs.bs <- lpmodel$beta.obs(data.bs)
      beta.obs <- lpmodel$beta.obs(data)
    } else if (class(lpmodel$beta.obs) == "list") {
      beta.obs.bs <- lpmodel$beta.obs[[i+1]]
      beta.obs <- lpmodel$beta.obs[[1]]
    }
   
    # Compute beta.bs.bar and test statistic
    beta.bs.bar <- beta.obs.bs - beta.obs + s.star
    T.bs.i <- dkqs.qlp(data.bs, lpmodel.bs, beta.tgt, beta.bs.bar, tau, "cone", 
                       nrow(data), solver)$objval
   
    # Append results
    T.bs <- data.frame(T.bs.i)
    beta.bs.bar.list <- data.frame(beta.bs.bar)
    list(T.bs, beta.bs.bar.list)
  }
  
  # ---------------- #
  # Step 4: Retrieve information from the output
  # ---------------- #
  T.bs <- as.vector(unlist(listans[[1]]))
  beta.bs.bar.list <- data.frame(matrix(unlist(listans[[2]]), 
                                        nrow = beta.obs.nr, 
                                        ncol = R,
                                        byrow = FALSE))
  
  # ---------------- #
  # Step 5: Return results
  # ---------------- #
  return(list(T.bs = T.bs,
              beta.bs.bar.list = beta.bs.bar.list,
              pb = pb))
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
pval <- function(T.bs, T.n, alpha = .05){
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
dkqs.check <- function(data, lpmodel, beta.tgt, R, tau, solver, cores, 
                       progress){
  # ---------------- #
  # Step 1: Check the arguments
  # ---------------- #
  # Check data
  data <- check.dataframe(data, "data")
  
  # Check lpmodel 
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel, 
                           name.var = "lpmodel",
                           A.tgt.cat = 1,
                           A.obs.cat = 1,
                           A.shp.cat = 0,
                           beta.obs.cat = c(2,3),
                           beta.shp.cat = 0,
                           R = R)
  
  # Check solver
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name
  
  # Check numerics 
  check.numeric(beta.tgt, "beta.tgt")
  check.positiveinteger(R, "R")
  cores <- check.positiveinteger(cores, "cores")
  
  # Check Boolean 
  check.boolean(progress, "progress")
  
  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              cores = cores))
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
  cat(sprintf("Test statistic: %s.             \n", round(x$T.n, digits = 5)))  
  cat(sprintf("p-value: %s.\n", round(x$pval, digits = 5)))
  cat(sprintf("Value of tau used: %s.\n", round(x$tau, digits = 5)))
  cat(sprintf("Solver used: %s.\n", x$solver))
  cat(sprintf("Number of cores used: %s.\n", x$cores))
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
  print(x)
}


