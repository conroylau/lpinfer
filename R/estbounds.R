#' Estimate bounds with shape restrictions
#'
#' @description This function computes the bound of the estimation problem
#'    subject to shape constraints. This function offers an option to estimate
#'    the shape constraints using a two-step procedure with the 1-norm or
#'    2-norm.
#'
#' @param lpmodel The \code{lpmodel} object used in the test that consists of
#'   the following objects:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{A.shp}}
#'      \item{\code{beta.obs}}
#'      \item{\code{beta.shp}}
#'    }
#' @param norm Norm used in the optimization problem. It can be either a 1-norm
#'   or a 2-norm.
#' @param kappa Tuning parameter used in the second step of the two-step
#'    procedure for obtaining the bounds subject to the shape constraints.
#'    It can be any nonnegative number.
#' @param estimate Boolean variable to indicate whether the bounds should be
#'    estimated or not.
#' @inheritParams dkqs
#' @inheritParams invertci
#'
#' @details The following input for \code{norm} will be interpreted as the
#'   1-norm:
#'   \itemize{
#'     \item{\code{1} (\code{numeric})}
#'     \item{\code{"1"} (\code{string})}
#'     \item{\code{"L1"}}
#'     \item{\code{"one"}}
#'     \item{\code{"o"}}
#'     \item{\code{"taxicab"}}
#'   }
#'   The following input for \code{norm} will be interpreted as the 2-norm:
#'   \itemize{
#'     \item{\code{2} (\code{numeric})}
#'     \item{\code{"2"} (\code{string})}
#'     \item{\code{"L2"}}
#'     \item{\code{"two"}}
#'     \item{\code{"t"}}
#'     \item{\code{"e"}}
#'     \item{\code{"euclidean"}}
#'   }
#'
#' @return Returns the bounds subject to the shape constraints.
#'   \item{ub}{Upper bound.}
#'   \item{lb}{Lower bound.}
#'   \item{mincriterion}{Objective value of the first-stage problem (i.e.
#'      \code{mincriterion}).}
#'   \item{est}{Indicator of whether estimation is involved in the
#'   estimation}
#'   \item{call}{The function that has been called.}
#'   \item{norm}{Norm used in the optimization problem.}
#'   \item{ub.status}{Status of the optimization problem for the upper bound.}
#'   \item{lb.status}{Status of the optimization problem for the lower bound.}
#'   \item{solver}{Name of the solver used.}
#'
#' @example ./inst/example/estbounds_example.R
#'
#' @export
#'
estbounds <- function(data = NULL, lpmodel, kappa = 0, norm = 2,
                      estimate = TRUE, solver = NULL) {
  # ---------------- #
  # Step 1: Obtain call, check and update input
  # ---------------- #
  # Obtain call information
  call = match.call()

  # Check and update
  estbounds.return <- estbounds.check(data, lpmodel, kappa, norm, solver,
                                      estimate)
  # Update the input
  data <- estbounds.return$data
  lpmodel <- estbounds.return$lpmodel
  solver <- estbounds.return$solver
  solver.name <- estbounds.return$solver.name
  norm <- estbounds.return$norm

  # ---------------- #
  # Step 2: Construct the bounds
  # ---------------- #
  # Default - Boolean variable of whether the answer to the scenario 1 is
  # feasible or not
  if (isFALSE(estimate)) {
    ### Scenario 1: Estimate = FALSE, i.e. solve the exact problem
    ub_shp0 <- estbounds.original(data, lpmodel, "max", solver)
    lb_shp0 <- estbounds.original(data, lpmodel, "min", solver)

    # Get the solution to the bounds
    ub <- ub_shp0$objval
    lb <- lb_shp0$objval

    # Status of the solutions
    ub.status <- ub_shp0$status
    lb.status <- lb_shp0$status

    # Return error message in case the identified set is empty
    non.error <- c("OPTIMAL", "UNBOUNDED")
    if (!((ub.status %in% non.error) & (lb.status %in% non.error))) {
      warning(paste0("The identified set is empty. Please check ",
                     "'ub.status' and 'lb.status' for the status codes."))
    }
    # Assign the parameters
    minc.objval <- NULL
    est <- FALSE
  } else {
    ### Scenario 2: Estimate = TRUE
    ## Solve model
    contra.msg <- paste0("The constraints in the estimation problem ",
                         "are contradictory. Please ensure that the ",
                         "constraints are correctly specified.")
    if (norm == 1) {
      ## L1-norm
      # Stage one of the problem
      minc <- mincriterion(data, lpmodel, norm, solver.name)

      # Return stop message if there is no feasible solution for stage one
      # of the problem
      if (is.numeric(minc$objval) == FALSE) {
        stop(contra.msg)
      }
      # Stage two of the problem
      estbounds_ub <- estbounds2.L1(data, minc, lpmodel, "max", kappa, solver)
      estbounds_lb <- estbounds2.L1(data, minc, lpmodel, "min", kappa, solver)
    } else if (norm == 2) {
      ## L2-norm
      # Stage one of the problem
      minc <- mincriterion(data, lpmodel, norm, solver.name)

      # Return stop message if there is no feasible solution for stage one
      # of the problem
      if (is.numeric(minc$objval) == FALSE) {
        stop(contra.msg)
      }
      # Stage two of the problem
      estbounds_ub <- estbounds2.L2(data, minc, lpmodel, "max", kappa, solver)
      estbounds_lb <- estbounds2.L2(data, minc, lpmodel, "min", kappa, solver)
    }

    # Store results
    ub <- estbounds_ub$objval
    lb <- estbounds_lb$objval
    ub.status <- estbounds_ub$status
    lb.status <- estbounds_lb$status

    est <- TRUE

    # Assign mincriterion
    minc.objval <- minc$objval
  }

  # ---------------- #
  # Step 3: Assign the return list and define class of output
  # ---------------- #
  output = list(ub = ub,
                lb = lb,
                mincriterion = minc.objval,
                est = est,
                call = call,
                norm = norm,
                ub.status = ub.status,
                lb.status = lb.status,
                solver = solver.name)

  attr(output, "class") = "estbounds"

  return(output)
}

#' Computes the true bounds with shape constraints
#'
#' @description This function computes the true bounds subject to the shape
#'    constraints without approximation.
#'
#' @param original.sense Sense of the constraints for the true bounds.
#' @inheritParams invertci
#' @inheritParams estbounds
#'
#' @return Returns the solution to the linear program.
#'  \item{objval}{Optimal objective value.}
#'  \item{status}{Status of the linear program.}
#'
#' @export
#'
estbounds.original <- function(data, lpmodel, original.sense, solver) {
  # ---------------- #
  # Step 1: Problem set-up
  # ---------------- #
  # Ensure A.tgt is matrix
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)
  if (!is.matrix(A.tgt.hat)) {
    A.tgt.matrix <- matrix(A.tgt.hat, nrow = 1)
  } else {
    A.tgt.matrix <- A.tgt.hat
  }
  A.tgt.nc <- ncol(A.tgt.matrix)

  # Ensure A.shp is matrix
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
  if (!is.matrix(A.shp.hat)) {
    A.shp.matrix <- matrix(A.shp.hat, nrow = 1)
  } else {
    A.shp.matrix <- A.shp.hat
  }

  # Matrices
  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
  A.original <- rbind(A.obs.hat, A.shp.matrix)
  if (!is.matrix(A.original)) {
    A.original <- matrix(A.original, nrow = 1)
  }

  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix" |
             class(lpmodel$beta.obs) == "list") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  }
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)
  beta.original <- c(beta.obs.hat, beta.shp.hat)

  # Sense constraints
  sense.original <- c(rep("=", nrow(A.original)))

  # Zero lower bound
  lb.zero <- rep(0, A.tgt.nc)

  # ---------------- #
  # Step 2: Formulate the argument for optimization
  # ---------------- #
  oarg <- list(Af = NULL,
               bf = A.tgt.matrix,
               nf = 1,
               A = A.original,
               rhs = beta.original,
               sense = sense.original,
               modelsense = original.sense,
               lb = lb.zero)

  # ---------------- #
  # Step 3: Solve the model
  # ---------------- #
  ans <- do.call(solver, oarg)

  # If the status is INF_OR_UNBD, set DualReductions as 0 and reoptimize
  if (ans$status == "INF_OR_UNBD") {
    oarg$DualReductions <- 0
    ans <- do.call(solver, oarg)
  }

  # If the status is UNBOUNDED, it tells nothing about feasibility so set
  # the objective as 0 to check feasibility
  if (ans$status == "UNBOUNDED") {
    oarg$nf <- 0
    ans2 <- do.call(solver, oarg)
  } else {
    ans2 <- list(status = "NA")
  }

  if (ans$status == "OPTIMAL") {
    # Case 1: Optimal
    objval <- ans$objval
    status <- ans$status
  } else if (ans$status == "UNBOUNDED" & ans2$status == "OPTIMAL") {
    # Case 2: Unbounded and the problem is feasible
    if (original.sense == "max") {
      objval <- Inf
    } else {
      objval <- -Inf
    }
    status <- ans$status
  } else {
    # Case 3: Infeasible
    if (original.sense == "max") {
      objval <- -Inf
    } else {
      objval <- Inf
    }
    if (ans2$status != "NA") {
      status <- ans2$status
    } else {
      status <- ans$status
    }
  }

  # ---------------- #
  # Step 4: Return result
  # ---------------- #
  invisible(list(objval = objval,
                 status = status))
}

#' Estimates the bounds with shape constraints (stage 2 with 1-norm)
#'
#' @description This function evaluates the solution to stage 2 of the
#'    two-step procedure to obtain the estimated bound with the 1-norm.
#'
#' @param firststepsoln List of solutions to the first step problem.
#' @inheritParams gurobi.optim
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the second step of the two-step procedure.
#'  \item{objval}{Optimal objective value.}
#'  \item{x}{Optimal point.}
#'  \item{status}{Status of the optimization problem.}
#'
#' @export
#'
estbounds2.L1 <- function(data, firststepsoln, lpmodel, modelsense, kappa,
                          solver) {
  # ---------------- #
  # Step 1: Initialization
  # ---------------- #
  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix" |
             class(lpmodel$beta.obs) == "list") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  }
  k <- length(beta.obs.hat)

  # ---------------- #
  # Step 2: Extract information from the first-stage solution
  # ---------------- #
  # Extract information from the first-stage solution
  Qhat <- firststepsoln$objval
  larg <- firststepsoln$larg

  # Ensure A.tgt is matrix
  A.tgt.matrix <- lpmodel.eval(data, lpmodel$A.tgt, 1)
  A.tgt.nr <- nrow(A.tgt.matrix)

  # ---------------- #
  # Step 3: Construct the inequality constraint
  # ---------------- #
  # Update the linear constraint
  c <- larg$bf
  A.step2 <- rbind(larg$A, c)
  if (!is.matrix(A.step2)) {
    A.step2 <- matrix(A.step2, nrow = 1)
  }
  b.step2 <- c(larg$rhs, Qhat * (1 + kappa))
  sense.step2 <- c(larg$sense, "<=")

  # Append the matrices to the list
  larg$A <- A.step2
  larg$rhs <- b.step2
  larg$sense <- sense.step2

  # ---------------- #
  # Step 4: Update objective function
  # ---------------- #
  # Update the objective matrix
  A.tgt.new <- cbind(A.tgt.matrix,
                     matrix(rep(0, 2 * k * A.tgt.nr), nrow = A.tgt.nr))
  larg$Af <- 0
  larg$bf <- A.tgt.new
  larg$nf <- 1

  # ---------------- #
  # Step 5: Update model sense based on max or min in step 3
  # ---------------- #
  larg$modelsense <- modelsense

  # ---------------- #
  # Step 6: Solve the model
  # ---------------- #
  step2.ans <- do.call(solver, larg)

  # ---------------- #
  # Step 7: Return results
  # ---------------- #
  return(list(objval = step2.ans$objval,
              x = step2.ans$x,
              status = step2.ans$status))
}

#' Estimates the bounds with shape constraints (Stage 2 with 2-norm)
#'
#' @description This function evaluates the solution to stage 2 of the
#'    two-step procedure to obtain the estimated bound with the 2-norm.
#'
#' @param firststepsoln List of solutions to the first step problem.
#' @inheritParams gurobi.optim
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the second step of the two-step procedure.
#'  \item{objval}{Optimal objective value.}
#'  \item{x}{Optimal point.}
#'  \item{status}{Status of the optimization problem.}
#'
#' @export
#'
estbounds2.L2 <- function(data, firststepsoln, lpmodel, modelsense, kappa,
                          solver) {
  # ---------------- #
  # Step 1: Extract information from the first-stage solution
  # ---------------- #
  Qhat <- firststepsoln$objval
  larg <- firststepsoln$larg

  # ---------------- #
  # Step 2: Construct the quadratic inequality constraint
  # ---------------- #
  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)

  if (class(lpmodel$beta.obs) == "function") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix" |
             class(lpmodel$beta.obs) == "list") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  }
  step2_qc <- list()
  if (is.null(A.obs.hat) == FALSE) {

    step2_qc$Qc <- t(A.obs.hat) %*% A.obs.hat
    step2_qc$q <- as.vector(-2 * t(A.obs.hat) %*% beta.obs.hat)
    step2_qc$rhs <- Qhat * (1 + kappa) - t(beta.obs.hat) %*% beta.obs.hat
    step2_qc$sense <- "<="
  } else {
    step2_qc <- NULL
  }
  qc_all <- list(step2_qc)

  # Update quadratic constraint
  larg$qc <- qc_all

  # ---------------- #
  # Step 3: Update objective function
  # ---------------- #
  larg$Af <- 0
  larg$bf <- A.tgt.hat
  larg$nf <- 1

  # ---------------- #
  # Step 4: Update model sense based on max or min in step 2
  # ---------------- #
  larg$modelsense <- modelsense

  # ---------------- #
  # Step 5: Solve the model
  # ---------------- #
  step2_ans <- do.call(solver, larg)

  # ---------------- #
  # Step 6: Return results
  # ---------------- #
  return(list(objval = step2_ans$objval,
              x = step2_ans$x,
              status = step2_ans$status))
}

#' Checks and updates the input in \code{estbounds}
#'
#' @description This function checks and updates the input from the user for
#'    the function \code{estbounds}. If there is any invalid input, this
#'    function will terminate the procedure and generate appropriate error
#'    messages.
#'
#' @inheritParams estbounds
#' @inheritParams invertci
#'
#' @return Returns the updated parameters and objects back to the function
#' \code{estbounds}. The following information are updated:
#'    \itemize{
#'       \item{\code{data}}
#'       \item{\code{lpmodel}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{norm}}
#'       \item{\code{kappa}}
#'    }
#'
#' @export
#'
estbounds.check <- function(data, lpmodel, kappa, norm, solver, estimate) {
  # ---------------- #
  # Step 1: Check the arguments
  # ---------------- #
  # Check lpmodel
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = c("matrix", "function_mat", "list"),
                           A.obs.cat = c("matrix", "function_mat", "list"),
                           A.shp.cat = c("matrix", "function_mat", "list"),
                           beta.obs.cat = c("function_mat",
                                            "list_vector",
                                            "matrix",
                                            "function_obs_var"),
                           beta.shp.cat = c("matrix", "function_mat", "list"),
                           R = 1,
                           is.estbounds = TRUE)

  # Check data
  if (!is.null(data)) {
    data <- check.dataframe(data)
  }

  # Check norm
  norm <- check.norm(norm, "norm")

  # Check solver
  if (norm == 1) {
    solver.return <- check.solver(solver, "solver", norm, FALSE)
  } else if (norm == 2) {
    solver.return <- check.solver(solver, "solver", norm, TRUE)
  }
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # Check kappa
  kappa <- check.nonnegative(kappa, "kappa")

  # Check Boolean
  check.boolean(estimate, "estimate")

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              norm = norm,
              kappa = kappa))
}

#' Print results from \code{estbounds}
#'
#' @description This function uses the print method on the return list of the
#'    function \code{estbounds}.
#'
#' @param x Object returned from \code{estbounds}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned. This function prints results from
#'    \code{estbounds}.
#'
#' @export
#'
print.estbounds <- function(x, ...) {
  # Print the estimated bounds, norm used and the solver used
  if (x$est == TRUE) {
    # Case 1: Report the estimated bounds
    cat(sprintf("Estimated bounds: [%s, %s] \n",
                round(x$lb, digits = 5), round(x$ub, digits = 5)))
  } else {
    # Case 2: Report the true bounds
    if (x$lb == Inf & x$ub == -Inf) {
      cat("The identified set is empty.\n")
    } else {
      cat(sprintf("True bounds: [%s, %s] \n", x$lb, x$ub))
    }
  }
}

#' Summary of results from \code{estbounds}
#'
#' @description This function prints a summary of the results obtained from
#'   \code{estbounds}.
#'
#' @param x Objects returned from \code{estbounds}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.estbounds <- function(x, ...) {
  # Print the estimated bounds, norm used and the solver used
  if (x$est == TRUE) {
    # Case 1: Report the estimated bounds
    cat(sprintf("Estimated bounds: [%s, %s] \n",
                round(x$lb, digits = 5), round(x$ub, digits = 5)))
    cat(sprintf("Norm used: %s \n", x$norm))
  } else {
    # Case 2: Report the true bounds
    print(x)
  }
  cat(sprintf("Solver: %s \n", x$solver))
}

#' First-stage estimation procedure for \code{estbounds}
#'
#' @description This function evaluates the solution to stage 1 of the
#'    two-step procedure to obtain the estimated bounds. This function can
#'    be used to evaluate both the estimation problem with the 1-norm or
#'    the 2-norm.
#'
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the first step of the two-step procedure
#'    and argument for the linear program.
#'  \item{objval}{Optimal objective value.}
#'  \item{x}{Optimal point.}
#'  \item{larg}{Arguments for the estimation program.}
#'  \item{norm}{Norm used in the estimation problem.}
#'  \item{solver}{The solver used in the estimation problem}
#'  \item{call}{The details of the function that has been called.}
#'
#' @example ./inst/example/mincriterion_example.R
#'
#' @export
#'
mincriterion <- function(data = NULL, lpmodel, norm = 2, solver = NULL) {
  # ---------------- #
  # Step 1: Obtain call, check and update the dependencies
  # ---------------- #
  # Obtain the call information
  call <- match.call()

  # Check the arguments
  checkupdate <- mincriterion.check(data, lpmodel, norm, solver)

  # Update the arguments
  data <- checkupdate$data
  lpmodel <- checkupdate$lpmodel
  solver <- checkupdate$solver
  solver.name <- checkupdate$solver.name
  norm <- checkupdate$norm

  # ---------------- #
  # Step 2: Obtain beta_obs and update solver
  # ---------------- #
  A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
  A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
  A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)

  # Count the dimension of matrices
  A.tgt.dim <- dim(A.tgt.hat)
  if (is.null(A.tgt.dim)) {
    A.tgt.nc <- length(A.tgt.hat)
  } else {
    A.tgt.nc <- A.tgt.dim[2]
  }

  # Count the dimension of matrices
  A.shp.dim <- dim(A.shp.hat)
  if (is.null(A.shp.dim)) {
    A.shp.nr <- 1
  } else {
    A.shp.nr <- A.shp.dim[1]
  }

  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix" |
             class(lpmodel$beta.obs) == "list") {
    beta.obs.hat <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)[[1]]
  }

  # ---------------- #
  # Step 3: Create common constraints for the problem with 1-norm and 2-norm
  # ---------------- #
  # Zero lower bound
  lb.zero <- rep(0, A.tgt.nc)
  # Generate the sense of models
  sense0 <- rep("=", A.shp.nr)

  # ---------------- #
  # Step 4: Set up argument for the optimizer
  # ---------------- #
  beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)

  if (norm == 1) {
    # Define the augmented matrices
    k <- length(beta.obs.hat)
    # Introduce slack variables into the matrix
    if (!is.matrix(A.shp.hat)) {
      A.shp.mat <- matrix(A.shp.hat, nrow = 1)
    } else {
      A.shp.mat <- A.shp.hat
    }
    A.aug <- cbind(A.shp.mat, matrix(rep(0, 2 * k * A.shp.nr),
                                         nrow = A.shp.nr))
    A.slack <- cbind(A.obs.hat, -diag(k), diag(k))
    # Combine the constraints
    A.new <- rbind(A.aug, A.slack)
    beta.new <- c(beta.shp.hat, beta.obs.hat)
    # New model sense
    sense.new <- c(sense0, rep("=", k))
    # New objective function
    c <- c(rep(0, dim(A.obs.hat)[2]), rep(1, k), rep(1, k))
    # New lower bound
    lb.new <- rep(0, length(c))
    # 1-norm
    optim.arg <- list(Af = NULL,
                      bf = c,
                      nf = 1,
                      A = A.new,
                      rhs = beta.new,
                      sense = sense.new,
                      modelsense = "min",
                      lb = lb.new)
  } else if (norm == 2) {
    if (!is.matrix(A.shp.hat)) {
      A.shp.new <- matrix(A.shp.hat, nrow = 1)
    } else {
      A.shp.new <- A.shp.hat
    }
    # 2-norm
    optim.arg <- list(Af = A.obs.hat,
                      bf = beta.obs.hat,
                      nf = 1,
                      A = A.shp.new,
                      rhs = beta.shp.hat,
                      sense = sense0,
                      modelsense = "min",
                      lb = lb.zero)
  }

  # ---------------- #
  # Step 5: Solve the model
  # ---------------- #
  ans <- do.call(solver, optim.arg)

  # ---------------- #
  # Step 6: Assign the return list and define class of output
  # ---------------- #
  output <- list(objval = ans$objval,
                 x = ans$x,
                 larg = optim.arg,
                 norm = norm,
                 solver = solver.name,
                 call = call)

  attr(output, "class") <- "mincriterion"

  return(output)
}

#' Checks and updates the input in \code{mincriterion}
#'
#' @description This function checks and updates the input from the user for
#'    the function \code{mincriterion}. If there is any invalid input, this
#'    function will terminate the procedure and generate appropriate error
#'    messages.
#'
#' @inheritParams estbounds
#' @inheritParams dkqs
#' @inheritParams invertci
#'
#' @return Returns the updated parameters back to the function
#' \code{mincriterion}. The following information are updated:
#'    \itemize{
#'       \item{\code{data}}
#'       \item{\code{lpmodel}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{norm}}
#'    }
#'
#' @export
#'
mincriterion.check <- function(data, lpmodel, norm, solver) {
  # ---------------- #
  # Step 1: Check the arguments
  # ---------------- #
  # Check lpmodel
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = c("matrix", "function_mat", "list"),
                           A.obs.cat = c("matrix", "function_mat", "list"),
                           A.shp.cat = c("matrix", "function_mat", "list"),
                           beta.obs.cat = c("function_mat",
                                            "list_vector",
                                            "matrix",
                                            "function_obs_var"),
                           beta.shp.cat = c("matrix", "function_mat", "list"),
                           R = 1,
                           is.estbounds = TRUE)

  # Check data
  if (!is.null(data)) {
    data <- check.dataframe(data)
  }

  # Check numerics
  norm <- check.norm(norm, "norm")

  # Check solver
  solver.return <- check.solver(solver, "solver", norm, FALSE)
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # ---------------- #
  # Step 2: Return results
  # ---------------- #
  return(list(data = data,
              lpmodel = lpmodel,
              solver = solver,
              solver.name = solver.name,
              norm = norm))
}

#' Print results from \code{mincriterion}
#'
#' @description This function prints a summary of the results obtained from
#'   \code{mincriterion}.
#'
#' @param x Object returned from \code{mincriterion}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.mincriterion <- function(x, ...) {
  # Print the minimum value
  cat(sprintf("Minimum value: %s \n", round(x$objval)))
}

#' Summary of results from \code{mincriterion}
#'
#' @description This function uses the summary method on the return list of the
#'    function \code{mincriterion}.
#'
#' @param x Objects returned from \code{mincriterion}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.mincriterion <- function(x, ...) {
  # Print the minimum value, norm used and solver
  cat(sprintf("Minimum value: %s \n", round(x$objval, digits = 5)))
  cat(sprintf("Norm used: %s \n", x$norm))
  cat(sprintf("Solver: %s \n", x$solver))
}
