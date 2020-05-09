#' Estimate bounds with shape restrictions
#'
#' @description This function computes the bound of the linear program
#'    subject to shape constraints. This function also offers an option
#'    to estimate the shape constraints using a two-step procedure and
#'    some tolerance level. \eqn{\ell^1}-norm and \eqn{\ell^2}-norm are
#'    supported in the estimation procedure.
#'
#' @import Matrix
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
#' @param norm Norm used in the optimization problem.
#' @param kappa Parameter used in the second step of the two-step procedure
#'    for obtaining the solution subject to the shape constraints.
#' @param estimate Boolean variable to indicate whether the estimated
#'    problem should be considered.
#' @inheritParams dkqs
#' @inheritParams invertci
#'
#' @return Returns the bounds subject to the shape constraints.
#'   \item{ub}{Upper bound with shape constraints}
#'   \item{lb}{Lower bound with shape constraints}
#'   \item{est}{Indicator of whether estimation is involved in the
#'   estimation}
#'   \item{call}{The function that has been called.}
#'   \item{norm}{Norm used in the optimization problem.}
#'
#' @export
#'
estbounds <- function(data, lpmodel, kappa = 1e-5, norm = 2, solver = NULL,
                      estimate = TRUE, progress = TRUE){

  # ---------------- #
  # Step 1: Obtain call, check and update input
  # ---------------- #
  # Obtain call information
  call = match.call()

  # Check and update
  estbounds.return = estbounds.check(data, lpmodel, kappa, norm, solver,
                                     estimate, progress)
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
  bound0infe <- FALSE

  ### Scenario 1: Estimate = FASLE, i.e. solve the exact problem
  if (estimate == FALSE){
    ub_shp0 <- estbounds.original(data, lpmodel, "max", solver)
    lb_shp0 <- estbounds.original(data, lpmodel, "min", solver)
    ub <- ub_shp0$objval
    lb <- lb_shp0$objval

    # Store indicator of whether the estimation procedure should be conducted
    if (is.numeric(ub) == FALSE | is.numeric(lb) == FALSE){
      bound0infe <- TRUE
      if (progress == TRUE){
        cat(paste("The original problem is infeasible. ",
                  "The estimated bounds will be displayed.", sep = ""))
      }
    } else {
      est = FALSE
    }
  }

  ### Scenario 2: Estimate = TRUE or scenario 1 is infeasible
  if (estimate == TRUE | bound0infe == TRUE){

    ## Solve model
    if (norm == 1){
      ## L1-norm
      # Stage one of the problem
      estbounds11 <- mincriterion(data, lpmodel, norm, solver.name)

      # Return stop message if there is no feasible solution for stage one
      # of the problem
      if (is.numeric(estbounds11$objval) == FALSE){
        stop("The constraints in the estimation problem are contradictory.
             Please ensure that the constraints are correctly specified.")
      }
      # Stage two of the problem
      estbounds_ub <- estbounds2.L1(data, estbounds11, lpmodel, "max", kappa,
                                    solver)
      estbounds_lb <- estbounds2.L1(data, estbounds11, lpmodel, "min", kappa,
                                    solver)
    } else if (norm == 2){
      ## L2-norm
      # Stage one of the problem
      estbounds21 <- mincriterion(data, lpmodel, norm, solver.name)

      # Return stop message if there is no feasible solution for stage one
      # of the problem
      if (is.numeric(estbounds21$objval) == FALSE){
        stop("The constraints in the estimation problem are contradictory.
             Please ensure that the constraints are correctly specified.")
      }
      # Stage two of the problem
      estbounds_ub <- estbounds2.L2(data, estbounds21, lpmodel, "max", kappa,
                                    solver)
      estbounds_lb <- estbounds2.L2(data, estbounds21, lpmodel, "min", kappa,
                                    solver)
    }

    # Store results
    ub <- estbounds_ub$objval
    lb <- estbounds_lb$objval

    est = TRUE
  }

  # ---------------- #
  # Step 3: Assign the return list and define class of output
  # ---------------- #
  output = list(ub = ub,
                lb = lb,
                est = est,
                call = call,
                norm = norm)

  attr(output, "class") = "estbounds"

  return(output)
}

#' Computes the true bounds with shape contraints
#'
#' @description The function computes the true bound subject to the shape
#'    constraints without approximation.
#'
#' @param original.sense Sense of the contraints to compute the true bound.
#' @inheritParams dkqs
#' @inheritParams invertci
#' @inheritParams estbounds
#'
#' @return Returns the solution to the linear program.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#'
#' @export
#'
estbounds.original <- function(data, lpmodel, original.sense, solver){

  # ---------------- #
  # Step 1: Problem set-up
  # ---------------- #
  # Matrices
  A.original <- rbind(lpmodel$A.obs, lpmodel$A.shp)

  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.hat <- lpmodel$beta.obs(data)
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix"){
    beta.obs.hat <- lpmodel$beta.obs
  }
  beta.original <- rbind(beta.obs.hat, lpmodel$beta.shp)
  # Sense contraints
  sense.original <- c(rep("=", nrow(A.original)))
  # Zero lower bound
  lb.zero <- rep(0, ncol(lpmodel$A.tgt))

  # ---------------- #
  # Step 2: Formulate the argument for optimization
  # ---------------- #
  oarg <- list(Af = NULL,
               bf = lpmodel$A.tgt,
               nf = NULL,
               A = A.original,
               rhs = beta.original,
               sense = sense.original,
               modelsense = original.sense,
               lb = lb.zero)

  # ---------------- #
  # Step 3: Solve the model
  # ---------------- #
  ans <- do.call(solver, oarg)

  # ---------------- #
  # Step 4: Return result
  # ---------------- #
  invisible(list(objval = ans$objval,
                 x = ans$x))
}

#' Estimates the bounds with shape contraints (Stage 2 with \eqn{\ell^1}-norm)
#'
#' @description This function evaluates the solution to stage 2 of the
#'    two-step procedure obtaining the estimated bound. \eqn{\ell^1}-norm
#'    is used in the constraint
#'
#' @param firststepsoln List of solutions to the first step problem.
#' @inheritParams gurobi.optim
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the second step of the two-step procedure.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#'
#' @export
#'
estbounds2.L1 <- function(data, firststepsoln, lpmodel, modelsense, kappa,
                          solver){

  # ---------------- #
  # Step 1: Initialization
  # ---------------- #
  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.hat <- lpmodel$beta.obs(data)
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix"){
    beta.obs.hat <- lpmodel$beta.obs
  }
  k <- length(beta.obs.hat)

  # ---------------- #
  # Step 2: Extract information from the first-stage solution
  # ---------------- #
  #### Step 1: Extract information from the first-stage solution
  Qhat <- firststepsoln$objval
  larg <- firststepsoln$larg

  # ---------------- #
  # Step 3: Construct the inequality constraint
  # ---------------- #
  # Update the linear constraint
  c <- larg$bf
  A.step2 <- rbind(larg$A, c)
  b.step2 <- c(larg$rhs, Qhat * (1+kappa))
  sense.step2 <- c(larg$sense, "<=")

  # Append the matrices to the list
  larg$A <- A.step2
  larg$rhs <- b.step2
  larg$sense <- sense.step2

  # ---------------- #
  # Step 4: Update objective function
  # ---------------- #
  # Update the objective matrix
  A.tgt.new <- cbind(lpmodel$A.tgt,
                     matrix(rep(0,
                                2*k*nrow(lpmodel$A.tgt)),
                            nrow = nrow(lpmodel$A.tgt)))
  larg$Af <- 0
  larg$bf <- A.tgt.new
  larg$nf <- 0

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
              x = step2.ans$x))
}

#' Estimates the bounds with shape contraints (Stage 2 with \eqn{\ell^2}-norm)
#'
#' @description This function evaluates the solution to stage 2 of the
#'    two-step procedure obtaining the estimated bound. \eqn{\ell^2}-norm
#'    is used in the constraint
#'
#' @param firststepsoln List of solutions to the first step problem.
#' @inheritParams gurobi.optim
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the second step of the two-step procedure.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#'
#' @export
#'
estbounds2.L2 <- function(data, firststepsoln, lpmodel, modelsense, kappa,
                          solver){
  # ---------------- #
  # Step 1: Extract information from the first-stage solution
  # ---------------- #
  Qhat <- firststepsoln$objval
  larg <- firststepsoln$larg

  # ---------------- #
  # Step 2: Construct the quadratic inequality constraint
  # ---------------- #
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.hat <- lpmodel$beta.obs(data)
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix"){
    beta.obs.hat <- lpmodel$beta.obs
  }
  step2_qc <- list()
  if (is.null(lpmodel$A.obs) == FALSE){
    step2_qc$Qc <- t(lpmodel$A.obs) %*% lpmodel$A.obs
    step2_qc$q <- -2*t(lpmodel$A.obs) %*% beta.obs.hat
    step2_qc$rhs <- Qhat * (1+kappa) - t(beta.obs.hat) %*% beta.obs.hat
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
  larg$bf <- lpmodel$A.tgt
  larg$nf <- 0

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
              x = step2_ans$x))
}

#' Checks and updates the input of the function \code{estbounds}
#'
#' @description This function checks and updates the input from the user for
#'    the function \code{estbounds}. If there is any invalid input, this
#'    function will terminate the procedure and generate appropriate error
#'    messages.
#'
#' @inheritParams estbounds
#' @inheritParams dkqs
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
estbounds.check <- function(data, lpmodel, kappa, norm, solver, estimate,
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
                           A.shp.cat = 1,
                           beta.obs.cat = c(2,3),
                           beta.shp.cat = 1,
                           R = 1)

  # Check solver
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # Check numerics
  norm <- check.norm(norm, "norm")
  kappa <- check.numrange(kappa, "kappa", "closed", 0, "closed", 1)

  # Check Boolean
  check.boolean(progress, "progress")
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
print.estbounds <- function(x, ...){
  cat("Call:\n")
  dput(x$call)
  cat("\n")

  if (x$est == TRUE){
    # Case 1: Report the estimated bounds
    if (is.numeric(x$norm) == TRUE){
      cat(sprintf("Norm used in optimization problem: L%s-norm \n", x$norm))
    } else {
      cat(sprintf("Norm used in optimization problem: %s-norm \n", x$norm))
    }
    cat(sprintf("Estimated bounds subject to shape constraints: [%s, %s] \n",
                round(x$lb, digits = 5), round(x$ub, digits = 5)))
  } else {
    # Case 2: Report the true bounds
    cat(sprintf("True bounds subject to shape constraints: [%s, %s] \n",
                x$lb, x$ub))
  }
}

#' Summary of results from \code{estbounds}
#'
#' @description This function uses the summary method on the return list of the
#'    function \code{estbounds}. This is a wrapper of the function
#'    \code{print.estbounds}.
#'
#' @param x Object returned from \code{estbounds}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned. This function prints results from
#'    \code{estbounds}.
#'
#' @export
#'
summary.estbounds <- function(x, ...){
  #### Call theprint function
  print(x)
}

#' First-stage of the estimation procedure for \code{estbounds}
#'
#' @description This function evaluates the solution to stage 1 of the
#'    two-step procedure obtaining the estimated bound. This function can
#'    be used to evaluate both the estimation problem with the 1-norm or
#'    the 2-norm.
#'
#' @inheritParams estbounds
#' @inheritParams dkqs
#'
#' @return Returns the solution to the first step of the two-step procedure
#'    and argument for the linear program.
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{larg}{Arguments for the estimation program.}
#'  \item{norm}{Norm used in the estimation problem.}
#'  \item{solver}{The solver used in the estimation problem}
#'  \item{call}{The details of the function that has been called.}
#'
#' @export
#'
mincriterion <- function(data, lpmodel, norm, solver){
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
  # Check if beta_obs is a function, then compute the
  if (class(lpmodel$beta.obs) == "function"){
    beta.obs.hat <- lpmodel$beta.obs(data)
  } else if (class(lpmodel$beta.obs) == "numeric" |
             class(lpmodel$beta.obs) == "matrix"){
    beta.obs.hat <- lpmodel$beta.obs
  }

  # ---------------- #
  # Step 3: Create common constraints for the problem with 1-norm and 2-norm
  # ---------------- #
  # Zero lower bound
  lb.zero <- rep(0,ncol(lpmodel$A.tgt))
  # Generate the sense of models
  sense0 <- rep("=", nrow(lpmodel$A.shp))

  # ---------------- #
  # Step 4: Set up argument for the optimizer
  # ---------------- #
  if (norm == 1){
    # Define the augmented matrices
    k <- length(beta.obs.hat)
    # Introduce slack variables into the matrix
    A.aug <- cbind(lpmodel$A.shp, matrix(rep(0, 2*k*nrow(lpmodel$A.shp)),
                                         nrow = nrow(lpmodel$A.shp)))
    A.slack <- cbind(lpmodel$A.obs, -diag(k), -diag(k))
    # Combine the constraints
    A.new <- rbind(A.aug, A.slack)
    beta.new <- c(lpmodel$beta.shp, beta.obs.hat)
    # New model sense
    sense.new <- c(sense0, rep("=", k))
    # New objective function
    c <- c(rep(0, dim(lpmodel$A.obs)[2]), rep(1, k), rep(1, k))
    # New lower bound
    lb.new <- rep(0, length(c))
    # 1-norm
    optim.arg <- list(Af = NULL,
                      bf = c,
                      nf = NULL,
                      A = A.new,
                      rhs = beta.new,
                      sense = sense.new,
                      modelsense = "min",
                      lb = lb.new)

  } else if (norm == 2){
    # 2-norm
    optim.arg <- list(Af = lpmodel$A.obs,
                      bf = beta.obs.hat,
                      nf = 1,
                      A = lpmodel$A.shp,
                      rhs = lpmodel$beta.shp,
                      sense = sense0,
                      modelsense = "min",
                      lb = lb.zero,
                      qc = NULL)
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

#' Checks and updates the input of the function \code{mincriterion}
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
mincriterion.check <- function(data, lpmodel, norm, solver){
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
                           A.shp.cat = 1,
                           beta.obs.cat = c(2,3),
                           beta.shp.cat = 1,
                           R = 1)

  # Check solver
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # Check numerics
  norm <- check.norm(norm, "norm")

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
#' @description This function uses the print method on the return list of the
#'    function \code{mincriterion}.
#'
#' @param x Object returned from \code{mincriterion}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned. This function prints results from
#'    \code{mincriterion}.
#'
#' @export
#'
print.mincriterion <- function(x, ...){
  cat(sprintf("Minimum value: %s \n", round(x$objval)))
  cat(sprintf("Norm: %s \n", x$norm))
  cat(sprintf("Solver: %s \n", x$solver))

}

#' Summary of results from \code{mincriterion}
#'
#' @description This function uses the summary method on the return list of the
#'    function \code{mincriterion}. This is a wrapper of the function
#'    \code{print.mincriterion}.
#'
#' @param x Object returned from \code{mincriterion}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned. This function prints results from
#'    \code{mincriterion}.
#'
#' @export
#'
summary.mincriterion <- function(x, ...){
  #### Call theprint function
  print(x)
}
