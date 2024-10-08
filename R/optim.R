#' LP and QP solver by \code{Gurobi}
#'
#' @description This function computes the solution to the linear or quadratic
#'    program using the \code{Gurobi} package. This function can have linear
#'    and/or quadratic constraints.
#'
#' @param Af The matrix that is involved in the objective function.
#' @param bf The vector that is involved in the objective function.
#' @param nf The number of observations in the data frame.
#' @param A The constraint matrix.
#' @param rhs The rhs vector for the linear constraints.
#' @param modelsense The indicator of whether the model is to max or min an
#'   objective function.
#' @param sense The sense of the linear constraints.
#' @param lb The lower bound vector.
#' @param qc The list of quadratic constraint(s). There can be multiple
#'   quadratic constraints. Each constraint has to be a list.
#' @param weight The weighting matrix.
#' @param ... List of options to be passed to the Gurobi solver. This part is
#'    optional.
#'
#' @return Returns the optimal point, the optimal value and the status of the
#'  solution.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'  \item{status}{The status of the optimization problem.}
#'
#' @export
#'
gurobi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb, qc = NULL,
                         weight = NULL, beta.r.program = FALSE, ...) {
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf, weight)

  # ---------------- #
  # Step 2: Gurobi set-up
  # ---------------- #
  model <- list()
  # Objective function - Quadratic / list
  model$Q <- objective_return$obj2
  model$obj <- smatrixconvert(objective_return$obj1)
  model$objcon <- as.numeric(objective_return$obj0)

  # Linear constraints
  model$A <- dmatrixconvert(A)
  model$rhs <- smatrixconvert(rhs)

  # Quadratic constraints
  model$quadcon <- qc

  # Model sense and lower bound
  model$sense <- smatrixconvert(sense)
  model$modelsense <- modelsense
  model$lb <- smatrixconvert(lb)

  # ---------------- #
  # Step 3: Result of the linear or quadratic program, and return result
  # ---------------- #
  params <- list(...)
  # Append parameters
  params <- append(params,
                   list(OutputFlag = 0,
                        PSDTol = Inf))
  if (length(params) == 0) {
    result <- gurobi::gurobi(model)
  } else {
    result <- gurobi::gurobi(model, params)
  }

  # ---------------- #
  # Step 4: Try through different options if the status is 'ITERATION_LIMIT'
  # ---------------- #
  if ((result$status == "ITERATION_LIMIT") |
      ((result$status != "OPTIMAL") & (isTRUE(beta.r.program)))) {
    gurobi.method <- c(0, 1, 2)
    gurobi.numericfocus <- c(999, 1, 2, 3)
    gurobi.scaleflag <- c(999, -1, 0, 1, 2, 3)
    for (i in seq_along(gurobi.scaleflag)) {
      params.new <- params
      if (gurobi.scaleflag[i] != 999) {
        params.new$ScaleFlag <- gurobi.scaleflag[i]
      }
      for (j in seq_along(gurobi.numericfocus)) {
        if (gurobi.numericfocus[j] != 999) {
          params.new$NumericFocus <- gurobi.numericfocus[j]
        }
        for (k in seq_along(gurobi.method)) {
          params.new$Method <- gurobi.method[k]
          result <- gurobi::gurobi(model, params.new)
          if (result$status == "OPTIMAL") {
            break()
          }
        }
        if (result$status == "OPTIMAL") {
          break()
        }
      }
      if (result$status == "OPTIMAL") {
        break()
      }
    }
  }

  return(list(objval = as.numeric(result$objval),
              x = as.numeric(result$x),
              status = result$status))
}

#' LP and QP solver by \code{cplexAPI}
#'
#' @description This function computes the solution to the quadratic and linear
#'    programs using the \code{cplexAPI} package.
#'
#' @inheritParams gurobi.optim
#'
#' @return Returns the optimal point, the optimal value and the status of the
#'  solution.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'  \item{status}{The status of the optimization problem.}
#'
#' @export
#'
cplexapi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb,
                           weight = NULL, ...) {
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf, weight)

  # ---------------- #
  # Step 2: Update the notations
  # ---------------- #
  # Model sense
  modelsense[modelsense == "min"] <- cplexAPI::CPX_MIN
  modelsense[modelsense == "max"] <- cplexAPI::CPX_MAX

  # Inequality/equality signs
  sense[sense == "<="] <- "L"
  sense[sense == ">="] <- "G"
  sense[sense == "=="] <- "E"
  sense[sense == "="] <- "E"

  # Bounds
  lb[lb == Inf] <- cplexAPI::CPX_INFBOUND
  lb[lb == -Inf] <- -cplexAPI::CPX_INFBOUND
  ub <- rep(cplexAPI::CPX_INFBOUND, length(lb))

  # ---------------- #
  # Step 3: cplexAPI environment
  # ---------------- #
  # Model environment
  env <- cplexAPI::openEnvCPLEX()
  cplexAPI::setDblParmCPLEX(env, 1016, 1e-06)
  prob <- cplexAPI::initProbCPLEX(env)
  cplexAPI::chgProbNameCPLEX(env, prob, "sample")

  # Constraint matrices
  cnt <- apply(A, MARGIN = 2, function(x) length(which(x != 0)))
  beg <- rep(0, ncol(A))
  beg[-1] <- cumsum(cnt[-length(cnt)])
  ind <- unlist(apply(A, MARGIN = 2, function(x) which(x != 0) - 1))
  val <- c(A)
  val <- val[val != 0]

  # ---------------- #
  # Step 4: Solve the problem
  # ---------------- #
  # A linear program is identified if obj2 == NULL
  if (is.null(objective_return$obj2) == TRUE) {
    # Solving linear program
    cplexAPI::copyLpwNamesCPLEX(env,
                                prob,
                                ncol(A),
                                nrow(A),
                                modelsense,
                                objective_return$obj1,
                                rhs,
                                sense,
                                beg,
                                cnt,
                                ind,
                                val,
                                lb,
                                ub)
    cplexAPI::lpoptCPLEX(env, prob)
    solution <- cplexAPI::solutionCPLEX(env, prob)
  } else {
    # Solving quadratic program
    stop("This version can only solve linear programs by CPLEX at the moment.
         Please use another solver for quadratic programs")
  }
  cplexAPI::closeEnvCPLEX(env)

  # Status code
  status <- solution$lpstat
  status.msg <- sprintf("Status code: %s", status)
  return(list(objval = as.numeric(solution$objval),
              x = as.numeric(solution$x),
              status = status.msg))
}

#' LP and QP solver by \code{Rcplex}
#'
#' @description This function computes the solution to the linear and quadratic
#'    programs using the \code{Rcplex} package.
#'
#' @inheritParams gurobi.optim
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the optimal point, the optimal value and the status of the
#'  solution.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'  \item{status}{The status of the optimization problem.}
#'
#' @export
#'
rcplex.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb,
                         weight = NULL, ...) {
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf, weight)

  # ---------------- #
  # Step 2: Update vectors and sense
  # ---------------- #
  # Update sense
  sense[sense == ">="] <- "G"
  sense[sense == "<="] <- "L"
  sense[sense == "="]  <- "E"

  # Define upper bound
  ub <- rep(Inf, length(lb))

  # Define Q matrix
  # - Keep Qmat as NULL for linear program
  # - Multiply obj2 by 2 for Q for quadratic program to offset the 1/2 factor
  if (is.null(objective_return$obj2) == TRUE) {
    Qmat <- objective_return$obj2
  } else {
    Qmat <- 2*objective_return$obj2
  }

  # ---------------- #
  # Step 3: Solve model
  # ---------------- #
  solution <- Rcplex::Rcplex(cvec = t(smatrixconvert(objective_return$obj1)),
                             Amat = A,
                             bvec = as.matrix(smatrixconvert(rhs)),
                             Qmat = Qmat,
                             lb = smatrixconvert(lb),
                             sense = sense,
                             ub = smatrixconvert(ub),
                             objsense = modelsense,
                             vtype = "C",
                             n = 1)

  # ---------------- #
  # Step 4: Update and return result
  # ---------------- #
  if (is.null(objective_return$obj0) == FALSE) {
    objval <- solution$obj + objective_return$obj0
  } else {
    objval <- solution$obj
  }
  status <- solution$status
  status.msg <- sprintf("Status code: %s", status)

  return(list(objval = as.numeric(objval),
              x = as.numeric(solution$xopt),
              status = status.msg))
}

#' LP and QP solver by \code{limSolve}
#'
#' @description This function computes the solution to linear and quadratic
#'    programs using the \code{limSolve} package.
#'
#' @inheritParams gurobi.optim
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the optimal point, the optimal value and the status of the
#'  solution.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'  \item{status}{The status of the optimization problem.}
#'
#' @export
#'
limsolve.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb,
                           weight = NULL, ...) {
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf, weight)

  # ---------------- #
  # Step 2: Update lower bounds
  # ---------------- #
  # Change the lower bounds to inequality constriants
  lb_Amat <- diag(length(lb))
  lb_bvec <- lb

  # Update constraint matrices
  A <- rbind(A, lb_Amat)
  rhs <- as.matrix(Reduce(rbind, c(rhs, lb_bvec)))
  sense <- c(sense, rep(">=", length(lb_bvec)))

  # ---------------- #
  # Step 3: Update constraints
  # ---------------- #
  # Objective function
  if (modelsense == "max") {
    fcost <- -objective_return$obj1
  } else if (modelsense == "min") {
    fcost <- objective_return$obj1
  }

  # Equality constraints
  Emat <- A[sense == "=",]
  Fvec <- rhs[sense == "="]

  # Inequality constraint >=
  Gmat1 <- A[sense == ">=",]
  Hvec1 <- rhs[sense == ">="]

  # Inequality constraint <=
  Gmat2 <- -A[sense == "<=",]
  Hvec2 <- -rhs[sense == "<="]

  # Combine G and h matrices
  Gmat <- rbind(Gmat1, Gmat2)
  Hvec <- Reduce(rbind, c(Hvec1, Hvec2))

  # ---------------- #
  # Step 4: Solve the model
  # ---------------- #
  # Linear solver is used if obj2 is a zero matrix (i.e. number of zeros equals
  # the total number of elements) or NULL
  if (is.null(objective_return$obj2) == TRUE |
      sum(objective_return$obj2 == 0) == length(objective_return$obj2)) {
    ### Linear program solver
    solution <- limSolve::linp(E = smatrixconvert(Emat),
                               F = smatrixconvert(Fvec),
                               G = smatrixconvert(Gmat),
                               H = smatrixconvert(Hvec),
                               Cost = smatrixconvert(fcost))
    # Obtain objective function, and add back the constant term, negate the
    # solution if it is a max problem
    if (modelsense == "max") {
      objval <- -solution$solutionNorm + objective_return$obj0
    } else if (modelsense == "min") {
      objval <- solution$solutionNorm + objective_return$obj0
    }
  } else {
    if (modelsense == "min") {
      ### Quadratic program solver
      # Formulate the two matrices
      Amat <- Af * sqrt(nf)
      Bvec <- bf * sqrt(nf)
      solution <- limSolve::lsei(A = Amat, B = Bvec, E = smatrixconvert(Emat),
                                 F = Fvec, G = smatrixconvert(Gmat), H = Hvec)

      # Obtain objective function
      objval <- solution$solutionNorm
    } else if (modelsense == "max") {
      stop("This package cannot be used to solve max problems.")
    }
  }
  # Optimal the optimal value of x
  x <- solution$X

  # Status of the optimization problem
  status <- solution$IsError
  if (status == TRUE) {
    status.msg <- "No error has occurred."
  } else {
    status.msg <- paste0("An error has occurred in solving the optimization ",
                         "problem by limSolve")
  }

  return(list(x = as.numeric(x),
              objval = as.numeric(objval),
              status = status.msg))
}


#' Computes the coefficient terms of the objective functions
#'
#' @description This function computes the matrices in the objective functions
#'    for linear programs. This function takes matrix \eqn{\bm{A}} and
#'    \eqn{\bm{\beta}} as input and computes the coefficients of the objective
#'    function.
#'
#' @importFrom Matrix t
#'
#' @param A The matrix \eqn{\bm{A}}.
#' @param b The column vector \eqn{\bm{\beta}}.
#' @param n The sample size \eqn{n}.
#' @inheritParams gurobi.optim
#'
#' @details
#' \itemize{
#'   \item{\strong{Quadratic programs} ---
#'      Given inputs \eqn{\bm{A} \in \mathbf{R}^{m\times m}},
#'      \eqn{\bm{b} \in \mathbf{R}^m} and
#'      \eqn{\bm{W} \in \mathbf{R}^{m\times m}}, the equation of the objective
#'      function of the quadratic program can be written as
#'      \deqn{n (\bm{A}\bm{x} -\bm{\beta})'\bm{W}(\bm{A}\bm{x}-\bm{\beta})
#'      = n\bm{x}'\bm{A}'\bm{W}\bm{A}\bm{x} - 2n\bm{\beta}'\bm{W}\bm{A}\bm{x} +
#'      n\bm{\beta}'\bm{W}\bm{\beta}.}
#'      If the \eqn{\bm{W}} matrix is not specified, then it will be taken as
#'      an identity matrix of order \eqn{m}.}
#'   \item{\strong{Linear programs} ---
#'      For all linear problems that are considered in this code, one of
#'      \eqn{\bm{A}} and \eqn{\bm{b}} is \code{NULL} or is a zero vector. The
#'      term that is nonzero and nonnull will be multiplied by \eqn{n} and
#'      used as \code{obj1}.}
#' }
#'
#' @return Returns the following three quantities: \code{obj2} is the
#'   coefficient of the quadratic term, \code{obj1} is the coefficient of the
#'   linear term and \code{obj0} is the constant term. More explicitly, their
#'   form are given as follows:
#'   \item{obj2}{The coefficient for the second-order term. It is
#'     returned as \code{NULL} for linear programs and
#'     \eqn{n\bm{A}'\bm{W}\bm{A}} for quadratic programs.}
#'   \item{obj1}{The coefficient term of the linear term. For quadratic
#'     programs, it is returned as \eqn{-2n\bm{\beta}'\bm{W}\bm{A}}.}
#'   \item{obj0}{The constant term of the linear program. For quadratic
#'     programs, it is returned as \eqn{n\bm{\beta}'\bm{W}\bm{\beta}}.}
#'
#' @export
#'
objective.function <- function(A, b, n, weight = NULL) {
  # If-else function to determine if it corresponds to a linear or quadratic
  # program. This is identified by whether one of A and b is null or nonzero
  # because it would be the case for linear programs that are considered in
  # this package.
  if (is.null(A) == TRUE | sum(A == 0) == length(A)) {
    # Linear program coefficients with nonzero vector b
    obj2 <- NULL
    obj1 <- b * n
    obj0 <- 0
  } else if (is.null(b) == TRUE | sum(b == 0) == length(b)) {
    # Linear program coefficients with nonzero matrix A
    obj2 <- NULL
    obj1 <- A * n
    obj0 <- 0
  } else {
    if (is.null(weight)) {
      weight <- diag(length(b))
    }

    # Quadratic program coefficients
    obj2 <- asmat(Matrix::t(A) %*% weight %*% A * n)
    rownames(obj2) <- 1:nrow(obj2)
    colnames(obj2) <- 1:ncol(obj2)
    obj1 <- -2 * asmat(Matrix::t(b) %*% weight %*% A * n)
    obj0 <- Matrix::t(b) %*% weight %*% b * n
  }

  # Return the above objective functions
  return(list(obj2 = obj2,
              obj1 = obj1,
              obj0 = obj0))
}

#' LP solver by \code{lpSolveAPI}
#'
#' @description This function computes the solution to the linear program
#'    using the \code{lpSolveAPI} package.
#'
#' @inheritParams gurobi.optim
#' @inheritParams dkqs
#' @inheritParams dkqs.qlp
#'
#' @return Returns the optimal point, the optimal value and the status of the
#'  solution.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'  \item{status}{The status of the optimization problem.}
#'
#' @details The package \code{lpSolveAPI} cannot be used to solve quadratic
#'   programs.
#'
#' @export
#'
lpsolveapi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb,
                             weight = NULL, ...) {
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf, weight)

  # ---------------- #
  # Step 2: Update the constraint matrices
  # ---------------- #
  # Change the lower bounds to inequality constriants
  lb_Amat <- diag(length(lb))
  lb_bvec <- lb
  # Update constraint matrices
  A <- smatrixconvert(rbind(A, lb_Amat))
  rhs <- smatrixconvert(Reduce(rbind, c(rhs, lb_bvec)))
  sense <- c(sense, rep(">=", length(lb_bvec)))

  # ---------------- #
  # Step 3: LP formulation
  # ---------------- #
  # solve object
  lprec <- lpSolveAPI::make.lp(nrow = nrow(A), ncol = ncol(A))
  # Model sense
  lpSolveAPI::lp.control(lprec, sense = modelsense)
  # Types of decision variables
  lpSolveAPI::set.type(lprec, 1:ncol(A), type = c("real"))
  lpSolveAPI::set.objfn(lprec, smatrixconvert(objective_return$obj1))
  # Define the constraints
  for (i in 1:nrow(A)) {
    lpSolveAPI::add.constraint(lprec, A[i, ], sense[i], rhs[i])
  }

  # ---------------- #
  # Step 4: Solve and obtain solution of LP
  # ---------------- #
  solve(lprec)
  x <- lpSolveAPI::get.variables(lprec)
  objval <- lpSolveAPI::get.objective(lprec)
  status <- lpSolveAPI::solve.lpExtPtr(lprec)
  status.msg <- sprintf("Status code: %s", status)

  # ---------------- #
  # Step 5: Return results
  # ---------------- #
  return(list(objval = objval,
              x = x,
              status = status.msg))
}
