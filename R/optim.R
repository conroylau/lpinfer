#' LP and QP solver by Gurobi
#'
#' @description This function computes the solution to the quadratic or linear
#'    program using the `\code{Gurobi}' package. This function can have linear
#'    and/or quadratic constraints.
#'
#' @param Af The matrix that is involved in the objective function.
#' @param bf The vector that is involved in the objective function.
#' @param nf The number of observations in the dataframe.
#' @param A The constraint matrix.
#' @param rhs The RHS vector for the linear constraints.
#' @param modelsense The indicator of whether the model is to max or min an
#'    objective function.
#' @param sense The sense of the linear constraints.
#' @param lb The lower lound vector.
#' @param qc List of quadratic constraint(s). There can be multiple quadratic
#'    constraints. Each constraint has to be a list.
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
#' 
gurobi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb, qc = NULL){
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf)

  # ---------------- #
  # Step 2: Gurobi set-up
  # ---------------- #
  model <- list()
  # Objective function - Quadratic / list
  model$Q <- objective_return$obj2
  model$obj <- objective_return$obj1
  model$objcon <- objective_return$obj0
  
  # Linear constraints
  model$A <- A
  model$rhs <- rhs
  
  # Quadrtaic constraints
  model$quadcon <- qc 
  
  # Model sense and lower bound
  model$sense <- sense 
  model$modelsense <- modelsense
  model$lb <- lb
  
  # ---------------- #
  # Step 3: Result of the linear or quadratic program, and return result
  # ---------------- #
  params <- list(OutputFlag=0, FeasibilityTol=1e-9)
  solution <- gurobi::gurobi(model, params)
  return(list(objval = as.numeric(solution$objval),
              x = as.numeric(solution$x)))
}

#' LP and QP solver by cplexAPI
#'
#' @description This function computes the solution to the quadratic and linear
#'    programs using the `\code{cplexAPI}' package.
#'    
#' @inheritParams cplexAPI
#' @inheritParams dkqs
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
#' 
cplexapi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf)
  
  # ---------------- #
  # Step 2: Update the notations
  # ---------------- #
  # Model sense
  modelsense[modelsense == "min"] <- CPX_MIN
  modelsense[modelsense == "max"] <- CPX_MAX
  
  # Inequality/equality signs
  sense[sense == "<="] <- "L"
  sense[sense == ">="] <- "G"
  sense[sense == "=="] <- "E"
  sense[sense == "="] <- "E"
  
  # Bounds
  lb[lb == Inf] <- cplexAPI::CPX_INFBOUND
  lb[lb == -Inf] <- -cplexAPI::CPX_INFBOUND
  ub <- rep(CPX_INFBOUND, length(lb))
  
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
  if (is.null(obj2) == TRUE){
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
         Please use another solver for quadratic progarms.")
  }
  cplexAPI::closeEnvCPLEX(env)
  return(list(objval = as.numeric(solution$objval),
              x = as.numeric(solution$x)))
}

#' LP and QP solver by Rcplex
#'
#' @description This function computes the solution to the linear and quadratic
#'    programs using the `\code{Rcplex}' package.
#'
#' @inheritParams gurobi_optim
#' @inheritParams dkqs
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'
#' @export
#' 
rcplex.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return = objective.function(Af, bf, nf)
  
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
  if (is.null(objective_return$obj2) == TRUE){
    Qmat <- objective_return$obj2
  } else {
    Qmat <- 2*objective_return$obj2
  }
  
  # ---------------- #
  # Step 3: Solve model
  # ---------------- #
  solution <- Rcplex::Rcplex(cvec = t(objective_return$obj1),
                             Amat = A, 
                             bvec = rhs,
                             Qmat = Qmat,
                             lb = lb,
                             sense = sense,
                             ub = ub,
                             objsense = modelsense,
                             vtype = "C",
                             n = 1)
  
  # ---------------- #
  # Step 4: Update and return result
  # ---------------- #
  if (is.null(objective_return$obj0) == FALSE){
    objval <- solution$obj + objective_return$obj0
  } else {
    objval <- solution$obj
  }
  return(list(objval = as.numeric(objval),
              x = as.numeric(solution$xopt)))
}

#' LP and QP solver by limSolve
#' 
#' @description This function computes the solution to linear and quadratic 
#'    programs using the `\code{limSolve}' package.
#'
#' @inheritParams gurobi_optim
#' @inheritParams dkqs
#' @inheritParams prog_cone
#'
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#'   
#' @export
#' 
limsolve.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf)
  
  # ---------------- #
  # Step 2: Update lower bounds
  # ---------------- #
  # Change the lower bounds to inequality constriants
  lb_Amat <- diag(length(lb))
  lb_bvec <- lb
  
  # Update constraint matrices
  A <- rbind(A, lb_Amat)
  rhs <- c(rhs, lb_bvec)
  sense <- c(sense, rep(">=", length(lb_bvec)))
  
  # ---------------- #
  # Step 3: Update constraints
  # ---------------- #
  # Objective function
  if (modelsense == "max"){
    fcost <- - objective_return$obj1 
  } else if (modelsense == "min"){
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
  Hvec <- as.matrix(c(c(Hvec1), c(Hvec2)), ncol = 1, byrow = TRUE)
  
  # ---------------- #
  # Step 4: Solve the model
  # ---------------- #
  # Linear solver is used if obj2 is a zero matrix (i.e. number of zeros equals 
  # the total number of elements) or NULL
  if (is.null(objective_return$obj2) == TRUE | 
      sum(objective_return$obj2 == 0) == length(objective_return$obj2)){
    ### Linear program solver
    solution <- limSolve::linp(E = Emat, F = Fvec, G = Gmat, H = Hvec, 
                              Cost = fcost)
    
    # Obtain objective function, and add back the constant term, negate the 
    # solution if it is a max problem
    if (modelsense == "max"){
      objval <- -solution$solutionNorm + objective_return$obj0      
    } else if (modelsense == "min"){
      objval <- solution$solutionNorm + objective_return$obj0      
    }
  } else {
    if (modelsense == "min"){
      ### Quadratic program solver
      # Formulate the two matrices
      Amat <- Af * sqrt(nf)
      Bvec <- bf * sqrt(nf)
      solution <- limSolve::lsei(A = Amat, B = Bvec, E = Emat, F = Fvec, 
                                G = Gmat, H = Hvec)
      
      # Obtain objective function
      objval <- solution$solutionNorm
    } else if (modelsense == "max"){
      stop("This package cannot be used to solve max problems.")
    }
  }
  # Optimal the optimal value of x
  x <- solution$X
  return(list(x = as.numeric(x),
              objval = as.numeric(objval)))
}


#' Auxiliary function to return the coefficient terms of the objective 
#' functions
#' 
#' @description This function computes the matrics in the objective functions 
#'    for linear programs. This function takes matrix \eqn{\bm{A}} and 
#'    \eqn{\bm{\beta}} as input and computes the coefficients of the objective 
#'    function. 
#'    
#' @param A The matrix \eqn{\bm{A}},
#' @param beta The column vector \eqn{\bm{\beta}}.
#' @param n The sample size \eqn{n}.
#' 
#' @details 
#' \itemize{
#'   \item{\strong{Quadratic programs} ---
#'      Given inputs \eqn{\bm{A} \in \mathbf{R}^{m\times m}} and 
#'      \eqn{\bm{b} \in \mathbf{R}^m}, the equation of the objective function 
#'      of the quadratic program can be written as
#'      \deqn{n (\bm{A}\bm{x} -\bm{\beta})'(\bm{A}\bm{x}-\bm{\beta})
#'      = n\bm{x}'\bm{A}'\bm{A}\bm{x} - 2n\bm{A}'\bm{\beta}\bm{x} + 
#'      n\bm{\beta}'\bm{\beta}.}}
#'   \item{\strong{Linear programs} ---
#'      For all linear problems that are considered in this code, one of 
#'      \eqn{\bm{A}} and \eqn{\bm{b}} is \code{NULL} or is a zero vector. The 
#'      term that is nonzero and nonnull will be used as \code{obj1}.}
#' }
#' 
#' @return Returns the following three quantities: \code{obj2} is the 
#'   coefficient of the quadratic term, \code{obj1} is the coefficient of the 
#'   linear term and \code{obj0} is the constant term. More explicitly, their
#'   form are given as follows:
#'   \item{obj2}{This is the coefficient for the second-order term. It is 
#'     returned as \code{NULL} for linear programs and \eqn{n\bm{A}'\bm{A}} for
#'     quadratic programs.}
#'   \item{obj1}{This is the coefficient term of the linear term. For quadratic 
#'     programs, it is returned as \eqn{-2n\bm{A}'\bm{\beta}}.}
#'   \item{obj0}{This is the constant term of the linear program. For quadratic
#'     programs, it is returned as \eqn{n\bm{\beta}'\bm{\beta}}.}
#'     
#' @export
#' 
objective.function <- function(A, b, n){
  # If-else function to determine if it corresponds to a linear or quadratic
  # program. This is identified by whether one of A and b is null or nonzero
  # because it would be the case for linear programs that are considered in
  # this package.
  if (is.null(A) == TRUE | sum(A == 0) == length(A)){
    # Linear program coefficients with nonzero vector b
    obj2 <- NULL
    obj1 <- b
    obj0 <- 0
  } else if (is.null(b) == TRUE | sum(b == 0) == length(b)){
    # Linear program coefficients with nonzero matrix A
    obj2 <- NULL
    obj1 <- A
    obj0 <- 0
  } else {
    # Quadratic program coefficients
    obj2 <- t(A) %*% A * n
    obj1 <- -2 * t(A) %*% b * n
    obj0 <- t(b) %*% b * n
  }
  # Return the above objective functions
  return(list(obj2 = obj2, obj1 = obj1, obj0 = obj0))
}

#' lpSolveAPI solver for linear programs
#'
#' @description This function computes the solution to the linear program
#'    using the `\code{lpsolveAPI}' package.
#'    
#' @inheritParams gurobi_optim
#' @inheritParams dkqs
#' @inheritParams prog_cone
#'
#' @returns Returns the optimal objective value and the corresponding argument
#'   to the linear program.
#'  \item{x}{Optimal point calculated from the linear program}
#'  \item{larg}{Arguments for the linear program.}
#'   
#' @details The package `\code{lpSolveAPI}' cannot be used to solve quadratic 
#'   programs. 
#'
#' @export
#' 
lpsolveapi.optim <- function(Af, bf, nf, A, rhs, sense, modelsense, lb){
  # ---------------- #
  # Step 1: Obtain the coefficients of the objective function
  # ---------------- #
  objective_return <- objective.function(Af, bf, nf)
  
  # ---------------- #
  # Step 2: Update the constraint matrices
  # ---------------- #
  # Change the lower bounds to inequality constriants
  lb_Amat <- diag(length(lb))
  lb_bvec <- lb
  # Update constraint matrices
  A <- rbind(A, lb_Amat)
  rhs <- c(rhs, lb_bvec)
  sense <- c(sense, rep(">=", length(lb_bvec)))
  
  # ---------------- #
  # Step 3: LP formulation
  # ---------------- #
  # solve object
  lprec <- make.lp(nrow = nrow(A), ncol = ncol(A))
  # Model sense
  lp.control(lprec, sense=modelsense)
  # Types of decision variables
  set.type(lprec, 1:ncol(A), type=c("real"))
  set.objfn(lprec, objective_return$obj1)
  #Define the constraints
  for (i in 1:nrow(A)){
    add.constraint(lprec, A[i, ], sense[i], rhs[i])
  }
  
  # ---------------- #
  # Step 4: Solve and obtain solution of LP
  # ---------------- #
  solve(lprec)
  x <- get.variables(lprec)
  objval <- get.objective(lprec)
  
  # ---------------- #
  # Step 5: Return results
  # ---------------- #
  return(list(objval = objval,
              x = x))
}