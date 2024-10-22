#' Conducts inference using the FSST procedure
#'
#' @description This module conducts inference in linear programs using the
#'   procedure by Fang, Santos, Shaikh and Torgovitsky (2020).
#'
#' @importFrom expm sqrtm
#' @importFrom Matrix t
#' @importFrom Matrix norm
#'
#' @inheritParams dkqs
#' @param data
#' @param lpmodel The \code{lpmodel} object.
#' @param lambda Parameter used to obtain the restricted estimator
#'   \eqn{\widehat{\bm{\beta}}^r_n}. A data-driven parameter \code{lambda} can
#'   be included if \code{NA} is included as part of the vector for
#'   \code{lambda}. For instance, if \code{lambda} is set as \code{c(0.1, NA)},
#'   then both 0.1 and the data-driven \code{lambda} will be applied in the
#'   \code{\link[lpinfer]{fsst}} test. The default is to use the data-driven
#'   \code{lambda}.
#' @param rho Parameter used in the studentization of matrices.
#' @param weight.matrix The option used in the weighting matrix. There are three
#'   options available: \itemize{ \item{\code{identity} --- identity matrix}
#'   \item{\code{diag} --- the diagonal matrix that takes the diagonal elements
#'   of the inverse of the variance matrix} \item{\code{avar} --- inverse of the
#'   variance matrix} }
#' @param beta.sigma.tol Tolerance level used to determine if the weight matrix
#'   is singular
#' @param beta.sigma.eps If weight.mat != "identity" and beta.sigma is
#'   determined to be singular, beta.sigma is replaced with beta.sigma +
#'   beta.sigma.eps * diag(nrow(beta.sigma))
#' @param sqrtm.method The method used to obtain the matrix square root in the
#'   \code{\link[lpinfer]{fsst}} procedure. This has to be a function that takes
#'   one argument that accepts a square matrix of size k x k and returns a
#'   square matrix of size k x k, where k can be the length of the
#'   \eqn{\beta(P)} vector, or the \code{beta.obs} component of the
#'   \code{lpinfer} object.
#' @param sqrtm.tol The absolute tolerance used to check whether the matrix
#'   square root is correct. This is done by checking whether the Frobenius norm
#'   is smaller than the tolerance level, i.e., when \eqn{A} is the give matrix,
#'   \eqn{B} is the matrix square root obtained from the given
#'   \code{sqrtm.method} function, and \eqn{\epsilon} is the tolerance level,
#'   the FSST test checks whether \eqn{||A - BB||_F < \epsilon}. If this does
#'   not hold, the FSST test will use the \code{\link[expm]{sqrtm}} function
#'   from the \code{expm} package to obtain the matrix square root.
#'
#' @return Returns the following information:
#'   \item{pval}{A table of \eqn{p}-values.}
#'   \item{cv.table}{A table of sample and bootstrap Cone and Range test
#'     statistics.}
#'   \item{call}{The matched call.}
#'   \item{range}{The sample range test statistic.}
#'   \item{cone}{The sample cone test statistic.}
#'   \item{test}{The sample test statistic.}
#'   \item{cone.n.list}{The list of bootstrap cone test statistics.}
#'   \item{range.n.list}{The list of bootstrap range test statistics.}
#'   \item{solver.name}{Name of the solver used.}
#'   \item{rho}{The value of \code{rho} provided by the user.}
#'   \item{rhobar.i}{The regularization parameter used for the Cone
#'     studentization matrix.}
#'   \item{lambda.data}{The value of the data-driven \code{lambda} (if
#'     applicable).}
#'   \item{var.method}{The method used in obtaining the asymptotic variance
#'     of \code{beta.obs}.}
#'   \item{test.logical}{An indicator variable for whether the computation has
#'     been conducted. If \code{test.logical} is 1, it refers to the case
#'     where \code{beta.tgt} is inside the logical bound. If
#'     \code{test.logical} is 0, it refers to the case where
#'     \code{beta.tgt} is outside the logical bound.}
#'   \item{logical.lb}{Logical lower bound.}
#'   \item{logical.ub}{Logical upper bound.}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'
#' @details The following components are required in the \code{lpmodel} for the
#'    \code{\link[lpinfer]{fsst}} procedure:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{A.shp}}
#'      \item{\code{beta.obs}}
#'      \item{\code{beta.shp}}
#'    }
#'
#' @section Example:
#' \preformatted{
#'   source("./example/dgp_missingdata.R") # Change directory if necessary
#'   J <- 5
#'   N <- 1000
#'   data <- missingdata_draw(J = J, n = N, seed = 1, prob.obs = .5)
#'   lpm <- missingdata_lpm(J = J, info = "full", data = data)
#'   fsst(data = data,
#'        lpmodel = lpm,
#'        beta.tgt = .2,
#'        R = 100,
#'        lambda = .2,
#'        rho = 1e-4,
#'        weight.matrix = "identity",
#'        solver = "gurobi")
#' }
#'
#' @section More examples:
#'   More examples can be found in the \code{fsst_example.R} file
#'   under the \code{example} subdirectory of the installation directory for
#'   the \code{lpinfer} package.
#'
#' @export
#'
fsst <- function(data = NULL, lpmodel, beta.tgt, R = 100, Rmulti = 1.25,
                 lambda = NA, rho = 1e-4, n = NULL, weight.matrix = "diag",
                 beta.sigma.tol = 1e-08, beta.sigma.eps = 1e-06,
                 solver = NULL, progress = TRUE,
                 sqrtm.method = function(m) pracma::sqrtm(m)$B,
                 sqrtm.tol = .Machine$double.eps^(1/2), previous.output = NA) {
   # ---------------- #
   # Step 1: Update call, check and update the arguments; initialize df.error
   # ---------------- #
   # Obtain call information
   call <- match.call()

   # Check the arguments
   fsst.return <- fsst.check(data, lpmodel, beta.tgt, R, Rmulti, lambda, rho,
                             n, weight.matrix, solver, progress,
                             sqrtm.method, sqrtm.tol, previous.output)

   # Update the arguments
   data <- fsst.return$data
   solver <- fsst.return$solver
   solver.name <- fsst.return$solver.name
   test.logical <- fsst.return$test.logical
   logical.lb <- fsst.return$logical.lb
   logical.ub <- fsst.return$logical.ub
   omega.i <- fsst.return$omega.i

   # Compute the maximum number of iterations
   maxR <- ceiling(R * Rmulti)

   # Initialize a table to contain the error messages
   df.error <- data.frame(matrix(vector(), ncol = 4))
   colnames(df.error) <- c("Iteration", "Step", "lambda", "Error message")

   ### Case 1: test.logical == 1. Proceed with the calculation because
   ### beta.tgt is inside the logical bounds
   if (test.logical == 1) {
      # The user must either provide the data or n
      if (is.null(data)) {
         n <- n
      } else {
         n <- nrow(data)
      }

      # Rearrange the lambda terms
      lambda.temp <- sort(lambda, decreasing = FALSE)

      # Check if NA was present in lambda
      if (NA %in% lambda) {
         lambda <- c(lambda.temp, NA)
      } else {
         lambda <- lambda.temp
      }

      # Define parameters and lists (-1 refers to the initial trial)
      R.succ <- -1
      i1 <- -1
      error.id0 <- NULL
      error.id <- NULL
      new.error.bs <- 0
      beta.obs.bs <- list()
      beta.n.bs <- list()

      # ---------------- #
      # Step 2: Obtain beta.obs, the list of bootstrap estimators and the
      # variance estimator
      # ---------------- #
      ### 2(a) Compute beta(P)
      beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)
      beta.obs.hat <- beta.obs.return[[1]]
      sigma.beta.obs <- beta.obs.return[[2]]
      beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)
      beta.n <- Reduce(rbind, c(unlist(beta.obs.hat), beta.shp.hat, beta.tgt))

      # Change maxR to the length of the list 'beta.obs' if it is a list
      if (inherits(lpmodel$beta.obs, "list")) {
         maxR <- length(lpmodel$beta.obs) - 1
      }

      # This while loop is used to re-draw the data if there are some
      # problematic draws in the bootstrap iterations
      while (((R.succ < R) & (i1 < maxR)) | (new.error.bs != 0)) {
         # Set the index set for the bootstrap replications
         if (R.succ == -1) {
            # -1 corresponds to the initial bootstrap replications
            i0 <- 1
            i1 <- R
            i1eqmax <- FALSE
            iseq <- 1:maxR
            eval.count <- 0
         } else {
            # Remove the problematic entries
            error.id.new <- error.id
            beta.obs.bs[error.id.new] <- NULL
            beta.n.bs[error.id.new] <- NULL

            # Save the current error.id
            error.id0 <- error.id

            # Update the sequence of indices
            i1eqmax <- (i1 == maxR)
            if (isFALSE(i1eqmax)) {
              i0 <- min(maxR, i1 + 1)
              i1 <- min(maxR, i0 + (R - R.succ) - 1)
            } else {
              i0 <- i1
              i1 <- i1
            }
            
            if (inherits(lpmodel$beta.obs, "list")) {
               i0 <- i1 + 1
               i1 <- i0 + (R - R.succ) - 1
            }
            iseq <- i0:i1
            eval.count <- eval.count + 1
         }

         ### 2(b) Estimate sigma.beta.obs and store the bootstrap estimates
         # If the user provided bootstrap estimates of beta, use it to compute
         # sigma
         if (inherits(lpmodel$beta.obs, "list")) {
            beta.obs.bs.new <- lpmodel$beta.obs[(i0 + 1):(i1 + 1)]
            if (!is.null(beta.obs.bs.new[[1]])) {
               beta.n.bs.new <- list()
               for (i in i0:i1) {
                  beta.n.bs.new[[i]] <- Reduce(rbind,
                                               c(beta.obs.bs.new[[i]],
                                                 beta.shp.hat,
                                                 beta.tgt))
               }
               beta.obs.bs <- c(beta.obs.bs, beta.obs.bs.new)
               beta.n.bs <- c(beta.n.bs, beta.n.bs.new)
            }
            var.method <- "list"
            if (is.null(sigma.beta.obs)) {
               sigma.beta.obs <- sigma.summation(n, lpmodel$beta.obs, progress,
                                                 eval.count)
               var.method <- "bootstrapped values of the input list"
            }
         } else {
            var.method <- "function"

            # Don't need to draw new obs if maxR is reached in last iteration
            if (isFALSE(i1eqmax)) {
              beta.obs.return <- fsst.beta.bs(n, data, beta.obs.hat, lpmodel,
                                              R, maxR, progress, df.error,
                                              iseq, eval.count)
              
              df.error <- beta.obs.return$df.error
              error.id <- beta.obs.return$error.id
              R.succ <- beta.obs.return$R.succ
              new.error.bs <- beta.obs.return$R.eval - R.succ
              if (new.error.bs != 0) {
                next
              }
              
              # Merge it with the beta.obs.bs that has been computed earlier
              beta.obs.bs.new <- beta.obs.return$beta.obs.bs
              beta.obs.bs <- c(beta.obs.bs, beta.obs.bs.new)
            }
            beta.obs.list <- c(list(beta.obs.hat), beta.obs.bs)

            # If sigma.beta is not provided by the function, compute it
            # from the bootstrap betas
            if (is.null(sigma.beta.obs)) {
               var.method <- paste0("bootstrapped 'beta.obs' ",
                                    "from the function")
               sigma.beta.obs <- sigma.summation(n, beta.obs.list, progress,
                                                 eval.count)
            }
            beta.n.bs <- full.beta.bs(lpmodel, beta.tgt, beta.obs.bs)
         }

         ### 2(c) Compute the beta.sigma
         n.beta1 <- nrow(sigma.beta.obs)
         n.beta23 <- length(beta.shp.hat) + length(beta.tgt)
         zero.12 <- matrix(rep(0, n.beta1 * n.beta23), nrow = n.beta1)
         zero.21 <- Matrix::t(zero.12)
         zero.22 <- matrix(rep(0, n.beta23^2), nrow = n.beta23)
         beta.sigma <- rbind(cbind(sigma.beta.obs, zero.12),
                             cbind(zero.21, zero.22))

         # ---------------- #
         # Step 3: Estimate beta.star, x.star and their bootstrap counterparts
         # ---------------- #
         p <- length(beta.n)
         if (!is.matrix(lpmodel$A.tgt)) {
            d <- length(lpmodel$A.tgt)
         } else {
            d <- ncol(lpmodel$A.tgt)
         }

         # Compute the weighting matrix
         weight.mat <- fsst.weight.matrix(weight.matrix,
                                          beta.obs.hat,
                                          sigma.beta.obs,
                                          beta.sigma.tol,
                                          beta.sigma.eps)

         # Compute the matrix square root of the weighting matrix
         weight.mat.root <- checkupdate.matrixroot(weight.mat,
                                                   "weighting matrix",
                                                   sqrtm.method,
                                                   sqrtm.tol)

         # Compute beta.star
         beta.star.return <- fsst.beta.star.bs(data, lpmodel, beta.n,
                                               beta.n.bs, beta.tgt, weight.mat,
                                               beta.obs.hat, beta.obs.bs, R,
                                               sigma.beta.obs, solver,
                                               df.error, p, d, progress,
                                               eval.count)
         beta.star <- beta.star.return$beta.star
         beta.star.bs <- beta.star.return$beta.star.bs
         x.star <- beta.star.return$x.star
         x.star.bs <- beta.star.return$x.star.bs

         # Update the list of errors and restart the loop if necessary
         df.error <- beta.star.return$df.error
         new.error.bs <- beta.star.return$new.error
         error.id <- beta.star.return$error.id
         R.succ <- length(beta.star.bs)
         if (new.error.bs != 0) {
            next
         }

         # Consolidate the bootstrap estimators into a list
         beta.star.l <- list(beta.star = beta.star)
         beta.star.list <- c(beta.star.l, beta.star.bs)

         # ---------------- #
         # Step 4: Studentization
         # ---------------- #
         # Obtain the star version of the sigma matrix
         if (d >= p) {
            sigma.star <- beta.sigma
         } else {
            sigma.star <- sigma.summation(n, beta.star.list, progress,
                                          eval.count)
         }

         # Compute the matrix square root
         rhobar.i <- Matrix::norm(sigma.star, type = "f") * rho

         # Compute the studentization matrix if 'omega.i' is NA and if d < p
         if (!is.matrix(omega.i) | d < p) {
            sigma.reg <- sigma.star + rhobar.i * diag(nrow(sigma.star))
            omega.i <- checkupdate.matrixroot(sigma.reg,
                                              "studentization matrix",
                                              sqrtm.method,
                                              sqrtm.tol)
         }

         # ---------------- #
         # Step 5: Test statistic
         # ---------------- #
         # Compute range.n
         if (d >= p) {
            range.n <- 0
            cone.n <- fsst.cone.lp(n, omega.i, beta.n, beta.star, lpmodel, 1,
                                   solver)
         } else {
            range.n <- fsst.range(n, beta.obs.hat, x.star, lpmodel,
                                  weight.mat.root)
            cone.n <- fsst.cone.lp(n, omega.i, beta.n, beta.star, lpmodel, 0,
                                   solver)
         }

         # ---------------- #
         # Step 6: Choosing a data-driven lambda (if applicable)
         # ---------------- #
         if (NA %in% lambda) {
            lambda.data <- fsst.lambda(n, omega.i, beta.n, beta.star, lpmodel,
                                       R.succ, beta.star.list, solver, progress,
                                       df.error, p, d, eval.count)
            lambda.dd <- lambda.data$lambda
            new.error.bs <- lambda.data$new.error
            df.error <- lambda.data$df.error
            R.succ <- lambda.data$R.succ
            error.id <- lambda.data$error.id
            if (new.error.bs != 0) {
               next
            } else {
               # Append the data-driven data if there is no error
               lambda <- c(lambda, lambda.dd)
            }
         } else {
            lambda.dd <- NULL
         }

         # Drop the NA term from lambda
         lambda <- lambda[!is.na(lambda)]

         # ---------------- #
         # Step 7: Compute bootstrap components of cone.n and range.n
         # ---------------- #
         if (d >= p) {
            # Compute the restricted estimator
            beta.r <- beta.r.compute(n, lpmodel, beta.obs.hat, beta.tgt,
                                     beta.n, beta.star, omega.i, 1, solver)$x

            # Compute range.n for bootstrap beta
            range.n.list <- rep(0, R.succ)

            # Compute cone.n for bootstrap beta
            cone.n.list <- list()
            for (i in 1:length(lambda)) {
               cone.return <- fsst.cone.bs(n, omega.i, beta.n, beta.star,
                                           lpmodel, R.succ, lambda[i], 1,
                                           beta.r, beta.star.list, solver,
                                           progress, df.error, eval.count,
                                           FALSE)
               cone.n.list[[i]] <- cone.return$cone.n.list

               # Update the list of errors and restart the loop if necessary
               df.error <- cone.return$df.error
               error.id <- cone.return$error.id
               new.error.bs <- cone.return$new.error
               R.succ <- length(cone.return$cone.n.list)
               if (new.error.bs != 0) {
                  next
               }
            }

         } else {
            # Compute the restricted estimator
            beta.r <- beta.r.compute(n, lpmodel, beta.obs.hat, beta.tgt,
                                     beta.n, beta.star, omega.i, 0, solver)$x

            # print(beta.r)
            
            # Compute range.n for bootstrap beta
            range.return <- fsst.range.bs(n, lpmodel, beta.obs.hat,
                                          beta.obs.bs, x.star, x.star.bs,
                                          weight.mat.root, R, progress,
                                          df.error, eval.count)

            range.n.list <- range.return$range.n.list

            # Update the list of errors and restart the loop if necessary
            df.error <- range.return$df.error
            new.error.bs <- range.return$new.error
            R.succ <- length(range.n.list)
            error.id <- range.return$error.id

            if (new.error.bs != 0) {
               break()
            }

            # Compute cone.n for bootstrap beta
            cone.n.list <- list()
            for (i in 1:length(lambda)) {
               cone.return <- fsst.cone.bs(n, omega.i, beta.n, beta.star,
                                           lpmodel, R.succ, lambda[i], 0,
                                           beta.r, beta.star.list, solver,
                                           progress, df.error, eval.count,
                                           FALSE)
               cone.n.list[[i]] <- cone.return$cone.n.list
               R.succ <- length(cone.return$cone.n.list)

               # Update the list of errors and break the for-loop if there is
               # any errors
               df.error <- cone.return$df.error
               new.error.bs <- cone.return$new.error
               error.id <- cone.return$error.id
               if (new.error.bs != 0) {
                  break()
               }
               new.error.bs <- 0
            }
         }

         # Restart the loop if necessary
         if (new.error.bs != 0) {
            next
         }
      }

      # Parameters
      n.lambda <- length(lambda)

      # Initialize the data frames
      df.pval <- data.frame(matrix(vector(), nrow = n.lambda, ncol = 2))
      colnames(df.pval) <- c("lambda", "p-value")
      df.pval$lambda <- lambda

      if (R.succ != 0) {
         # ---------------- #
         # Step 8: Compute decision, p-value and the quantiles of the test
         # statistics
         # ---------------- #
         T.bs <- list()
         for (i in 1:n.lambda) {
            # Compute the p-values
            pval.return <- fsst.pval(range.n, cone.n$objval, range.n.list,
                                     cone.n.list[[i]], R.succ)
            df.pval[i, 2] <- pval.return$pval
            T.bs[[i]] <- pval.return$T.bs
         }

         # Compute cv.table
         cv.table <- fsst.cv.table(lambda, "lambda",
                                   rep(cone.n$objval, n.lambda),
                                   range.n, cone.n.list, range.n.list, T.bs)

         # ---------------- #
         # Step 9: Assign the return list and return output
         # ---------------- #
         # Assign the list of objects returned
         output <- list(pval = df.pval,
                        cv.table = cv.table,
                        range = range.n,
                        cone = cone.n$objval,
                        test = max(range.n, cone.n$objval),
                        cone.n.list = cone.n.list,
                        range.n.list = range.n.list,
                        rhobar.i = rhobar.i,
                        lambda.data = lambda.dd,
                        var.method = var.method,
                        omega.i = omega.i,
                        R.succ = R.succ,
                        df.error = df.error)
      } else {
         df.pval[, 2] <- NA
         output <- list(pval = df.pval,
                        cv.table = NA,
                        rhobar.i = NA,
                        R.succ = R.succ,
                        df.error = df.error)
      }
   } else {
      ### Case 2: test.logical == 0. Set the p-value as 0 directly because
      ### beta.tgt is outside the logical bounds
      output <- list(pval = 0)

      # Print warning message
      infeasible.betatgt.warning()
   }

   # Assign the common objects in the output list
   output <- append(output,
                    list(solver = solver.name,
                         call = call,
                         rho = rho,
                         test.logical = test.logical,
                         logical.lb = logical.lb,
                         logical.ub = logical.ub))

   # Turn returned lists for bounds as vectors if parameter is not multivalued
   if (length(output$cone.n.list) == 1) {
      output$cone.n.list <- unlist(output$cone.n.list)
   }

   # Return NULL df.error if no error
   if (length(nrow(output$df.error)) == 0) {
      output$df.error <- NULL
   } else if (nrow(output$df.error) == 0) {
      output$df.error <- NULL
   }

   # Assign class
   attr(output, "class") <- "fsst"

   # Return output
   return(output)
}

#' Construct the full beta vector in the \code{\link[lpinfer]{fsst}} procedure
#'
#' @description This function concatenate the \code{beta.obs}, \code{beta.shp}
#'   and \code{beta.tgt} components to form the full \code{beta} vector. In
#'   particular, the full vector is defined as
#'   \eqn{\widehat{\bm{\beta}}_n} vector that is defined as
#'   \eqn{\widehat{\bm{\beta}}_n \equiv (\widehat{\bm{\beta}}_{{\rm obs},n},
#'   \bm{\beta}_{{\rm shp},n}, \bm{\beta}_{{\rm tgt}})'}.
#'
#' @inheritParams fsst
#' @param beta.obs.bs The bootstrap estimates of \code{beta.obs}.
#'
#' @return List of full beta vectors
#'   \item{beta.bs}{The list of full beta vector
#'   \eqn{\{\widehat{\bm{\beta}}_{n,b}\}^B_{b=1}}.}
#'
#' @export
#'
full.beta.bs <- function(lpmodel, beta.tgt, beta.obs.bs) {
   beta.bs <- list()
   for (i in 1:length(beta.obs.bs)) {
      beta.bs[[i]] <- Reduce(rbind,
                             c(beta.obs.bs[[i]], lpmodel$beta.shp, beta.tgt))
   }

   return(beta.bs)
}

#' Computing the bootstrap estimates of \code{beta.obs}
#'
#' @description This function computes the bootstrap estimates of
#'   \eqn{\hat{\beta}_{{\rm obs}, n}}.
#'
#' @import furrr progressr
#' @importFrom utils tail
#'
#' @inheritParams dkqs.bs
#' @inheritParams dkqs.bs.fn
#' @inheritParams fsst
#' @param iseq The list of indices or betas to iterate over.
#' @param beta.obs.hat The sample estimator
#'   \eqn{\widehat{\bm{\beta}}_{\mathrm{obs}, n}} based on the given
#'   information in \code{lpmodel} (and \code{data} if applicable).
#' @param df.error A table showing the id of the bootstrap replication(s)
#'   with error(s) and the corresponding error message(s).
#'
#' @return Returns the bootstrap estimators.
#'   \item{beta.obs.bs}{A list of bootstrap estimators
#'   \eqn{\{\hat{\beta}_{{\rm obs}, n, b}\}^B_{b=1}}.}
#'   \item{df.error}{An updated table showing the id of the bootstrap
#'     replication(s) with error(s) and the corresponding error message(s).}
#'   \item{R.eval}{The number of bootstrap replications that has been
#'     conducted.}
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'
#' @export
#'
fsst.beta.bs <- function(n, data, beta.obs.hat, lpmodel, R, maxR, progress,
                         df.error, iseq, eval.count) {

   # ---------------- #
   # Step 1: Initialization
   # ---------------- #
   R.succ <- 0
   R.eval <- 0
   beta.obs.bs <- list()
   error.list <- list()

   # Check if there is any list objects in 'lpmodel'
   any.list <- lpmodel.anylist(lpmodel)

   # Assign the number of replications required
   if (identical(iseq, 1:maxR)) {
      Rcomp <- R
   } else {
      Rcomp <- length(iseq)
   }
   # ---------------- #
   # Step 2: Bootstrap replications
   # ---------------- #
   # This function will not be called if 'beta.obs' is a list. Hence, it
   # suffices to consider the case where it is a function.
   while ((R.succ < Rcomp) & (R.eval != maxR)) {
      # Compute the list of indices to be passed to 'future_map'
      if (identical(iseq, 1:maxR)) {
         # Evaluate the list of indices to be passed to 'future_map'
         bs.ind <- bs.index(R, R.eval, R.succ, maxR)
         i0 <- bs.ind$i0
         i1 <- bs.ind$i1
      } else {
         bs.list <- iseq
         i0 <- bs.list[1]
         i1 <- utils::tail(bs.list, n = 1)
      }

      # Set the default for progress bar
      progressr::handlers("progress")

      # Obtain results from the bootstrap replications
      progressr::with_progress({
         if (isTRUE(progress)) {
            pbar <- progressr::progressor(along = i0:i1)
         } else {
            pbar <- NULL
         }

         beta.obs.return <- furrr::future_map(i0:i1,
                                              .f = fsst.beta.bs.fn,
                                              data = data,
                                              lpmodel = lpmodel,
                                              pbar = pbar,
                                              progress = progress,
                                              eval.count = eval.count,
                                              n.bs = i1 - i0 + 1,
                                              any.list = any.list,
                                              .options =
                                                 furrr::furrr_options(seed = TRUE))
         eval.count <- eval.count + 1
      })

      # Update the list and parameters
      post.return <- post.bs(beta.obs.return, i0, i1, R.eval, T.list = NULL,
                             beta.list = beta.obs.bs, error.list = error.list)
      beta.obs.bs <- post.return$beta.list
      error.list <- post.return$error.list
      R.succ <- post.return$R.succ
      R.eval <- post.return$R.eval
   }

   # ---------------- #
   # Step 3: Consolidate the error messages
   # ---------------- #
   if (R.eval != R.succ) {
      # Update df.error if it is coming from this step
      if (length(unlist(error.list)) != 0) {
         # Create data.frame for error messages
         df.error <- data.frame(id = NA,
                                lambda = NA,
                                message = unlist(error.list))

         # Match the id of the error messages
         df.error <- error.id.match(error.list, df.error)
      }
   }

   return(list(beta.obs.bs = beta.obs.bs,
               df.error = df.error,
               R.eval = R.eval,
               R.succ = R.succ))
}

#' Computes one bootstrap estimates for \code{beta.obs}.
#'
#' @description This function carries out one bootstrap replication for
#'   getting \code{beta.obs} in the \code{\link[lpinfer]{fsst}} procedure.
#'   This function is used in the \code{\link[lpinfer]{fsst.beta.bs}} function
#'   via the \code{future_map} function.
#'
#' @inheritParams dkqs.bs.fn
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @param x This is either the list of indices that represent the bootstrap
#'   replications, or the list of bootstrap components of the \code{lpmodel}
#'   object passed from the user.
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{beta}{A bootstrap estimator of \eqn{\hat{\beta}_{n,b}}.}
#'   \item{msg}{An error message (if applicable).}
#'
#' @export
#'
fsst.beta.bs.fn <- function(x, data, lpmodel, pbar, progress, eval.count,
                            n.bs, any.list) {
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
   # Step 2: Initialize lpmodel
   # ---------------- #
   # Replace lpmodel by x if x is a list
   lpm <- lpmodel
   if (!is.null(any.list)) {
      for (nm in any.list$name) {
         lpm[[nm]] <- lpmodel[[nm]][[x + 1]]
      }
   }

   # ---------------- #
   # Step 3: Conduct one bootstrap replication
   # ---------------- #
   # Draw data
   data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
   rownames(data.bs) <- 1:nrow(data.bs)

   # Bootstrap estimator
   result <- tryCatch({
      beta.obs.return <- lpmodel.beta.eval(data.bs, lpm$beta.obs, 1)[[1]]
      list(beta.obs.return = beta.obs.return)
   }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
   }, error = function(e) {
      return(list(status = "error",
                  msg = e))
   })

   # Consolidate the information to be returned
   if (is.null(result$status)) {
      beta <- result$beta
      msg <- NULL
   } else {
      beta <- NULL
      msg <- result$msg$message
   }

   return(list(beta = beta,
               msg = msg))
}

#' Computes the weighting matrix in the \code{\link[lpinfer]{fsst}} procedure
#'
#' @description This function returns the weighting matrix in the
#'   \code{\link[lpinfer]{fsst}} procedure.
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @param beta.sigma The variance estimator of \eqn{\bm{\beta}_{\rm obs}}.
#'
#' @return The weighting matrix is returned.
#'    \item{weight.max}{The weighting matrix.}
#'
#' @export
#'
fsst.weight.matrix <- function(weight.matrix, beta.obs.hat, beta.sigma,
                               beta.sigma.tol, beta.sigma.eps) {
   # ---------------- #
   # Step 1: Convert the string to lower case.
   # ---------------- #
   weight.matrix <- tolower(weight.matrix)

   # ---------------- #
   # Step 2: Create the matrix
  # ---------------- #
  if (weight.matrix == "identity") {
    weight.mat <- diag(nrow(asmat(beta.obs.hat)))
  } else if (weight.matrix == "diag") {
    if (min(eigen(beta.sigma)$values) < beta.sigma.tol) {
      warning("beta.sigma is singular. A small identity matrix has been added to beta.sigma to calculate your selected weight.matrix.")
      beta.sigma <- beta.sigma + beta.sigma.eps * diag(nrow(beta.sigma))
    }
    weight.mat <- diag(diag(solve(beta.sigma)))     
  } else if (weight.matrix == "diag2") {
    weight.mat <- 1/diag(diag(beta.sigma))
  } else if (weight.matrix == "avar") {
    if (min(eigen(beta.sigma)$values) < beta.sigma.tol) {
      warning("beta.sigma is singular. A small identity matrix has been added to beta.sigma to calculate your selected weight.matrix.")
      beta.sigma <- beta.sigma + beta.sigma.eps * diag(nrow(beta.sigma))
    }
    weight.mat <- solve(beta.sigma)
  } else {
    stop("'weight.matrix' has to be one of 'avar', 'diag', 'diag2', and 'identity'.")
  }

   return(weight.mat)
}

#' Computes the asymptotic variance estimator
#'
#' @import furrr progressr
#'
#' @description Based on the bootstrap estimates
#'   \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}, this function computes the
#'   asymptotic variance estimator of the bootstrap estimator, i.e.
#'   \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b -
#'   \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b -
#'   \widehat{\bm{\beta}}\right)'.}
#'   This function supports parallel programming via the \code{furrr}
#'   package.
#'
#' @inheritParams dkqs
#' @inheritParams dkqs.bs.fn
#' @param n Sample size.
#' @param beta.bs.list A list of bootstrap estimators
#'    \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.
#'
#' @return Returns the estimator of the asymptotic variance.
#'     \item{sigma.mat}{The estimator of the asymptotic variance.}
#'
#' @usage sigma.summation(n, beta.bs.list, progress, eval.count)
#'
#' @export sigma.summation
#'
sigma.summation <- function(n, beta.bs.list, progress, eval.count) {
   beta.obs.hat <- beta.bs.list[[1]]

   progressr::with_progress({
      if (isTRUE(progress)) {
         pbar <- progressr::progressor(along = seq_along(beta.bs.list[-1]))
      } else {
         pbar <- NULL
      }

      beta.prod.return <- furrr::future_map(beta.bs.list[-1],
                                            .f = beta.product,
                                            beta.obs.hat = beta.obs.hat,
                                            pbar = pbar,
                                            progress = progress,
                                            eval.count = eval.count)
   })

   sigma.mat <- Reduce("+", beta.prod.return) * n / (length(beta.bs.list) - 1)

   return(sigma.mat)
}

#' Computes vector products
#'
#' @description This function computes the product of the two vectors. This is
#'   used in the \code{\link[lpinfer]{sigma.summation}} function that computes
#'   the asymptotic variance estimator.
#'
#' @importFrom Matrix t
#'
#' @details Denote \eqn{\bm{\beta}} and \eqn{\hat{\bm{\beta}}_{\rm obs}} as
#'   the \eqn{n \times 1} vectors \code{beta} and \code{beta.obs.hat}
#'   respectively. This function returns the \eqn{n \times n} matrix
#'   \eqn{\widetilde{\bm{\beta}}\widetilde{\bm{\beta}}'} where
#'   \eqn{\widetilde{\bm{\beta}}} is defined as
#'   \eqn{\widetilde{\bm{\beta}} \equiv \bm{\beta} -
#'   \widehat{\bm{\beta}}_{\rm obs}}.
#'
#' @param beta A bootstrap estimator \code{beta.obs}.
#' @inheritParams fsst.beta.bs
#' @inheritParams dkqs.bs.fn
#'
#' @return Returns an \eqn{n \times n} matrix.
#'     \item{beta.prod}{An \eqn{n \times n} matrix.}
#'
#' @export
#'
beta.product <- function(beta, beta.obs.hat, pbar, progress, eval.count) {
   # ---------------- #
   # Step 1: Print progress bar
   # ---------------- #
   if (eval.count == 0) {
      compute <- "Computing"
   } else {
      compute <- "Recomputing"
   }

   if (isTRUE(progress)) {
      pbar(sprintf("(%s asymptotic variance estimator)", compute))
   }

   beta.diff <- asmat(beta - beta.obs.hat)
   if (nrow(beta.diff) == 1) {
      beta.prod <- Matrix::t(beta.diff) %*% beta.diff
   } else {
      beta.prod <- beta.diff %*% Matrix::t(beta.diff)
   }

   return(beta.prod)
}

#' Computes the starred components of \eqn{\widehat{\bm{\beta}}}
#'
#' @description This function computes the vector
#'   \eqn{\widehat{\bm{\beta}}_n^\star} for the case where \eqn{d<p} in the
#'   \code{\link[lpinfer]{fsst}} procedure.
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @inheritParams fsst.weight.matrix
#' @inheritParams sigma.summation
#' @param weight.mat The weighting matrix for the \code{\link[lpinfer]{fsst}}
#'   procedure.
#'
#' @details This corresponding to solving the following quadratic program
#' \deqn{
#'   \min_{x \in \mathbf{R}^d} \,\,
#'   \left(\hat{\beta}_{{\rm obs}, n} - A_{\rm obs} x\right)' \hat{\Xi}
#'   \left(\hat{\beta}_{{\rm obs}, n} - A_{\rm obs} x\right)
#'   \quad \mathrm{s.t.
#'   } \quad
#'   A_{\rm shp} x = \beta_{\rm shp}
#'   \quad \mathrm{ and } \quad
#'   A_{\rm tgt} x = \beta_{\rm tgt}
#' }
#' in the \code{\link[lpinfer]{fsst}} procedure.
#'
#' @return Returns the following objects:
#'   \item{beta.star}{The vector \eqn{\widehat{\bm{\beta}}_n^\star}.}
#'   \item{x}{The optimal point.}
#'   \item{status}{The status of the optimization problem.}
#'
#' @export
#'
beta.star.qp <- function(data, lpmodel, beta.tgt, weight.mat, beta.obs.hat,
                         beta.sigma, solver) {
   # ---------------- #
   # Step 1: Solve the quadratic program
   # ---------------- #
   # Define the A matrices
   A.obs.hat <- lpmodel.eval(data, lpmodel$A.obs, 1)
   A.shp.hat <- lpmodel.eval(data, lpmodel$A.shp, 1)
   A.tgt.hat <- lpmodel.eval(data, lpmodel$A.tgt, 1)

   # Define the parameters
   d <- ncol(A.tgt.hat)

   # Constraints matrix
   A.mat1 <- A.shp.hat
   A.mat2 <- A.tgt.hat
   A.mat <- asmat(rbind(A.mat1, A.mat2))

   # RHS matrix
   b <- matrix(c(lpmodel$beta.shp, c(beta.tgt)), ncol = 1)

   # Objective function
   A.obj <- A.obs.hat

   # Assign the optimization argument
   optim.arg <- list(Af = A.obj,
                     bf = beta.obs.hat,
                     nf = 1,
                     A = A.mat,
                     rhs = b,
                     sense = "=",
                     modelsense = "min",
                     lb = rep(-Inf, d),
                     weight = weight.mat,
                     Method = 0)

   # Solve the model
   ans <- do.call(solver, optim.arg)

   # ---------------- #
   # Step 2: Compute beta.star
   # ---------------- #
   # Compute x.star
   x.star <- ans$x

   # Compute beta.star
   A <- asmat(rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt))
   beta.star <- A %*% x.star

   return(list(beta.star = beta.star,
               x.star = x.star,
               status = ans$status))
}

#' Computes the solution to the cone problem
#'
#' @importFrom methods as
#'
#' @description This function computes the solution to the cone problem.
#'
#' @importFrom Matrix t
#' @importFrom Matrix Matrix
#'
#' @inheritParams fsst
#' @inheritParams fsst.cone.bs
#' @inheritParams fsst.beta.star.bs
#' @param omega.i The matrix \eqn{\widehat{\bm{\Omega}}^i_n}, i.e. the
#'   regularized matrix for
#'   \eqn{\widehat{\bm{\Sigma}}^{\beta^\star}_{n,\bar{\rho}}}.
#' @param indicator A binary variable that equals to 1 for \eqn{d \geq p} and
#'   equals to 0 for \eqn{d < p}.
#'
#' @return Returns the optimal point and optimal value.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'
#' @export
#'
fsst.cone.lp <- function(n, omega.i, beta.n, beta.star, lpmodel, indicator,
                         solver) {
   # ---------------- #
   # Step 1: Construct the linear program
   # ---------------- #
   p <- length(beta.n)
   d <- ncol(lpmodel$A.obs)
   ones.p <- Matrix::Matrix(rep(1, p), nrow = 1, ncol = p, sparse = TRUE)
   zero.p <- Matrix::Matrix(data = 0, nrow = 1, ncol = p, sparse = TRUE)
   zero.p1 <- Matrix::Matrix(data = 0, nrow = p, ncol = 1, sparse = TRUE)
   zero.d <- Matrix::Matrix(data = 0, nrow = 1, ncol = d, sparse = TRUE)
   zero.d1 <- Matrix::Matrix(data = 0, nrow = d, ncol = 1, sparse = TRUE)
   zero.dp <- Matrix::Matrix(data = 0, nrow = d, ncol = p, sparse = TRUE)

   # Update the objective function
   obj <- rbind(beta.star, zero.p1, zero.p1)

   # Construct the lower bound
   lb <- c(rep(-Inf, p), rep(0, 2 * p))

   # Construct the constraints matrix
   A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
   A.mat1 <- methods::as(cbind(omega.i, -diag(p), diag(p)), "sparseMatrix")
   A.mat2 <- cbind(zero.p, ones.p, ones.p)
   A.mat3 <- methods::as(cbind(Matrix::t(A), zero.dp, zero.dp), "sparseMatrix")
   A.mat <- rbind(A.mat1, A.mat2, A.mat3)

   # Construct RHS vector
   rhs.mat <- Reduce(rbind, c(zero.p1, 1, zero.d1))

   # Sense
   sense.mat <- c(rep("=", p), rep("<=", d + 1))

   # Construct the arguments
   if (indicator == 1) {
      optim.arg <- list(Af = NULL,
                        bf = obj,
                        nf = 1,
                        A = A.mat,
                        rhs = rhs.mat,
                        sense = sense.mat,
                        modelsense = "max",
                        lb = lb)
   } else {
      zero.pp <- Matrix::Matrix(rep(0, p * p), ncol = p, nrow = p,
                                sparse = TRUE)
      zero.Am <- Matrix::Matrix(rep(0, (p + d + 1) * d),
                                nrow = p + d + 1,
                                ncol = d,
                                sparse = TRUE)

      # Update objective function
      obj.ext <- rbind(obj, zero.d1)

      # Update constraints matrix
      A.mat.ext1 <- asmat(cbind(A.mat, zero.Am))
      A.mat.ext2 <- asmat(cbind(diag(p), zero.pp, zero.pp, -A))
      A.mat.ext <- methods::as(rbind(A.mat.ext1, A.mat.ext2), "sparseMatrix")

      # Update RHS vector
      rhs.ext <- Reduce(rbind, c(rhs.mat, zero.p1))
      sense.ext <- c(sense.mat, rep("=", p))

      # Update lower bound
      lb.ext <- c(lb, rep(-Inf, d))

      # Set the arguments
      optim.arg <- list(Af = NULL,
                        bf = obj.ext,
                        nf = 1,
                        A = A.mat.ext,
                        rhs = rhs.ext,
                        sense = sense.ext,
                        modelsense = "max",
                        lb = lb.ext)
   }

   # ---------------- #
   # Step 2: Solve the linear program
   # ---------------- #
   ans <- do.call(solver, optim.arg)
   objval <- ans$objval * sqrt(n)

   return(list(x = ans$x,
               objval = objval,
               status = ans$status))
}

#' Computes the bootstrap estimates of the starred version of \code{beta.obs}
#'
#' @description This function computes the bootstrap estimates of
#'   \eqn{\widehat{\bm{\beta}}^\star_n}.
#'
#' @import furrr progressr
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @inheritParams beta.star.qp
#' @inheritParams full.beta.bs
#' @inheritParams fsst.weight.matrix
#' @param beta.n The sample \eqn{\widehat{\bm{\beta}}_n} vector that is defined
#'   as \eqn{\widehat{\bm{\beta}}_n \equiv (\widehat{\bm{\beta}}_{{\rm obs},n},
#'   \bm{\beta}_{{\rm shp},n}, \bm{\beta}_{{\rm tgt}})'}.
#' @param beta.n.bs The bootstrap estimates of \code{beta.n}.
#' @param sigma.beta.obs An estimator of the asymptotic variance for
#'   \code{beta.obs}.
#' @param p The length of the beta vector.
#' @param d The number of columns of the \code{A} matrix.
#'
#' @return Return the following list of objects:
#'   \item{beta.star}{This corresponds to \eqn{\widehat{\bm{\beta}}^\star}.}
#'   \item{beta.star.bs}{This corresponds to the bootstrap estimates of
#'     \eqn{\widehat{\bm{\beta}}^\star}.}
#'   \item{x.star}{This corresponds to \eqn{\bm{x}^\star}, which is the
#'     solution to the quadratic program \code{\link[lpinfer]{beta.star.qp}}.}
#'   \item{x.star.bs}{This corresponds to the bootstrap estimates of
#'     \eqn{\bm{x}^\star}.}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{new.error}{The number of new errors.}
#'   \item{error.id}{The indices of the bootstrap replications that give
#'     errors.}
#'
#' @export
#'
fsst.beta.star.bs <- function(data, lpmodel, beta.n, beta.n.bs, beta.tgt,
                              weight.mat, beta.obs.hat, beta.obs.bs, R,
                              sigma.beta.obs, solver, df.error, p, d,
                              progress, eval.count) {
   # ---------------- #
   # Step 1: Initialization
   # ---------------- #
   # Original number of problematic bootstrap draws
   error0 <- nrow(df.error)

   # List of error indices
   error.id <- NULL

   # ---------------- #
   # Step 2: Compute 'beta.star', 'x.star' and their bootstrap estimates
   # ---------------- #
   if (d >= p) {
      # Case 1: d >= p - All star components are the same as the non-star
      # components and the x.star objects are not used in the latter parts
      beta.star <- beta.n
      beta.star.bs <- beta.n.bs
      x.star <- NULL
      x.star.bs <- NULL

      # There is no error in simply defining new objects
      df.error <- df.error
      new.error <- 0
   } else {
      # Case 2: d < p - Need to compute beta.star, x.star and their bootstrap
      # counterparts
      # Solve the quadratic program to get beta.star and x.star
      qp.return <- beta.star.qp(data, lpmodel, beta.tgt, weight.mat,
                                beta.obs.hat, sigma.beta.obs, solver)
      beta.star <- qp.return$beta.star
      x.star <- qp.return$x.star

      # Construct bootstrap estimates of beta.star
      beta.star.bs <- list()
      x.star.bs <- list()
      n.bs <- length(beta.obs.bs)

      progressr::with_progress({
         if (isTRUE(progress)) {
            pbar <- progressr::progressor(along = 1:n.bs)
         } else {
            pbar <- NULL
         }

         fsst.star.return <-
            furrr::future_map(beta.obs.bs,
                              .f = fsst.beta.star.bs.fn,
                              data = data,
                              lpmodel = lpmodel,
                              beta.tgt = beta.tgt,
                              weight.mat = weight.mat,
                              sigma.beta.obs = sigma.beta.obs,
                              solver = solver,
                              pbar = pbar,
                              progress = progress,
                              n.bs = n.bs,
                              eval.count = eval.count)
      })

      # Get beta.star and x.star and exclude those with 'NULL'
      beta.star.bs <- sapply(fsst.star.return, "[", "beta.star")
      beta.star.bs <- beta.star.bs[lengths(beta.star.bs) != 0]
      x.star.bs <- sapply(fsst.star.return, "[", "x.star")
      x.star.bs <- x.star.bs[lengths(x.star.bs) != 0]

      # Retrieve error list
      error.list <- sapply(fsst.star.return, "[", "msg")

      # Evaluate R.eval and R.succ
      R.eval <- length(beta.obs.bs)
      R.succ <- length(beta.star.bs)

      # ---------------- #
      # Step 3: Consolidate the error messages
      # ---------------- #
      if (R.eval != R.succ) {
         # Create data.frame for error messages
         df.error1 <- data.frame(id = NA,
                                 lambda = NA,
                                 message = unlist(error.list))

         # Match the id of the error messages
         df.error1 <- error.id.match(error.list, df.error1)
         if (!is.null(df.error)) {
            # Case 1: There are errors in this function and the earlier parts
            error0 <- nrow(df.error)
            df.error <- rbind(df.error, df.error1)
            rownames(df.error) <- 1:nrow(df.error)
         } else {
            # Case 2: There are no errors in the earlier steps and there are
            # errors in this step
            df.error <- df.error1
         }
         new.error <- R.eval - R.succ
         error.id <- df.error1$id
      } else {
         # Case 3: There are no new errors in this procedure
         df.error <- df.error
         new.error <- 0
         error.id <- NULL
      }
   }

   return(list(beta.star = beta.star,
               beta.star.bs = beta.star.bs,
               x.star = x.star,
               x.star.bs = x.star.bs,
               df.error = df.error,
               new.error = new.error,
               error.id = error.id))
}

#' Computes one bootstrap estimates of \code{beta.star} and \code{x.star}
#'
#' @description This function computes one bootstrap estimate for
#'   \code{beta.star} and \code{x.star}. This function is used
#'   in the \code{\link[lpinfer]{fsst.beta.star.bs}} function via the
#'   \code{future_map} command.
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.star.bs
#' @inheritParams full.beta.bs
#' @inheritParams dkqs.bs.fn
#'
#' @return Returns a list of output that are obtained from the subsampling
#'   procedure:
#'   \item{beta.star}{The bootstrap estimates of \code{beta.star}.}
#'   \item{x.star}{The bootstrap estimates of \code{x.star}.}
#'   \item{msg}{An error message (if applicable).}
#'
#' @export
#'
fsst.beta.star.bs.fn <- function(beta.obs.bs, data, lpmodel, beta.tgt,
                                 weight.mat, sigma.beta.obs, solver,
                                 pbar, progress, n.bs, eval.count) {
   # ---------------- #
   # Step 1: Print progress bar
   # ---------------- #
   if (eval.count == 0) {
      compute <- "Computing"
   } else {
      compute <- "Recomputing"
   }

   if (isTRUE(progress)) {
      pbar(sprintf("(%s %s bootstrap starred estimates)", compute, n.bs))
   }

   # ---------------- #
   # Step 2: Solve the quadratic program
   # ---------------- #

   # Evaluate the quadratic program
   qp.result <- tryCatch({
      qp.return <- beta.star.qp(data, lpmodel, beta.tgt, weight.mat,
                                beta.obs.bs, sigma.beta.obs, solver)
      list(beta.star = qp.return$beta.star,
           x.star = qp.return$x.star)
   }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
   }, error = function(e) {
      return(list(status = "error",
                  msg = e))
   })

   if (is.null(qp.result$status)) {
      beta.star <- qp.result$beta.star
      x.star <- qp.result$x.star
      msg <- NULL
   } else {
      beta.star <- NULL
      x.star <- NULL
      msg <- qp.result$msg$message
   }

   return(list(beta.star = beta.star,
               x.star = x.star,
               msg = msg))
}

#' Computes the restricted estimator in the \code{\link[lpinfer]{fsst}}
#' procedure
#'
#' @description This function computes the restricted estimator
#'   \eqn{\widehat{\bm{\beta}}^r_n} in the \code{\link[lpinfer]{fsst}}
#'   procedure.
#'
#' @importFrom Matrix t
#' @importFrom Matrix Matrix
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @inheritParams fsst.cone.lp
#'
#' @return Returns the optimal point and optimal value.
#'  \item{objval}{The optimal value.}
#'  \item{x}{The optimal point.}
#'
#' @export
#'
beta.r.compute <- function(n, lpmodel, beta.obs.hat, beta.tgt, beta.n,
                           beta.star, omega.i, indicator, solver) {
   # ---------------- #
   # Step 1: Construct the linear program
   # ---------------- #
   # Construct the zero matrices and vectors and a vector of ones
   p <- length(beta.n)
   if (!is.matrix(beta.obs.hat)) {
      q <- length(beta.obs.hat)
   } else {
      q <- nrow(beta.obs.hat)
   }
   d <- ncol(lpmodel$A.obs)
   ones.1p <- Matrix::Matrix(rep(1, p), nrow = 1, ncol = p, sparse = TRUE)
   ones.p1 <- Matrix::Matrix(rep(1, p), nrow = p, ncol = 1, sparse = TRUE)
   zero.1p <- Matrix::Matrix(rep(0, p), nrow = 1, ncol = p, sparse = TRUE)
   zero.p1 <- Matrix::Matrix(rep(0, p), nrow = p, ncol = 1, sparse = TRUE)
   zero.1d <- Matrix::Matrix(rep(0, d), nrow = 1, ncol = d, sparse = TRUE)
   zero.1q <- Matrix::Matrix(rep(0, q), nrow = 1, ncol = q, sparse = TRUE)
   zero.q1 <- Matrix::Matrix(rep(0, q), nrow = q, ncol = 1, sparse = TRUE)
   zero.pd <- Matrix::Matrix(rep(0, p * d), nrow = p, ncol = d, sparse = TRUE)
   zero.pp <- Matrix::Matrix(rep(0, p * p), nrow = p, ncol = p, sparse = TRUE)
   zero.pq <- Matrix::Matrix(rep(0, p * q), nrow = p, ncol = q, sparse = TRUE)
   zero.qq <- Matrix::Matrix(rep(0, q * q), nrow = q, ncol = q, sparse = TRUE)
   zero.pqq <- Matrix::Matrix(rep(0, (p - q) * q), nrow = q, ncol = p - q,
                              sparse = TRUE)
   diagm.pq <- Matrix::Matrix(diag(p - q), sparse = TRUE)
   iden.beta <- rbind(cbind(zero.qq, zero.pqq),
                      cbind(Matrix::t(zero.pqq), diagm.pq))

   # Construct the constraints matrix
   A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
   A.list <- list()
   A.list[[1]] <- asmat(cbind(sqrt(n) * diag(p), zero.pd, zero.pq, -omega.i,
                              zero.pp, A, zero.pd, zero.p1, zero.p1, zero.p1))
   A.list[[2]] <- asmat(cbind(sqrt(n) * diag(p), zero.pd, zero.pq, zero.pp,
                              -omega.i, zero.pd, -A, zero.p1, zero.p1, zero.p1))
   if (indicator == 0) {
      # Multiply A.mat1 and A.mat2 by t(A) if indicator == 0 (i.e. d < p)
      for (idx in 1:2) {
         A.list[[idx]] <- Matrix::t(A) %*% A.list[[idx]]
      }
   }
   A.list[[3]] <- asmat(cbind(diag(p), -A, zero.pq, zero.pp, zero.pp, zero.pd,
                              zero.pd, zero.p1, zero.p1, zero.p1))
   A.list[[4]] <- asmat(cbind(diag(p), zero.pd,
                              rbind(-diag(q), Matrix::t(zero.pqq)),
                              zero.pp, zero.pp, zero.pd, zero.pd, zero.p1,
                              zero.p1, zero.p1))
   A.list[[5]] <- asmat(cbind(zero.pp, zero.pd, zero.pq, diag(p), zero.pp,
                              zero.pd, zero.pd, ones.p1, zero.p1, zero.p1))
   A.list[[6]] <- asmat(cbind(zero.pp, zero.pd, zero.pq, - diag(p), zero.pp,
                              zero.pd, zero.pd, ones.p1, zero.p1, zero.p1))
   A.list[[7]] <- asmat(cbind(zero.pp, zero.pd, zero.pq, zero.pp, diag(p),
                              zero.pd, zero.pd, zero.p1, ones.p1, zero.p1))
   A.list[[8]] <- asmat(cbind(zero.pp, zero.pd, zero.pq, zero.pp, - diag(p),
                              zero.pd, zero.pd, zero.p1, ones.p1, zero.p1))
   A.list[[9]] <- asmat(cbind(zero.1p, zero.1d, zero.1q, zero.1p, zero.1p,
                              zero.1d, zero.1d, -1, 0, 1))
   A.list[[10]] <- asmat(cbind(zero.1p, zero.1d, zero.1q, zero.1p, zero.1p,
                               zero.1d, zero.1d, 0, -1, 1))
   A.mat <- Reduce(rbind, A.list)

   # Construct the rhs vector and the sense vector
   rhs.2 <- Matrix::Matrix(rep(0, p * 4 + 2), ncol = 1, sparse = TRUE)
   sense.2 <- rep(">=", 4 * p + 2)
   if (indicator == 0) {
      rhs.1 <- Reduce(rbind,
                      c(sqrt(n) * Matrix::t(A) %*% beta.star,
                        sqrt(n) * Matrix::t(A) %*% beta.star,
                        zero.p1,
                        zero.q1,
                        lpmodel$beta.shp,
                        beta.tgt))
      sense.1 <- rep("=", 2 * (d + p))
   } else {
      rhs.1 <- Reduce(rbind,
                      c(sqrt(n) * beta.n,
                        sqrt(n) * beta.n,
                        zero.p1,
                        zero.q1,
                        lpmodel$beta.shp,
                        beta.tgt))
      sense.1 <- rep("=", 4 * p)
   }

   # Combine the RHS and sense vectors
   rhs.mat <- rbind(rhs.1, rhs.2)
   sense.mat <- c(sense.1, sense.2)

   # Construct the objective function
   obj <- c(rep(0, 3 * (p + d) + q + 2), 1)

   # Lower bound
   lb <- c(rep(-Inf, p),
           rep(0, d),
           rep(-Inf, 2 * p + q),
           rep(0, 2 * d + 2),
           -Inf)

   # Set the arguments
   optim.arg <- list(Af = NULL,
                     bf = obj,
                     nf = 1,
                     A = A.mat,
                     rhs = rhs.mat,
                     sense = sense.mat,
                     modelsense = "min",
                     lb = lb,
                     beta.r.program = TRUE)

   # ---------------- #
   # Step 2: Solve the linear program
   # ---------------- #
   ans <- do.call(solver, optim.arg)
   objval <- ans$objval

   # ---------------- #
   # Step 3: Compute beta.r
   # ---------------- #
   x = ans$x[1:p]

   return(list(x = x,
               objval = objval))
}

#' Computes the range component of the test statistics
#'
#' @description This function computes the range component of the test
#'   statistics in the \code{\link[lpinfer]{fsst}} procedure.
#'
#' @importFrom Matrix norm
#'
#' @inheritParams fsst
#' @inheritParams beta.star.qp
#' @inheritParams fsst.beta.bs
#' @inheritParams sigma.summation
#' @param x.star The optimal value from \code{beta.star.qp}.
#' @param weight.mat.root The matrix square root of \code{weight.mat}.
#'
#' @return Returns the range test statistic.
#'  \item{range}{The sample range test statistic.}
#'
#' @export
#'
fsst.range <- function(n, beta.obs.hat, x.star, lpmodel, weight.mat.root) {
   # ---------------- #
   # Step 1: Compute the matrix inside the norm
   # ---------------- #
   A.obs.hat <- lpmodel$A.obs
   beta.obs.star <- A.obs.hat %*% x.star
   range.arg <- sqrt(n) * weight.mat.root %*% (beta.obs.hat - beta.obs.star)

   # ---------------- #
   # Step 2: Compute the range component
   # ---------------- #
   range <- Matrix::norm(range.arg, type = "I")

   return(range)
}

#' Bootstrap procedure of computing the range component
#'
#' @description This function computes the bootstrap estimates of the range
#'   components.
#'
#' @import furrr progressr
#'
#' @inheritParams fsst
#' @inheritParams fsst.beta.bs
#' @inheritParams fsst.range
#' @inheritParams fsst.weight.matrix
#' @inheritParams full.beta.bs
#' @param x.star.bs The bootstrap estimates of \code{x.star}.
#'
#' @return Return the following list of objects:
#'   \item{range.n.list}{The list of bootstrap range components.}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{new.error}{The number of new errors.}
#'   \item{error.id}{The indices of the bootstrap replications that give
#'     errors.}
#'
#' @export
#'
fsst.range.bs <- function(n, lpmodel, beta.obs.hat, beta.obs.bs, x.star,
                          x.star.bs, weight.mat.root, R, progress,
                          df.error, eval.count) {

   # ---------------- #
   # Step 1: Compute the bootstrap range estimates
   # ---------------- #
   # Combine the bootstrap replications of 'beta.obs' and 'x.star'
   beta.x.star <- mapply(list, beta.obs.bs, x.star.bs, SIMPLIFY = FALSE)

   # Compute the bootstrap range estimates
   n.bs <- length(beta.x.star)
   progressr::with_progress({
      if (isTRUE(progress)) {
         pbar <- progressr::progressor(along = 1:n.bs)
      } else {
         pbar <- NULL
      }

      range.return <- furrr::future_map(beta.x.star,
                                        .f = fsst.range.bs.fn,
                                        n = n,
                                        lpmodel = lpmodel,
                                        beta.obs.hat = beta.obs.hat,
                                        x.star = x.star,
                                        weight.mat.root = weight.mat.root,
                                        pbar = pbar,
                                        progress = progress,
                                        n.bs = n.bs,
                                        eval.count = eval.count)
   })

   # Extract the test statistics
   range.n.list <- unlist(sapply(range.return, "[", "Ts"))
   error.list <- sapply(range.return, "[", "msg")
   R.succ <- length(range.n.list)
   R.eval <- length(beta.obs.bs)
   new.error <- R.eval - R.succ

   # ---------------- #
   # Step 2: Consolidate the error messages
   # ---------------- #
   if (R.eval != R.succ) {
      # Create data.frame for error messages
      df.error1 <- data.frame(id = NA,
                              lambda = NA,
                              message = unlist(error.list))

      # Match the id of the error messages
      df.error1 <- error.id.match(error.list, df.error1)
      if (!is.null(df.error)) {
         # Case 1: There are errors in this function and the earlier parts
         df.error <- rbind(df.error, df.error1)
         rownames(df.error) <- 1:nrow(df.error)
      } else {
         # Case 2: There are no errors in the earlier parts
         df.error <- df.error1
      }
      new.error <- R.eval - R.succ
      error.id <- df.error1$id
   } else {
      # Case 3: There are no errors in this part but there are errors in the
      # earlier parts
      df.error <- df.error
      new.error <- 0
      error.id <- NULL
   }

   return(list(range.n.list = range.n.list,
               new.error = new.error,
               df.error = df.error,
               error.id = error.id))
}

#' Computes one bootstrap estimate for the range component of the test
#' statistic
#'
#' @description This function computes the one bootstrap estimate for the
#'   range component of the test statistic.
#'
#' @inheritParams fsst
#' @inheritParams dkqs.bs
#' @inheritParams fsst.beta.bs
#' @inheritParams fsst.range
#' @inheritParams fsst.range.bs
#' @inheritParams fsst.weight.matrix
#' @inheritParams dkqs.bs.fn
#' @param beta.x.star A list that combines the bootstrap estimates of
#'   \code{beta.star} and \code{x.star}.
#'
#' @return Return the following list of objects:
#'   \item{Ts}{The bootstrap range test statistic.}
#'   \item{msg}{An error message (if applicable).}
#'
#' @export
#'
fsst.range.bs.fn <- function(beta.x.star, n, lpmodel, beta.obs.hat, x.star,
                             weight.mat.root, pbar, progress, n.bs,
                             eval.count) {
   # ---------------- #
   # Step 1: Print progress bar
   # ---------------- #
   if (eval.count == 0) {
      compute <- "Computing"
   } else {
      compute <- "Recomputing"
   }

   if (isTRUE(progress)) {
      pbar(sprintf(paste0("(%s %s bootstrap estimates of the range test ",
                          "statistic)"), compute, n.bs))
   }

   # Compute the differences in 'beta.bs' and 'x.star'
   beta.bs.1 <- beta.x.star[[1]] - beta.obs.hat
   x.star.bs.1 <- beta.x.star[[2]] - x.star

   # Compute one range estimate
   result.range <- tryCatch({
      range.n.return <- fsst.range(n, beta.bs.1, x.star.bs.1, lpmodel,
                                   weight.mat.root)
      list(range = range.n.return)
   }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
   }, error = function(e) {
      return(list(status = "error",
                  msg = e))
   })

   if (is.null(result.range$status)) {
      Ts <- result.range$range
      msg <- NULL
   } else {
      Ts <- NULL
      msg <- result.range$msg$message
   }

   return(list(Ts = Ts,
               msg = msg))
}

#' Computes the bootstrap estimates of the cone component of the test
#'   statistics
#'
#' @description This function computes the bootstrap estimates of the cone
#'   components.
#'
#' @import furrr progressr
#'
#' @inheritParams fsst
#' @inheritParams fsst.cone.lp
#' @inheritParams fsst.beta.star.bs
#' @param beta.star The starred version of the \code{beta.n} vector, i.e.
#'   \eqn{\widehat{\bm{\beta}}^\star_n}.
#' @param beta.star.list This corresponds to the bootstrap estimates of
#'   \code{beta.star}.
#' @param R.succ The number of successful bootstrap replications.
#' @param beta.r The restricted estimator.
#' @param data.driven A boolean variable. This indicates whether the
#'   data-driven problem is being considered.
#'
#' @return Return the following list of objects:
#'   \item{cone.n.list}{The list of bootstrap cone components.}
#'   \item{new.error}{The number of new errors.}
#'   \item{df.error}{Table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{error.id}{An indices of the bootstrap replications that give errors.}
#'   \item{R.succ}{A number of successful bootstrap replications.}
#'
#' @export
#'
fsst.cone.bs <- function(n, omega.i, beta.n, beta.star, lpmodel, R.succ,
                         lambda, indicator, beta.r, beta.star.list, solver,
                         progress, df.error, eval.count, data.driven = FALSE) {
   # ---------------- #
   # Step 1: Compute the bootstrap cone estimates
   # ---------------- #
   n.bs <- length(beta.star.list[-1])
   progressr::with_progress({
      if (isTRUE(progress)) {
         pbar <- progressr::progressor(along = 1:n.bs)
      } else {
         pbar <- NULL
      }
      cone.return <- furrr::future_map(beta.star.list[-1],
                                        .f = fsst.cone.bs.fn,
                                        n = n,
                                        omega.i = omega.i,
                                        beta.n = beta.n,
                                        beta.star = beta.star,
                                        lpmodel = lpmodel,
                                        lambda = lambda,
                                        indicator = indicator,
                                        beta.r = beta.r,
                                        solver = solver,
                                        pbar = pbar,
                                        progress = progress,
                                        n.bs = n.bs,
                                        eval.count = eval.count,
                                        data.driven = data.driven)
   })

   cone.n.list <- unlist(sapply(cone.return, "[", "Ts"), use.names = FALSE)
   error.list <- sapply(cone.return, "[", "msg")
   R.succ <- length(cone.n.list)
   R.eval <- length(beta.star.list) - 1
   new.error <- R.eval - R.succ

   # ---------------- #
   # Step 2: Consolidate the error messages
   # ---------------- #
   if (R.eval != R.succ) {
      # Create data.frame for error messages
      df.error1 <- data.frame(id = NA,
                              lambda = NA,
                              message = unlist(error.list))

      # Match the id of the error messages
      df.error1 <- error.id.match(error.list, df.error1)
      if (!is.null(df.error)) {
         # Case 1: There are errors in this function and the earlier parts
         df.error <- rbind(df.error, df.error1)
         rownames(df.error) <- 1:nrow(df.error)
      } else {
         # Case 2: There are no errors in the earlier parts
         df.error <- df.error1
      }
      new.error <- R.eval - R.succ
      error.id <- df.error1$id
   } else {
      # Case 3: There are no errors in this part but there are errors in the
      # earlier parts
      df.error <- df.error
      new.error <- 0
      error.id <- NULL
   }

   return(list(cone.n.list = cone.n.list,
               new.error = new.error,
               df.error = df.error,
               error.id = error.id,
               R.succ = R.succ))
}

#' Computes one bootstrap estimate of the cone component of the test
#'   statistics
#'
#' @description This function computes one bootstrap estimate of the cone
#'   component of the test statistics.
#'
#' @inheritParams fsst
#' @inheritParams dkqs.bs
#' @inheritParams fsst.beta.bs
#' @inheritParams fsst.cone.bs
#' @inheritParams fsst.weight.matrix
#' @inheritParams dkqs.bs.fn
#' @param beta.star.bs One bootstrap estimate of the \code{beta.star} object.
#'
#' @return Return the following list of objects:
#'   \item{Ts}{The bootstrap cone test statistic.}
#'   \item{msg}{An error message (if applicable).}
#'
#' @export
#'
fsst.cone.bs.fn <- function(beta.star.bs, n, omega.i, beta.n, beta.star,
                            lpmodel, lambda, indicator, beta.r, solver,
                            pbar, progress, n.bs, eval.count, data.driven) {
   # ---------------- #
   # Step 1: Print progress bar
   # ---------------- #
   if (eval.count == 0) {
      compute <- "Computing"
   } else {
      compute <- "Recomputing"
   }

   if (isTRUE(data.driven)) {
      what.compute <- "data-driven lambda"
   } else {
      what.compute <- paste0(n.bs,
                             " bootstrap estimates of the cone test statistic")
   }

   if (isTRUE(progress)) {
      pbar(sprintf("(%s %s)", compute, what.compute))
   }

   # Compute one bootstrap cone estimate
   result.cone <- tryCatch({
      beta.new <- beta.star.bs - beta.star + lambda * beta.r
      cone.return <- fsst.cone.lp(n, omega.i, beta.n, beta.new, lpmodel,
                                  indicator, solver)
      list(objval = cone.return$objval,
           msg = cone.return$status)
   }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
   }, error = function(e) {
      return(list(status = "error",
                  msg = e))
   })

   if (is.null(result.cone$status)) {
      Ts <- result.cone$objval
      msg <- NULL
   } else {
      Ts <- NULL
      msg <- result.cone$msg$message
   }

   if (length(Ts) == 0) {
      msg <- result.cone$msg
   }

   return(list(Ts = Ts,
               msg = msg))
}

#' Calculates the \eqn{p}-value for the \code{\link[lpinfer]{fsst}} procedure
#'
#' @description This function computes the \eqn{p}-value of the test based on
#'    the bootstrap estimates for the \code{\link[lpinfer]{fsst}} procedure.
#'
#' @param range.n The range component of the sample test statistic.
#' @param cone.n The cone component of the sample test statistic.
#' @param range.n.list The bootstrap estimates of the range component in the
#'   test statistics.
#' @param cone.n.list The bootstrap estimates of the cone component in the
#'   test statistics.
#' @param alpha The significance level.
#' @inheritParams fsst
#'
#' @return Returns the \eqn{p}-value and the decision.
#'   \item{pval}{The \eqn{p}-value.}
#'   \item{decision}{The decision to reject or not.}
#'   \item{T.bs}{List of bootstrap test statistics.}
#'
#' @export
#'
fsst.pval <- function(range.n, cone.n, range.n.list, cone.n.list, R,
                      alpha = .05) {
   # ---------------- #
   # Step 1: Compute the test statistic and the bootstrap test statistics
   # ---------------- #
   T.n <- max(range.n, cone.n)
   T.bs <- NULL
   for (i in 1:R) {
      T.bs <- c(T.bs, max(range.n.list[i], cone.n.list[i]))
   }

   # ---------------- #
   # Step 2: Decision
   # ---------------- #
   # Compute p-value
   pval <- mean(T.bs >= T.n)

   # Decision
   if (pval > alpha) {
      decision <- 1
   } else {
      decision <- 0
   }
   return(list(pval = pval,
               decision = decision,
               T.bs = T.bs))
}

#' Checks and updates the input in \code{fsst}
#'
#' @importFrom methods is
#'
#' @description This function checks and updates the input from the user in the
#'    \code{\link[lpinfer]{fsst}} function. If there is any invalid input,
#'    the function will be terminated and error messages will be printed.
#'
#' @inheritParams fsst
#' @inheritParams dkqs
#'
#' @return Returns the updated parameters back to the function \code{fsst}.
#'   The following information are updated:
#'    \itemize{
#'       \item{\code{data}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{test.logical}}
#'       \item{\code{logical.lb}}
#'       \item{\code{logical.ub}}
#'    }
#'
#' @export
#'
fsst.check <- function(data, lpmodel, beta.tgt, R, Rmulti, lambda, rho, n,
                       weight.matrix, solver, progress,
                       sqrtm.method, sqrtm.tol, previous.output) {

   # ---------------- #
   # Step 1: Check data
   # ---------------- #
   # Check data. If data is NULL, check if n is a positive integer
   if (!is.null(data)) {
      data <- check.dataframe(data)
   } else {
      check.samplesize(n, "n")
   }

   # ---------------- #
   # Step 2: Check lpmodel
   # ---------------- #
   lpmodel <- check.lpmodel(data = data,
                            lpmodel = lpmodel,
                            name.var = "lpmodel",
                            A.tgt.cat = "matrix",
                            A.obs.cat = "matrix",
                            A.shp.cat = "matrix",
                            beta.obs.cat = c("function_mat",
                                             "list",
                                             "function_obs_var"),
                            beta.shp.cat = "matrix",
                            R = R)
   # Length of the \beta(P) and beta.obs vector
   nbobs <- length(lpmodel.beta.eval(data, lpmodel$beta.obs, 1)$beta.obs)
   p <- nbobs + length(lpmodel$beta.shp) + 1

   # ---------------- #
   # Step 3: Check solver
   # ---------------- #
   solver.return <- check.solver(solver, "solver")
   solver <- solver.return$solver
   solver.name <- solver.return$solver.name

   # ---------------- #
   # Step 4: Check numerics
   # ---------------- #
   check.numeric(beta.tgt, "beta.tgt")
   check.positiveinteger(R, "R")
   check.positive(rho, "rho")

   # Check Rmulti
   check.numrange(Rmulti, "Rmulti", "closed", 1, "open", Inf)

   # ---------------- #
   # Step 5: Check lambda
   # ---------------- #
   # If NA is present in lambda, drop it before checking it because NA refers
   # to using the data-driven lambda in FSST.
   if (NA %in% lambda) {
      lambda.temp <- lambda[!is.na(lambda)]
   } else {
      lambda.temp <- lambda
   }

   # Check lambda without NA. If the user did not provide a lambda, then the
   # length of lambda.temp is 0, and there is nothing to check.
   if (length(lambda.temp) == 1) {
      check.numrange(lambda.temp, "lambda", "closed", 0, "closed", 1)
   } else if (length(lambda.temp) > 1) {
      for (i in 1:length(lambda.temp)) {
         check.numrange(lambda.temp[i], "lambda", "closed", 0, "closed", 1)
      }
   }

   # ---------------- #
   # Step 6: Check whether beta.tgt is within the logical bounds
   # ---------------- #
   test.return <- check.betatgt(data, lpmodel, beta.tgt, solver)
   test.logical <- test.return$inout
   logical.lb <- test.return$lb
   logical.ub <- test.return$ub

   # ---------------- #
   # Step 7: Check the arguments related to the `sqrtm` function
   # ---------------- #
   # Check if `sqrtm.method` is correct
   sqrtm.errmsg <- paste0("'sqrtm.method' has to be a function that takes ",
                          "one argument that accepts a square matrix of size ",
                          "k x k and returns a square matrix of size k x k, ",
                          "where k is the length of the beta(P) vector or ",
                          "the length of the 'beta.obs' component.")
   if (!inherits(sqrtm.method, "function")) {
      # Check if `sqrtm.method` is a function
      stop(sqrtm.errmsg)
   } else {
      # If it is a function, check whether it is correct
      if (length(formals(sqrtm.method)) != 1) {
         # Return an error if it has more than one argument
         stop(sqrtm.errmsg)
      } else {
         # Check if it accepts a matrix of size p x p and nbobs x nbobs
         for (idx in c(p, nbobs)) {
            tmp.result <- tryCatch({
               tmp.return <- sqrtm.method(diag(idx))
               list(mat = tmp.return)
            }, warning = function(w) {
               return(list(status = "warning",
                           msg = w))
            }, error = function(e) {
               return(list(status = "error",
                           msg = e))
            })

            if (!is.null(tmp.result$status)) {
               # If tmp.result$status is not NULL, then there's some problem
               # with the function
               stop(sqrtm.errmsg)
            } else {
               # Check if it returns a matrix of size p x p
               tmp.nr <- nrow(tmp.result$mat)
               tmp.nc <- ncol(tmp.result$mat)
               if (!(is.numeric(tmp.nr) & is.numeric(tmp.nc))) {
                  # If the number of rows or columns is not numeric, return error
                  stop(sqrtm.errmsg)
               }
               if (tmp.nr != tmp.nc) {
                  # Return an error if the matrix obtained is not a square
                  stop(sqrtm.errmsg)
               } else {
                  # Return an error if the matrix is not of size p x p
                  # or nbobs x nbobs
                  if (tmp.nr != idx) {
                     stop(sqrtm.errmsg)
                  }
               }
            }
         }
      }
   }

   # Check if the tolerance level is positive
   check.positive(sqrtm.tol, "sqrtm.tol")

   # ---------------- #
   # Step 8: Check previous.output
   # ---------------- #
   prev.out.msg1 <- paste0("'previous.output' has to be either 'NA' or ",
                           "a list of objects returned from a previous ",
                           "evaluation of the FSST test that contains the ",
                           "'omega.i' matrix.")
   prev.out.msg2 <- "Therefore, the 'omega.i' matrix will be computed. "
   if (!is.list(previous.output)) {
      if (is.null(previous.output)) {
         warning(paste0(prev.out.msg1, " ", prev.out.msg2), immediate. = TRUE)
      } else if (!is.na(previous.output)) {
         warning(paste0(prev.out.msg1, " ", prev.out.msg2), immediate. = TRUE)
      }
      omega.i <- NA
   } else {
      if (is.null(previous.output$omega.i)) {
         # It has to contain 'omega.i' if previous.output is a list.
         omega.i <- NA
         warning(paste0(prev.out.msg1, " ", prev.out.msg2), immediate. = TRUE)
      } else {
         # Check if the omega.i object is correct
         # It can be a square 'data.frame', 'matrix' or a 'sparseMatrix'.
         omega.i <- previous.output$omega.i
         if (!(is.matrix(omega.i) | is.data.frame(omega.i) |
             methods::is(omega.i, "sparseMatrix"))) {
            omega.i <- NA
            warning(paste0("The class of the 'omega.i' matrix in the list ",
                           "'previous.output' has to be one of the ",
                           "followings: data.frame, matrix, or sparseMatrix. ",
                           prev.out.msg2),
                    immediate. = TRUE)
         } else {
            # Check if the dimension of the omega.i object is correct,
            # i.e., it has to be a square matrix; the number of rows and
            # columns have to equal the length of \beta(P)
            omega.i.dim <- dim(omega.i)
            if (!((omega.i.dim[1] == omega.i.dim[2]) & omega.i.dim[1] == p)) {
               warning(paste0("The dimension of the 'omega.i' matrix is ",
                              "incorrect. ", prev.out.msg2),
                       immediate. = TRUE)
               omega.i <- NA
            }
         }
      }
   }

   # ---------------- #
   # Step 9: Return updated items
   # ---------------- #
   return(list(solver = solver,
               solver.name = solver.name,
               data = data,
               test.logical = test.logical,
               logical.lb = logical.lb,
               logical.ub = logical.ub,
               omega.i = omega.i))
}

#' Print results from \code{\link[lpinfer]{fsst}}
#'
#' @description This function prints the \eqn{p}-values from
#'   \code{\link[lpinfer]{fsst}}.
#'
#' @param x The output objects returned from \code{\link[lpinfer]{fsst}}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.fsst <- function(x, ...) {

   if (x$test.logical == 1) {
      # Case 1: 'beta.tgt' is within the logical bound
      # Print the p-values
      df.pval <- x$pval
      if (nrow(df.pval) == 1) {
         pv <- "p-value"
         # Indicates if it is obtained by a data-driven lambda
         if (!is.null(x$lambda.data)) {
            pv <- paste(pv, "(by data-driven 'lambda')")
         }
         cat(sprintf("%s: %s\n", pv, round(df.pval[1, 2], digits = 5)))
      } else {
         cat("p-values:\n")
         # Label the data-driven lambda with a "*" if it is used
         dfl <- fsst.label.lambda(round(df.pval$`lambda`, digits = 5),
                                  x$lambda.data)
         df.pval$`lambda` <- dfl$lambdas
         print(df.pval, row.names = FALSE)

         # Print the message for data-driven lambda if necessary
         if (!is.null(dfl$msg)) {
            cat(paste0("   ", dfl$msg))
         }
      }
   } else {
      # Case 2: 'beta.tgt' is outside the logical bound
      infeasible.pval.msg()
   }
}

#' Summary of results from \code{\link[lpinfer]{fsst}}
#'
#' @description This function prints a summary of the results obtained from
#'   \code{\link[lpinfer]{fsst}}.
#'
#' @inheritParams print.fsst
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.fsst <- function(x, ...) {

   if (x$test.logical == 1) {
      # Case 1: 'beta.tgt' is within the logical bound
      # Print the p-values
      df.pval <- x$pval
      n.pval <- nrow(df.pval)
      if (n.pval == 1) {
         print.fsst(x)
      } else {
         cat("p-values:\n")
         df.pval.2 <- data.frame(matrix(vector(), nrow = 1, ncol = n.pval + 1))
         # Label the data-driven lambda with a "*" if it is used
         dfl <- df.pval$lambda
         dfl <- fsst.label.lambda(dfl, x$lambda.data)
         colnames(df.pval.2) <- c("    lambda    ", dfl$lambdas)
         df.pval.2[1, ] <- c("    p-value   ", df.pval[, 2])
         print(df.pval.2, row.names = FALSE)
      }

      # Print the sample and bootstrap test statistics
      if (!is.null(nrow(x$cv.table))) {
         cat("Sample and quantiles of bootstrap test statistics: \n")
         cv.tab <- x$cv.table
         cv.tab[is.na(cv.tab)] <- ""
         cv.tab[, 1] <- paste0("   ", cv.tab[, 1], " ")
         cv.tab[, 2] <- paste0(cv.tab[, 2], "  ")
         colnames(cv.tab)[2] <- paste0(colnames(cv.tab)[2], "  ")
         # Label the data-driven lambda with a "*" if it is used
         cvlambda <- as.numeric(colnames(cv.tab)[-c(1, 2)])
         if (is.null(x$lambda.data)) {
            x.lambda.data <- x$lambda.data
         } else {
            x.lambda.data <- round(x$lambda.data, digits = 5)
         }
         cvlambda <- fsst.label.lambda(cvlambda, x.lambda.data)
         colnames(cv.tab)[-c(1, 2)] <- cvlambda$lambdas
         print(cv.tab, row.names = FALSE)
      }

      # Regularization parameters
      cat("Regularization parameters: \n")
      cat(sprintf("   - Input value of rho: %s\n",
                  round(x$rho, digits = 5)))
      cat(sprintf(paste0("   - Regularization parameter for the Cone ",
                         "studentization matrix: %s\n"),
                  round(x$rhobar.i, digits = 5)))

      # Print solver
      cat(sprintf("Solver: %s\n", x$solver))

      # Number of successful bootstrap replications
      cat(sprintf("Number of successful bootstrap replications: %s\n",
                  x$R.succ))

      # Number of failed bootstrap replications
      if (!is.null(x$df.error)) {
         if (nrow(x$df.error) != 0) {
            nerr <- nrow(x$df.error)
            errstring <- "Number of failed bootstrap"
            if (nerr == 1) {
               cat(sprintf(paste(errstring, "replication: %s\n"), nerr))
            } else {
               cat(sprintf(paste(errstring, "replications: %s\n"), nerr))
            }
         }
      }

      # How the variance matrix is estimated
      cat(sprintf(paste0("The asymptotic variance of the observed component ",
                         "of the beta vector is approximated from the %s.\n"),
                  x$var.method))

      # Print the message for data-driven lambda if necessary
      if (!is.null(nrow(x$cv.table))) {
         if (!is.null(cvlambda$msg)) {
            cat(cvlambda$msg)
         }
      }
   } else if (x$test.logical == 0) {
      # Case 2: 'beta.tgt' is outside the logical bound
      infeasible.pval.msg()
      cat(sprintf("\nSolver used: %s\n", x$solver.name))
   }
}

#' Data-driven choice of \code{lambda} in the \code{\link[lpinfer]{fsst}}
#' procedure
#'
#' @description This function provides a data-driven choice of \code{lambda}
#'   in the \code{\link[lpinfer]{fsst}} procedure.
#'
#'
#' @inheritParams fsst.cone.bs
#' @inheritParams fsst.beta.star.bs
#' @param p The length of the vector \eqn{\bm{\beta}}.
#' @param d The number of columns of the \eqn{\bm{A}} matrix.
#'
#' @details \eqn{\alpha_n} is set as \eqn{1} if the number of observations in
#'   the data set is less than 16.
#'
#' @return Returns the data-driven \code{lambda} and the error messages.
#'   \item{lambda}{The data-driven \code{lambda}}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{new.error}{The number of new errors.}
#'   \item{error.id}{The indices of the bootstrap replications that give
#'     errors.}
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'
#' @export
#'
fsst.lambda <- function(n, omega.i, beta.n, beta.star, lpmodel, R.succ,
                        beta.star.list, solver, progress, df.error, p, d,
                        eval.count) {
   # ---------------- #
   # Step 1: Compute the bootstrap cone estimates with lambda = 0, beta.r = 0
   # ---------------- #
   beta.r <- rep(0, length(beta.n))
   if (d >= p) {
      indicator <- 1
   } else {
      indicator <- 0
   }
   fsst.cone.return <- fsst.cone.bs(n, omega.i, beta.n, beta.star, lpmodel,
                                    R.succ, 0, indicator, beta.r,
                                    beta.star.list, solver, progress, df.error,
                                    eval.count, data.driven = TRUE)

   # ---------------- #
   # Step 2: Consolidate the error messages
   # ---------------- #
   new.error <- fsst.cone.return$new.error
   df.error <- fsst.cone.return$df.error
   error.id <- fsst.cone.return$error.id
   R.succ <- fsst.cone.return$R.succ

   # ---------------- #
   # Step 3: Consolidate the results
   # ---------------- #
   if (new.error != 0) {
      # If there is an error, return the new error messages
      lambda <- NULL
   } else {
      # Otherwise, compute the data-driven lambda
      if (n < 16) {
         alpha <- 1
      } else {
         alpha <- 1/sqrt(log(log(n)))
      }
      Tn <- sort(fsst.cone.return$cone.n.list)
      q <- Tn[ceiling(length(Tn) * (1 - alpha))]
      lambda <- min(1, 1/q)
   }

   return(list(lambda = lambda,
               new.error = new.error,
               df.error = df.error,
               error.id = error.id,
               R.succ = R.succ))
}

#' Indicates the data-driven \code{lambda} in the output
#'
#' @description This function labels the data-driven \code{lambda} in the
#'   \code{print} or \code{summary} output and returns an indicative message.
#'   This only affects the output messages, not the objects returned.
#'
#' @param lambdas The vector of \code{lambda}.
#' @param lambda.data The data-driven lambda.
#'
#' @return Returns the updated vector of \code{lambda} where the data-driven
#'   \code{lambda} is labeled with a star and the corresponding message
#'   \item{lambdas}{The updated \code{lambda}.}
#'   \item{msg}{An indicative message.}
#'
#' @export
#'
fsst.label.lambda <- function(lambdas, lambda.data) {
   if (!is.null(lambda.data)) {
      lambdas[lambdas %in% lambda.data] <-
         paste(round(lambda.data, digits = 5), "*")
      msg <- "\n* refers to the data-driven 'lambda' parameter.\n"
   } else {
      msg <- NULL
   }
   return(list(lambdas = lambdas,
               msg = msg))
}

#' Checks whether the matrix square root is correct
#'
#' @description This function is used to check whether the matrix square
#'   root is correct. This is done by checking whether the Frobenius
#'   norm is smaller than the tolerance level, i.e., when \eqn{A} is the
#'   give matrix, \eqn{B} is the matrix square root obtained from the
#'   given \code{sqrtm.method} function, and \eqn{\epsilon} is the
#'   tolerance level, the FSST test checks whether
#'   \eqn{||A - BB||_F < \epsilon}. If this does not hold, the FSST test will
#'   use the \code{\link[expm]{sqrtm}} function from the \code{expm} package
#'   to obtain the matrix square root.
#'
#' @import expm
#' @importFrom Matrix norm
#'
#' @inheritParams fsst
#' @param mat The matrix that we wish to obtain the matrix square root.
#' @param mat.name The name of the matrix for \code{mat}.
#'
#' @return Returns the matrix square root.
#'   \item{mat.sqrt}{The matrix square root.}
#'
checkupdate.matrixroot <- function(mat, mat.name, sqrtm.method, sqrtm.tol) {
   # Obtain the square root by the `sqrtm.method`
   sqrtm.tmp <- sqrtm.method(as.matrix(mat))

   # Check whether the matrix square root is correct
   if (Matrix::norm(mat - sqrtm.tmp %*% sqrtm.tmp, type = "f") >= sqrtm.tol) {
      warning(paste0("In computing the square root of the ", mat.name,  ", ",
                     "the difference bewteen the matrix square root obtained ",
                     "from the function passed to 'sqrtm.method' and the ",
                     "true matrix square root is larger than the tolerance ",
                     "level. The matrix square root will now be computed ",
                     "by the expm::sqrtm function."),
              immediate. = TRUE)
      mat.sqrt <- expm::sqrtm(mat)
   } else {
      mat.sqrt <- sqrtm.tmp
   }
   return(mat.sqrt)
}
