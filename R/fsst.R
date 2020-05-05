#' Conducts inference using the `\code{fsst}` procedure
#'
#' @description This module conducts inference in linear programs using the 
#'   `\code{fsst}` procedure by Fang, Santos, Shaikh and Torgovitsky (2020).
#' 
#' @inheritParams dkqs
#' @param weight.matrix The option used in the weighting matrix.
#'   
#' @return Return a list of statistics for the `\code{fsst}` procedure:
#'   \item{pval}{\eqn{p}-value}
#'   \item{cores}{Number of cores used}
#'   \item{call}{Information used to call the function}
#'   \item{range}{The range test statistic}
#'   \item{cone}{The cone test statistic}
#'   \item{range.n.list}{The list of bootstrap range test statistics}
#'   \item{cone.n.list}{The list of bootstrap cone test statistics}
#'   \item{solver}{Name of the solver used}
#'   \item{test}{Test statistic used}
#'   
#' @export
#'
fsst <- function(data = NULL, lpmodel, beta.tgt, R, alpha = .05, lambda, rho, 
                 n = NULL, weight.matrix = "diag", solver, cores, progress){
   
   # ---------------- #
   # Step 1: Update call, check and update the arguments
   # ---------------- #
   # Obtain call information
   call <- match.call()
   
   # Check the arguments
   fsst.return <- fsst.check(data, lpmodel, beta.tgt, R, lambda, rho, n,
                             weight.matrix, solver, cores, progress)
   
   # Update the arguments
   solver <- fsst.return$solver
   solver.name <- fsst.return$solver.name
   cores <- fsst.return$cores
   
   # The user must either provide the data or n
   if (is.null(data)){
      n <- n
   } else {
      n <- nrow(data)
   }
   
   # ---------------- #
   # Step 2: Estimate beta.star
   # ---------------- #
   ### 2(a) Compute beta(P)
   beta.obs.return <- lpmodel.beta.eval(data, lpmodel$beta.obs, 1)
   beta.obs.hat <- beta.obs.return[[1]]
   beta.shp.hat <- lpmodel.eval(data, lpmodel$beta.shp, 1)
   beta.n <- c(unlist(beta.obs.hat), beta.shp.hat, beta.tgt)
   
   ### 2(b) Estimate sigma.beta.obs
   # If the user provided bootstrap estimates of beta, use it to compute sigma
   if (class(lpmodel$beta.obs) == "list"){
      if (class(lpmodel$beta.obs[[1]]) == "list"){
         beta.obs.hat <- lpmodel$beta.obs[[1]][[1]]
         sigma.beta.obs <- as.matrix(lpmodel$beta.obs[[1]][[2]])
      } else {
         beta.obs.hat <- lpmodel$beta.obs
         sigma.beta.obs <- sigma.summation(n, lpmodel$beta.obs) 
      }
      beta.obs.bs <- lpmodel$beta.obs[2:(R+1)]
      beta.n.bs <- list()
      for (i in 1:R){
         beta.n.bs[[i]] <- c(beta.obs.bs[[i]], beta.shp.hat, beta.tgt)
      }
   } else {
      if (cores == 1){
         sigma.return <- sigma.est(n, data, beta.obs.hat, lpmodel, R)
         sigma.beta.obs <- sigma.return$sigma.hat
         beta.obs.bs <- sigma.return$beta.obs.bs
         beta.n.bs <- full.beta.bs(lpmodel, beta.tgt, beta.obs.bs, R)
      } else {
         sigma.return <- sigma.est.parallel(n, data, beta.obs.hat, lpmodel, R, 
                                            cores, progress)
         sigma.beta.obs <- sigma.return$sigma.hat
         beta.obs.bs <- sigma.return$beta.obs.bs
         beta.n.bs <- full.beta.bs(lpmodel, beta.tgt, beta.obs.bs, R)
      }
   }
   
   ### 2(c) Compute the beta.sigma
   n.beta1 <- nrow(sigma.beta.obs)
   n.beta23 <- length(c(beta.shp.hat, beta.tgt))
   zero.12 <- matrix(rep(0, n.beta1*n.beta23), nrow = n.beta1)
   zero.21 <- t(zero.12)
   zero.22 <- matrix(rep(0, n.beta23^2), nrow = n.beta23)
   beta.sigma <- rbind(cbind(sigma.beta.obs, zero.12),
                       cbind(zero.21, zero.22))
   
   ### 2(d) Compute beta.star
   p <- length(beta.n)
   d <- ncol(lpmodel$A.tgt)
   if (d >= p){
      beta.star <- beta.n
      
      # Construct the bootstrap estimates, which are the same for all of them
      beta.star.bs <- beta.n.bs
   } else {
      # Solve the quadratic program
      beta.star <- beta.star.qp(data, lpmodel, beta.tgt, weight.matrix,
                                   beta.obs.hat, sigma.beta.obs, solver)
      
      # Construct bootstrap estimates of beta.star
      beta.star.bs <- list()
      for (i in 1:R){
         beta.star.bs[[i]] <- beta.star.qp(data, lpmodel, beta.tgt,
                                              weight.matrix, beta.obs.bs[[i]],
                                              sigma.beta.obs, solver)
      }
   }
   
   # Consolidate the bootstrap estimators into a list
   beta.star.l <- list(beta.star = beta.star)
   beta.star.list <- c(beta.star.l, beta.star.bs)
   
   # ---------------- #
   # Step 3: Studentization
   # ---------------- #
   # Obtain bootstrap star
   if (d >= p){
      sigma.star <- beta.sigma
      sigma.star.diff <- matrix(rep(0, p^2), nrow = p)
   } else {
      if (cores == 1){
         sigma.star <- sigma.summation(n, beta.star.list)
      } else {
         sigma.star <- sigma.summation.parallel(n, beta.star.list, cores,
                                                progress, 1)
      }
      beta.diff.bs <- list()
      beta.diff.bs[[1]] <- beta.n - beta.star
      for (i in 1:R){
         beta.diff.bs[[i+1]] <- beta.n.bs[[i]] - beta.star.bs[[i]]
      }
      
      if (cores == 1){
         sigma.star.diff <- sigma.summation(n, beta.diff.bs)
      } else {
         sigma.star.diff <- sigma.summation.parallel(n, beta.diff.bs, cores,
                                                     progress, 2)
      }
   }
   
   # Compute the matrix square root
   #rho.bar <- max(abs(sigma.star)) * rho
   rho.bar <- base::norm(sigma.star, type = "f") * rho
   rho.sup <- base::norm(sigma.star.diff, type = "f") * rho
   
   omega.i <- expm::sqrtm(sigma.star + rho.bar * diag(nrow(sigma.star)))
   if (d >= p){
      omega.e <- sigma.star.diff
   } else {
      omega.e <- expm::sqrtm(sigma.star.diff + rho.sup * 
                                diag(nrow(sigma.star)))
   }

   # ---------------- #
   # Step 4: Test statistic
   # ---------------- #
   # Compute range.n
   if (d >= p){
      range.n <- list(objval = 0,
                      x = NA)
      cone.n <- cone.n.lp(n, omega.i, beta.n, beta.star, lpmodel, 1, solver)
   } else {
      range.n <- rangen.lp(n, omega.e, beta.n, beta.star, solver)
      cone.n <- cone.n.lp(n, omega.i, beta.n, beta.star, lpmodel, 0, solver)
   }

   # ---------------- #
   # Step 5: Compute bootstrap components of cone.n and range.n
   # ---------------- #
   if (d >= p){
      # Compute the restricted estimator
      beta.r <- beta.r.compute(n, lpmodel, beta.obs.hat, beta.tgt, beta.n,
                               beta.star, omega.i, 1, solver)$x

      # Compute range.n for bootstrap beta
      range.n.list <- rep(0, R)
      
      cone.n.list <- list()
      # Compute cone.n for bootstrap beta
      for (i in 1:length(lambda)){
         cone.n.temp <- cone.n.bs(n, omega.i, beta.n, beta.star, lpmodel, 
                                  lambda[i], 1, beta.star.bs, beta.r, 
                                  beta.star.list, solver, cores, progress,
                                  length(lambda), i)
         cone.n.list[[i]] <- cone.n.temp
      }
      
   } else {
      # Compute the restricted estimator
      beta.r <- beta.r.compute(n, lpmodel, beta.obs.hat, beta.tgt, beta.n,
                               beta.star, omega.i, 0, solver)$x

      # Compute range.n for bootstrap beta
      range.n.list <- rangen.bs(n, omega.e, beta.n, beta.star, R, beta.n.bs, 
                                 beta.star.bs, solver, cores, progress)

      cone.n.list <- list()
      # Compute cone.n for bootstrap beta
      for (i in 1:length(lambda)){
         cone.n.temp <- cone.n.bs(n, omega.i, beta.n, beta.star, lpmodel, 
                                  lambda[i], 0, beta.star.bs, beta.r, 
                                  beta.star.list, solver, cores, progress,
                                  length(lambda), i)
         cone.n.list[[i]] <- cone.n.temp
      }
   }

   # ---------------- #
   # Step 6: Compute decision and p-value
   # ---------------- #
   pval <- data.frame(matrix(vector(), nrow = length(lambda), ncol = 2))
   colnames(pval) <- c("lambda", "p-value")
   for (i in 1:length(lambda)){
      pval.return <- fsst.pval(range.n$objval, cone.n$objval, range.n.list, 
                               cone.n.list[[i]], R, alpha)  
      pval[i,1] <- lambda[i]
      pval[i,2] <- pval.return$pval
   }

   # ---------------- #
   # Step 7: Close the progress bar that is used in the bootstrap procedure
   # ---------------- #
   if (progress == TRUE){
      cat("\n\b                                                ")
   }
   
   # ---------------- #
   # Step 8: Assign the return list and return output
   # ---------------- #
   # Assign the list of objects returned
   output <- list(pval = pval, 
                  cores = cores,
                  call = call,
                  range = range.n,
                  cone = cone.n,
                  test = max(range.n$objval, cone.n$objval),
                  cone.n.list = cone.n.list,
                  range.n.list = range.n.list,
                  solver = solver.name)
   
   # Assign class
   attr(output, "class") <- "fsst"
   
   # Return output
   return(output)
}

#' Computing bootstrap and variance estimators
#' 
#' @description Based on the function \code{beta.obs} from the object 
#'   \code{lpmodel}, this function computes the bootstrap estimates and 
#'   the asympototic variance estimator of the boostrap estimator. The 
#'   asymptotic variance estimator is computed as
#'     \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)'}
#'   for some suitable estimators \eqn{\widehat{\bm{\beta}}} and the 
#'   corresponding bootstrap estimators 
#'   \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.
#'   
#' @inheritParams fsst
#' @param beta.obs.hat Estimator of \eqn{\widehat{\bm{\beta}_{\mathrm{obs}, n}}}
#'   based on the given information in \code{lpmodel}.
#'   
#' @return Returns the bootstrap estimators and the estimator of the 
#'   asymptotic variance.
#'     \item{sigma.mat}{Estimator of the asympototic variance}
#'     \item{beta.obs.bs}{List of bootstrap betas}
#' 
#' @export
#' 
sigma.est <- function(n, data, beta.obs.hat, lpmodel, R){
   
   # ---------------- #
   # Step 1: Initialization
   # ---------------- #
   sigma.sum <- matrix(rep(0, length(beta.obs.hat)^2), 
                       nrow = length(beta.obs.hat))
   beta.obs.bs <- list()
   for (i in 1:R){
      # ---------------- #
      # Step 2: Obtain the bootstrap beta
      # ---------------- #
      # Obtain the bootstrap  beta
      data.bs <- as.data.frame(data[sample(1:nrow(data), replace = TRUE),])
      rownames(data.bs) <- 1:nrow(data.bs)
      beta.obs.bs[[i]] <- lpmodel$beta.obs(data.bs)
      beta.product <- (beta.obs.bs[[i]] - beta.obs.hat)
      
      # ---------------- #
      # Step 3: Compute the inner product
      # ---------------- #
      sigma.sum <- sigma.sum + t(beta.product) %*% (beta.product)
   }
   
   # Compute sigma hat
   sigma.hat <- sigma.sum*n/R
   return(list(sigma.hat = sigma.hat,
               beta.obs.bs = beta.obs.bs))
}

#' Computing bootstrap and variance estimators with parallel programming
#' 
#' @description Based on the function \code{beta.obs} from the object 
#'   \code{lpmodel}, this function computes the bootstrap estimates and 
#'   the asympototic variance estimator of the boostrap estimator. The 
#'   asymptotic variance estimator is computed as
#'     \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)'}
#'   for some suitable estimators \eqn{\widehat{\bm{\beta}}} and the 
#'   corresponding bootstrap estimators 
#'   \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.
#'   
#' @inheritParams fsst
#' @param beta.obs.hat Estimator of \eqn{\widehat{\bm{\beta}_{\mathrm{obs}, n}}}
#'   based on the given information in \code{lpmodel}.
#'   
#' @return Returns the bootstrap estimators and the estimator of the 
#'   asymptotic variance.
#'     \item{sigma.mat}{Estimator of the asympototic variance}
#'     \item{beta.obs.bs}{List of bootstrap betas}
#'     \item{pb}{Progress bar object}
#' 
#' @export
#' 
sigma.est.parallel <- function(data, beta.obs.hat, lpmode, R, cores,
                               progress){
   # ---------------- #
   # Step 1: Register the number of cores and extract information
   # ---------------- #
   options(warn = -1)
   
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
      pb <- utils::txtProgressBar(max = R, initial = 0, style = 3, width = 20) 
      cat("\r")
      progress <- function(n){
         utils::setTxtProgressBar(pb, n/10) 
         if (n < R){
            cat("\r\r") 
         } else {
            cat("\r\r")     
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
   listans <- foreach(i = 1:R, 
                      .multicombine = TRUE, 
                      .combine = "comb",
                      .options.snow = opts,
                      .packages = "lpinfer") %dorng% {
                         
       lpmodel.bs <- lpmodel
       
       # Re-sample the data 
       data.bs <- as.data.frame(Momocs::sample_n(data, 
                                                 nrow(data), 
                                                 replace = TRUE))
       rownames(data.bs) <- 1:nrow(data.bs)
       
       # Compute the bootstrap test statistic
       beta.obs.bs <- lpmodel$beta.obs(data.bs)
       beta.product <- (beta.obs.bs - beta.obs.hat)
       sigma.mat <- t(beta.product) %*% (beta.product)
       
       list(sigma.mat)
    }
   
   # ---------------- #
   # Step 4: Retrieve information from the output
   # ---------------- #
   sigma.hat <- Reduce('+', listans[[1]])*n/R
   
   # ---------------- #
   # Step 5: Return results
   # ---------------- #
   return(list(sigma.hat = sigma.hat,
               beta.obs.bs = beta.obs.bs))
}

#' Construct the full beta vector
#' 
#' @descripton This function concatenate the \code{beta.obs}, \code{beta.shp}
#'   and \code{beta.tgt} components to form the full \code{beta} vector.
#'   
#' @inheritParams fsst
#' @param beta.obs.bs List of bootstrap components of \code{beta.obs}.
#' 
#' @return List of full beta vectors
#'   \item{beta.bs}{List of full beta vector 
#'   \eqn{\{\widehat{\bm{\beta}_{n,b}}}\}^B_{b=1}}
#' 
#' @export
#' 
full.beta.bs <- function(lpmodel, beta.tgt, beta.obs.bs, R){
   beta.bs <- list()
   for (i in 1:R){
      beta.bs[[i]] <- c(beta.obs.bs[[i]], lpmodel$beta.shp, beta.tgt)
   }
   return(beta.bs)
}

#' Computing the variance estimator
#' 
#' @description Based on the bootstrap estimators
#'   \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}, this function computes the 
#'   asympototic variance estimator of the boostrap estimator, i.e.
#'     \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)'.}
#'   
#' @inheritParams fsst
#' @param beta.bs.list List of bootstrap estimators  
#'    \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.
#'   
#' @return Returns the bootstrap estimators and the estimator of the 
#'   asymptotic variance.
#'     \item{sigma.mat}{Estimator of the asympototic variance}
#' 
#' @export
#' 
sigma.summation <- function(n, beta.bs.list){
   # ---------------- #
   # Step 1: Estimate beta.obs.hat
   # ---------------- #
   beta.obs.hat <- as.matrix(beta.bs.list[[1]])
   beta.dim <- nrow(beta.obs.hat) * ncol(beta.obs.hat)
   
   # ---------------- #
   # Step 2: Obtain the variance
   # ---------------- #
   sigma.sum <- matrix(rep(0, beta.dim^2), nrow = beta.dim)
   for (i in 1:(length(beta.bs.list)-1)){
      beta.diff <- as.matrix(beta.bs.list[[i+1]] - beta.obs.hat)
      if (nrow(beta.diff) == 1){
         sigma.sum <- sigma.sum + t(beta.diff) %*% (beta.diff)
      } else {
         sigma.sum <- sigma.sum + (beta.diff) %*% t(beta.diff)
      }
   }
   sigma.hat <- sigma.sum/(length(beta.bs.list)-1)*n
   return(sigma.hat)
}

#' Computing the variance estimator using parallel programming
#' 
#' @description Based on the bootstrap estimators
#'   \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}, this function computes the 
#'   asympototic variance estimator of the boostrap estimator, i.e.
#'     \deqn{\frac{n}{B} \sum^B_{i=1} \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)  \left(\widehat{\bm{\beta}}_b - 
#'     \widehat{\bm{\beta}}\right)'.}
#'   
#' @inheritParams fsst
#' @param beta.bs.list List of bootstrap estimators  
#'    \eqn{\{\widehat{\bm{\beta}}_b\}^B_{b=1}}.
#' @param ind.times Indicator variable showing whether this is called for the 
#'    first or the second time
#'   
#' @return Returns the bootstrap estimators and the estimator of the 
#'   asymptotic variance.
#'     \item{sigma.mat}{Estimator of the asympototic variance}
#' 
#' @export
#' 
sigma.summation.parallel <- function(n, beta.bs.list, cores, progress, 
                                     ind.times){
   # ---------------- #
   # Step 1: Register the number of cores and extract information
   # ---------------- #
   options(warn = -1)
   
   # Register core
   doMC::registerDoMC(cores)
   R <- length(beta.bs.list) - 1
   
   # ---------------- #
   # Step 2: Initialize progress bar, comb function and assign doRNG
   # ---------------- #
   if (progress == TRUE){
      # Initialize the counter
      cl <- PtProcess::makeSOCKcluster(8)
      doSNOW::registerDoSNOW(cl)
      
      # Set the counter and progress bar
      if (ind.times == 1){
         pb <- utils::txtProgressBar(max = R, initial = 0, style = 3, 
                                     width = 20) 
         initial.bar <- 0
      } else if (ind.times == 2){
         pb <- utils::txtProgressBar(max = R, initial = R/10, style = 3, 
                                     width = 20) 
         initial.bar <- R/10
      }
      cat("\r")
      progress <- function(n){
         utils::setTxtProgressBar(pb, n/10 + initial.bar) 
         if (n < R){
            cat("\r\r") 
         } else {
            cat("\r\r")     
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
   # Step 3: Estimate beta.obs.hat
   # ---------------- #
   beta.obs.hat <- as.matrix(beta.bs.list[[1]])
   beta.dim <- nrow(beta.obs.hat) * ncol(beta.obs.hat)
   
   
   p <- length(beta.bs.list[[1]])
   
   # ---------------- #
   # Step 4: Bootstrap procedure
   # ---------------- #
   listans <- foreach(i = 1:R, 
                      .multicombine = TRUE, 
                      .combine = "comb", 
                      .options.snow = opts, 
                      .packages = "lpinfer") %dorng% {

      beta.diff <- as.matrix(beta.bs.list[[i+1]] - beta.obs.hat)
      if (nrow(beta.diff) == 1){
         listans <- t(beta.diff) %*% (beta.diff)
      } else {
         listans <- (beta.diff) %*% t(beta.diff)
      }
      listans <- as.matrix(listans, nrow = p)
      list(x=listans)
   }
   
   # ---------------- #
   # Step 5: Retrieve information from the output
   # ---------------- #
   sigma.hat1 <- matrix(unlist(listans[[1]][1:(p^2)]), nrow = p, byrow = TRUE)
   sigma.hat2 <- Reduce('+', listans[[1]][(p^2+1):(R+p^2-1)])
   sigma.hat <- (sigma.hat1 + sigma.hat2)*n/R
   
   return(sigma.hat)
}

#' Computing \eqn{\widehat{\bm{\beta}}_n^\star}
#' 
#' @description This function computes the vector 
#'   \eqn{\widehat{\bm{\beta}}_n^\star} for the case where \eqn{d<p}.
#'   
#' @inheritParams fsst
#' @inheritParams sigma.est
#' 
#' @details This corresponding to solving problem (3) of Torgovitsky (2020).
#'   Three options for \code{weight.matrix} are available:
#'   \itemize{
#'     \item{\code{identity} --- identity matrix}
#'     \item{\code{diag} --- the diagonal matrix that takes the diagonal 
#'        elements of the inverse of the variance matrix}
#'     \item{\code{avar} --- inverse of the variance matrix}
#' }
#' 
#' @return Returns the vector \eqn{\widehat{\bm{\beta}}_n^\star}.
#'   \item{beta.star}{The vector \eqn{\widehat{\bm{\beta}}_n^\star}.}
#' 
#' @export
#' 
beta.star.qp <- function(data, lpmodel, beta.tgt, weight.matrix, beta.obs.hat, 
                            beta.sigma, solver){
   # ---------------- #
   # Step 1: Choose the weighting matrix
   # ---------------- #
   weight.matrix <- tolower(weight.matrix)
   
   if (weight.matrix == "identity"){
      weight.mat = diag(nrow(as.matrix(beta.obs.hat)))
   } else if (weight.matrix == "diag"){
      weight.mat <- diag(diag(solve(beta.sigma)))
   } else if (weight.matrix == "avar"){
      weight.mat <- solve(beta.sigma)
   } else {
      stop("'weight.matrix' has to be one of 'identity', 'diag' and 'avar'.")
   }
   
   # ---------------- #
   # Step 2: Solve the quadratic program
   # ---------------- #
   # Define the parameters
   d <- ncol(as.matrix(lpmodel$A.tgt))
   
   # Constraints matrix
   A.mat1 <- as.matrix(lpmodel$A.shp)
   A.mat2 <- as.matrix(lpmodel$A.tgt)
   A.mat <- as.matrix(rbind(A.mat1, A.mat2))
   
   # RHS matrix
   b <- as.matrix(cbind(lpmodel$beta.shp, c(beta.tgt)))
   
   # Objective function
   A.obj <- as.matrix(cbind(lpmodel$A.obs))
   
   # Assign the optimization argument
   optim.arg <- list(Af = as.matrix(A.obj),
                     bf = as.matrix(beta.obs.hat),
                     nf = 1,
                     A = A.mat,
                     rhs = b,
                     sense = "=",
                     modelsense = "min",
                     lb = rep(-Inf, d),
                     weight = weight.mat)
   
   # Solve the model
   ans <- do.call(solver, optim.arg)
   
   # ---------------- #
   # Step 3: Compute beta.star
   # ---------------- #
   # Compute x.star
   x.star <- ans$x
   
   # Compute beta.star
   A <- as.matrix(rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt))
   beta.star <- A %*% x.star
   return(beta.star)
}

#' Computes the solution to the range problem
#' 
#' @description This function computes the solution to the range problem.
#' 
#' @import Matrix
#' 
#' @inheritParams fsst
#' @param omega.e The matrix \eqn{\widehat{\bm{\Omega}}^e_n}, i.e. the 
#'   regularized matrix for 
#'   \eqn{\widehat{\bm{\Sigma}}^{\beta^\star}_{n,\bar{\rho}}}.
#' 
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#' 
#' @export
#' 
rangen.lp <- function(n, omega.e, beta.n, beta.star, solver){
   # ---------------- #
   # Step 1: Compute the norms
   # ---------------- #
   # beta.norm
   beta.norm <- base::norm(beta.n - beta.star, type = "2")
   
   # Operator norm
   omega.operator <- Matrix::norm(omega.e, type = "2")
   
   # ---------------- #
   # Step 2: Construct the linear program
   # ---------------- #
   # Predefine a vector of zeros and ones
   p <- length(beta.n)
   ones.p <- rep(1, p)
   zero.p <- rep(0, p)
   
   # Update the objective function
   beta.diff <- (beta.n - beta.star)/beta.norm
   obj <- c(beta.diff, zero.p, zero.p)
   
   # Construct the lower bound
   lb <- c(rep(-Inf, p), rep(0, 2*p))
   
   # constraints matrix
   A.mat1 <- cbind(omega.e/omega.operator, -diag(p), diag(p))
   A.mat2 <- c(zero.p, ones.p, ones.p)
   
   A.mat <- rbind(A.mat1, A.mat2)
   rhs.mat <- c(zero.p, 1/omega.operator)
   
   # Sense
   sense.mat <- c(rep("=", p), rep("<=", 1))
   
   # Construct the arguments
   optim.arg <- list(Af = NULL,
                     bf = obj,
                     nf = 1,
                     A = A.mat,
                     rhs = rhs.mat,
                     sense = sense.mat,
                     modelsense = "max",
                     lb = lb)
   
   # ---------------- #
   # Step 3: Solve the quadartic program
   # ---------------- #
   ans <- do.call(solver, optim.arg)
   
   objval <- ans$objval * sqrt(n) * beta.norm
   
   return(list(x = ans$x,
               objval = objval))
}

#' Computes the solution to the cone problem
#' 
#' @description This function computes the solution to the cone problem.
#' 
#' @import Matrix
#' 
#' @inheritParams fsst
#' @param omega.i The matrix \eqn{\widehat{\bm{\Omega}}^i_n}, i.e. the 
#'   regularized matrix for 
#'   \eqn{\widehat{\bm{\Sigma}}^{\beta^\star}_{n,\bar{\rho}}}.
#' @param indicator A binary variable that equals to 1 for \eqn{d \geq p} and 
#'   equals to 0 for \eqn{d < p}.
#'   
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#' 
#' @export 
#' 
cone.n.lp <- function(n, omega.i, beta.n, beta.star, lpmodel, indicator,
                      solver){
   # ---------------- #
   # Step 1: Construct the linear program
   # ---------------- #
   p <- length(beta.n)
   d <- ncol(lpmodel$A.obs)
   ones.p <- rep(1, p)
   zero.p <- rep(0, p)
   zero.d <- rep(0, d)
   zero.dp <- matrix(rep(0, d*p), nrow = d)
   
   # Update the objective function
   obj <- c(beta.star, zero.p, zero.p)
   
   # Construct the lower bound
   lb <- c(rep(-Inf, p), rep(0, 2*p))
   
   # Construct the constraints matrix
   A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
   ones <- rep(1, p)
   
   A.mat1 <- cbind(omega.i, -diag(p), diag(p))
   A.mat2 <- c(zero.p, ones.p, ones.p)
   A.mat3 <- cbind(t(A), zero.dp, zero.dp)
   A.mat <- rbind(A.mat1, A.mat2, A.mat3)
   
   # Construct RHS vector
   rhs.mat <- c(zero.p, 1, zero.d)
   
   # Sense
   sense.mat <- c(rep("=", p), rep("<=", d+1))
   
   # Construct the arguments
   if (indicator == 1){
      optim.arg <- list(Af = NULL,
                        bf = obj,
                        nf = 1,
                        A = A.mat,
                        rhs = rhs.mat,
                        sense = sense.mat,
                        modelsense = "max",
                        lb = lb)
   } else {
      zero.pp <- matrix(rep(0, p*p), nrow = p)
      zero.Am <- matrix(rep(0, (p+d+1)*(d)), nrow = (p+d+1))
      
      # Update objective function
      obj.ext <- c(obj, zero.d)
      
      # Update constraints matrix
      A.mat.ext1 <- as.matrix(cbind(A.mat, zero.Am))
      A.mat.ext2 <- as.matrix(cbind(diag(p), zero.pp, zero.pp, -A))
      A.mat.ext <- rbind(A.mat.ext1, A.mat.ext2)
      
      # Update RHS vector
      rhs.ext <- c(rhs.mat, zero.p)
      sense.ext <- c(sense.mat, rep("=", p))
      
      # Update lower bound
      lb.ext <- c(lb, rep(-Inf,d))
      
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
               objval = objval))
}

#' Computes \eqn{\widehat{\bm{\beta}}^r_n}
#' 
#' @description This function computes the restricted estimator 
#'   \eqn{\widehat{\bm{\beta}}^r_n}.
#' 
#' @inheritParams fsst
#' @inheritParams cone.n.lp
#'   
#' @return Returns the optimal point and optimal value.
#'  \item{x}{Optimal point calculated from the optimizer.}
#'  \item{objval}{Optimal value calculated from the optimizer.}
#' 
#' @export 
#' 
beta.r.compute <- function(n, lpmodel, beta.obs.hat, beta.tgt, beta.n, 
                           beta.star, omega.i, indicator, solver){
   # ---------------- #
   # Step 1: Construct the linear program
   # ---------------- #
   # Construct the zero matrices and vectors and a vector of ones
   p <- length(beta.n)
   q <- nrow(beta.obs.hat)
   d <- ncol(lpmodel$A.obs)
   ones.1p <- rep(1, p)
   ones.p1 <- matrix(ones.1p, nrow = p)
   zero.1p <- rep(0, p)
   zero.p1 <- matrix(zero.1p, nrow = p)
   zero.1d <- rep(0, d)
   zero.pd <- matrix(rep(0, p*d), nrow = p)
   zero.pp <- matrix(rep(0, p*p), nrow = p)
   zero.qq <- matrix(rep(0, q*q), nrow = q)
   zero.pqq <- matrix(rep(0, (p-q)*q), nrow = q)
   iden.beta <- rbind(cbind(zero.qq, zero.pqq), 
                      cbind(t(zero.pqq), diag(p-q)))
   
   # Construct the constraints matrix
   A <- rbind(lpmodel$A.obs, lpmodel$A.shp, lpmodel$A.tgt)
   A.mat1 <- as.matrix(cbind(sqrt(n)*diag(p), zero.pd, -omega.i, A, zero.p1))
   if (indicator == 0){
      # Multiply A.mat 1 by t(A) if indicator == 0 (i.e. d < p)
      A.mat1 <- t(A) %*% A.mat1
   }
   A.mat2 <- as.matrix(cbind(diag(p), -A, zero.pp, zero.pd, zero.p1))
   A.mat3 <- as.matrix(cbind(zero.pp, zero.pd, -diag(p), zero.pd, ones.p1))
   A.mat4 <- as.matrix(cbind(zero.pp, zero.pd, diag(p), zero.pd, ones.p1))
   A.mat5 <- as.matrix(cbind(iden.beta, zero.pd, zero.pp, zero.pd, zero.p1))

   # Construct the rhs vector and the sense vector
   A.mat <- rbind(A.mat1, A.mat2, A.mat3, A.mat4, A.mat5)
   if (indicator == 0){
      rhs.mat <- c(sqrt(n) * t(A) %*% beta.star, zero.1p, zero.1p, zero.1p, 
                   rep(0,q), lpmodel$beta.shp, beta.tgt)
      sense.mat <- c(rep("=", nrow(A.mat1) + nrow(A.mat2)), 
                     rep(">=", 2*p), rep("=", p))
   } else {
      rhs.mat <- c(sqrt(n) * beta.n, zero.1p, zero.1p, zero.1p, rep(0,q), 
                   lpmodel$beta.shp, beta.tgt)
      sense.mat <- c(rep("=", 2*p), rep(">=", 2*p), rep("=", p))
   }
   
   # Construct the objective function
   obj <- c(rep(0, 2*(p+d)), 1)
   
   # Lower bound
   lb <- c(rep(-Inf,p), rep(0, d), rep(-Inf, p), rep(0, d+1))
   
   # Set the arguments
   optim.arg <- list(Af = NULL,
                     bf = obj,
                     nf = 1,
                     A = A.mat,
                     rhs = rhs.mat,
                     sense = sense.mat,
                     modelsense = "min",
                     lb = lb)
   
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

#' Computing bootstrap components of `Range` 
#' 
#' @description This function computes the bootstrap components of the 
#'   `Range` test statistics.
#'   
#' @inheritParams fsst
#' @inheritParams rangen.lp
#' 
#' @return A list of bootstrap Range statistics.
#'   \item{range.n.list}{A list of bootstrap Range statistics 
#'     \eqn{\{\mathrm{Range}_{n,b}\}^B_{b=1}}}
#'     
#' @export
#' 
rangen.bs <- function(n, omega.e, beta.n, beta.star, R, beta.n.bs, 
                       beta.star.bs, solver, cores, progress){
   range.n.list <- NULL
   
   if (cores == 1){
      for (i in 1:R){
         # ---------------- #
         # Step 1: Compute the replacements
         # ---------------- #
         beta.bs.1 <- beta.n.bs[[i]] - beta.n
         beta.bs.2 <- beta.star.bs[[i]] - beta.star
         
         # ---------------- #
         # Step 2: Solve the linear program and extract the solution
         # ---------------- #
         range.n.return <- rangen.lp(n, omega.e, beta.bs.1, beta.bs.2, 
                                         solver)
         range.n.list <- c(range.n.list, range.n.return$objval)
      } 
   } else {
      if (progress == TRUE){
         # Initialize the counter
         cl <- PtProcess::makeSOCKcluster(8)
         doSNOW::registerDoSNOW(cl)
         
         # Set the counter and progress bar
         pb <- utils::txtProgressBar(max = R, initial = R/5, style = 3, 
                                     width = 20) 
         cat("\r")
         progress <- function(n){
            utils::setTxtProgressBar(pb, (4*n)/10+R/5)
            if (n < R){
               cat("\r\r") 
            } else {
               cat("\r\r")     
            }
         }
         
         opts <- list(progress = progress) 
      } else {
         pb <- NULL
         opts <- NULL
      }
      
      # Assign doRnG
      `%dorng%` <- doRNG::`%dorng%`
      
      listans <- foreach(i = 1:R,
                         .multicombine = TRUE,
                         .options.snow = opts, 
                         .packages = "lpinfer") %dorng% {
                            
       # ---------------- #
       # Step 1: Compute the replacements
       # ---------------- #
       beta.bs.1 <- beta.n.bs[[i]] - beta.n
       beta.bs.2 <- beta.star.bs[[i]] - beta.star     
       
       # ---------------- #
       # Step 2: Solve the linear program and extract the solution
       # ---------------- #
       range.n.return <- rangen.lp(n, omega.e, beta.bs.1, beta.bs.2, 
                                       solver)
       list(range.n.return$objval)
      }
      range.n.list <- unlist(listans)
   }
   
   return(range.n.list)
   
}

#' Computing bootstrap components of `Cone` 
#' 
#' @description This function computes the bootstrap components of the 
#'   `Cone` test statistics.
#'   
#' @inheritParams fsst
#' @inheritParams cone.n.lp
#' 
#' @return A list of bootstrap Range statistics.
#'   \item{range.n.list}{A list of bootstrap Range statistics 
#'     \eqn{\{\mathrm{Cone}_{n,b}\}^B_{b=1}}}
#'     
#' @export
#' 
cone.n.bs <- function(n, omega.i, beta.n, beta.star, lpmodel, lambda,
                      indicator, beta.star.bs, beta.r, beta.star.list, 
                      solver, cores, progress, length.lambda, lambda.i){
   cone.n.list <- NULL
   
   if (cores == 1){
      for (i in 1:R){
         # ---------------- #
         # Step 1: Compute the replacements
         # ---------------- #
         beta.new <- beta.star.list[[i+1]] - beta.star + lambda*beta.r
         
         # ---------------- #
         # Step 2: Solve the linear program and extract the solution
         # ---------------- #
         cone.n.return <- cone.n.lp(n, omega.i, beta.n, beta.new, lpmodel, 
                                       indicator, solver)
         cone.n.list <- c(cone.n.list, cone.n.return$objval)
      }
   } else {
      options(warn=-1)
      
      # Register core
      doMC::registerDoMC(cores)
      
      if (progress == TRUE){
         
         # Initialize the counter
         cl <- PtProcess::makeSOCKcluster(8)
         doSNOW::registerDoSNOW(cl)
         
         # Set the counter and progress bar
         lambda.bar.i0 <- .4*(lambda.i-1)/length.lambda
         pb <- utils::txtProgressBar(max = R, 
                                     initial = R*.6 + R*lambda.bar.i0,
                                     style = 3, 
                                     width = 20) 
         cat("\r")
         progress <- function(n){
            utils::setTxtProgressBar(pb, R*.6 + R*lambda.bar.i0 + 
                                         n*.4/length.lambda) 
            if (n < R){
               cat("\r\r") 
            } else {
               cat("\r\r")     
            }
         }
         
         opts <- list(progress = progress) 
      } else {
         pb <- NULL
         opts <- NULL
      }
      
      # Assign doRnG
      `%dorng%` <- doRNG::`%dorng%`
      
      listans <- foreach(i = 1:R,
                         .multicombine = TRUE,
                         .options.snow = opts, 
                         .packages = "lpinfer") %dorng% {
         # ---------------- #
         # Step 1: Compute the replacements
         # ---------------- #
         beta.new <- beta.star.list[[i+1]] - beta.star + lambda*beta.r 
          
         # ---------------- #
         # Step 2: Solve the linear program and extract the solution
         # ---------------- #
         cone.n.return <- cone.n.lp(n, omega.i, beta.n, beta.new, lpmodel, 
                                       indicator, solver)
         list(cone.n.return$objval)
      }
      cone.n.list <- unlist(listans)
   }
   
   return(cone.n.list)
}

#' Auxiliary function to calculate the p-value of `\code{fsst}`
#'
#' @description This function computes the \eqn{p}-value of the test based on
#'    the bootstrap estimates.
#'
#' @param range.n The point estimate for the `Range` test statistics.
#' @param cone.n The point estimate for the `Cone` test statistics.
#' @param range.n.list The bootstrap test estimates for the `Range` test 
#'   statistics.
#' @param cone.n.list The bootstrap test estimates for the `Cone` test 
#'   statistics.
#' @inheritParams fsst
#'
#' @return Returns the \eqn{p}-value and the decision.
#'   \item{pval}{\eqn{p}-value.}
#'   \item{decision}{Decision to reject or not.}
#'
#' @export
#' 
fsst.pval <- function(range.n, cone.n, range.n.list, cone.n.list, R, 
                      alpha = .05){
   # ---------------- #
   # Step 1: Compute the test statistic and the bootstrap test statistics
   # ---------------- #
   T.n <- max(range.n, cone.n)
   T.bs <- NULL
   for (i in 1:R){
      T.bs <- c(T.bs, max(range.n.list[[i]], cone.n.list[[i]]))
   }
   
   # ---------------- #
   # Step 2: Decision
   # ---------------- #
   # Compute p-value
   pval <- mean(T.bs >= T.n)
   
   # Decision
   if (pval > alpha){
      decision <- 1
   } else {
      decision <- 0
   }
   return(list(pval = pval,
               decision = decision))
}

#' Checks and updates the input
#'
#' @description This function checks and updates the input of the user. If 
#'    there is any invalid input, this function will be terminated and 
#'    generates appropriate error messages.
#'
#' @inheritParams fsst
#' 
#' @return Returns the list of updated parameters as follows:
#'   \item{data}{Upated data in class \code{data.frame}}
#'   \item{lpmodel}{A list of linear programming objects in this 
#'      `\code{lpinfer}` package.}
#'   \item{solver}{Solver to be used.}
#'   \item{solver.name}{Updated name of solver in lower case.}
#'   \item{cores}{Updated number of cores to be used in the parallelization
#'      of the for-loop in the bootstrap procedure.}
#' 
#' @export
#' 
fsst.check <- function(data, lpmodel, beta.tgt, R, lambda, rho, n,
                       weight.matrix, solver, cores, progress){
   
   # ---------------- #
   # Step 1: Check the arguments
   # ---------------- #
   # Check data
   if (is.null(n)){
      data <- check.dataframe(data, "data")  
   } else {
      check.positiveinteger(n)
   }
   
   # Check lpmodel 
   lpmodel <- check.lpmodel(data = data,
                            lpmodel = lpmodel, 
                            name.var = "lpmodel",
                            A.tgt.cat = 1,
                            A.obs.cat = 1,
                            A.shp.cat = 1,
                            beta.obs.cat = c(2,3),
                            beta.shp.cat = 1,
                            R = R)
   # Check solver
   solver.return <- check.solver(solver, "solver")
   solver <- solver.return$solver
   solver.name <- solver.return$solver.name
   
   # Check numerics 
   check.numeric(beta.tgt, "beta.tgt")
   check.positiveinteger(R, "R")
   cores <- check.positiveinteger(cores, "cores")
   
   # Check lambda - can be vector of scalar
   if (length(lambda) == 1){
      check.numeric(lambda, "lambda")
   } else {
      if (class(lambda) != "numeric"){
         stop("The class of the variable 'lambda' has to be numeric.",
              call. = FALSE)
      }
   }
   
   # ---------------- #
   # Step 2: Return results
   # ---------------- #
   return(list(solver = solver,
               solver.name = solver.name,
               cores = cores))
}

#' Print results from \code{fsst}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{fsst}.
#'    
#' @param x Object returned from \code{fsst}.
#' @param ... Additional arguments.
#' 
#' @details The following information are printed
#'   \itemize{
#'     \item{Test statistic (Max of Cone and Range)}
#'     \item{Test statistic (Cone)}
#'     \item{Test statistic (Range)}
#'     \item{\eqn{p}-value}
#'     \item{Solver used}
#'     \item{Number of cores used}
#'   }
#'    
#' @export
#'    
print.fsst <- function(x, ...){
   cat("\r\r")
   cat(sprintf("Test statistic: %s.             \n", round(x$test, digits = 5)))
   cat(sprintf("   - Range component: %s\n", round(x$range$objval, digits = 5)))  
   cat(sprintf("   - Cone component: %s\n", round(x$cone$objval, digits = 5)))
   if (nrow(x$pval) == 1){
      cat(sprintf("p-value: %s\n", x$pval[1,2]))
   } else {
      cat("p-values:\n")
      cat("     lambda\tp-value\n")
      for (i in 1:nrow(x$pval)){
         cat(sprintf("     %s\t%s\n", 
                     round(x$pval[i,1], digits = 5),
                     round(x$pval[i,2], digits = 5)))
      }
   }
   cat(sprintf("Solver used: %s\n", x$solver.name))
   cat(sprintf("Number of cores used: %s\n", x$cores))
}

#' Summary of results from \code{fsst}
#' 
#' @description This function uses the summary method on the return list of the
#'    function \code{fsst}. This is a wrapper of the \code{print} command.
#'    
#' @param x Object returned from \code{fsst}.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{fsst}.
#' 
#' @export
#' 
summary.fsst <- function(x, ...){
   print(x)
}