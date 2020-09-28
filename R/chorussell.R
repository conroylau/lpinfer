#' Conducts inference using the Cho-Russell procedure
#'
#' @description This function conducts inference using the procedure
#'   proposed by Cho and Russell (2019).
#'
#' @inheritParams estbounds
#' @inheritParams dkqs
#' @param ci A boolean variable that indicates whether a \eqn{p}-value or a
#'   \eqn{(1-\alpha)}-confidence interval is returned. If \code{ci} is
#'   \code{TRUE}, then a confidence interval is returned. Otherwise, a
#'   \eqn{p}-value is returned.
#' @param alpha The significance level. This can be a vector.
#' @param tol The tolerance level in the bisection procedure.
#' @param kappa The tuning parameter used in the second step of the two-step
#'    procedure for obtaining the bounds subject to the shape constraints.
#'    It can be any nonnegative number or a vector of nonnegative numbers.
#' @param remove.const A boolean variable. This determine whether the
#'    constraints are to be removed.
#'
#' @return Returns the following list of output:
#'   \item{ub}{The upper bound from original data.}
#'   \item{lb}{The lower bound from original data.}
#'   \item{ub.bs}{The list of upper bounds from bootstrap data.}
#'   \item{lb.bs}{The list of lower bounds from bootstrap data.}
#'   \item{test.logical}{An indicator variable for whether the computation has
#'     been conducted. If \code{test.logical} is 1, it refers to the case
#'     where \code{beta.tgt} is inside the logical bound. If
#'     \code{test.logical} is 0, it refers to the case where
#'     \code{beta.tgt} is outside the logical bound.}
#'   \item{logical.lb}{The logical lower bound.}
#'   \item{logical.ub}{The logical upper bound.}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'   \item{ci}{A boolean variable that indicates whether a \eqn{p}-value or a
#'     \eqn{(1-\alpha)}-confidence interval is returned.}
#'   \item{pval}{\eqn{p}-value (if \code{ci} is set as \code{FALSE}).}
#'   \item{c.ub}{The upper bound of the \eqn{(1-\alpha)}-confidence interval.}
#'   \item{c.lb}{The lower bound of the \eqn{(1-\alpha)}-confidence interval.}
#'   \item{alpha}{The significance level.}
#'   \item{iter}{The total number of iterations (if \code{ci} is \code{FALSE}.)}
#'   \item{unique}{A boolean variable showing whether the solution is unique.}
#'
#' @details
#' \itemize{
#'  \item{See the details section of the \code{\link[lpinfer]{estbounds}}
#'     function for a list of strings acceptable for the option \code{norm}.}
#'  \item{The following components are required in the \code{lpmodel} for the
#'    Cho-Russell procedure:
#'    \itemize{
#'      \item{\code{A.tgt}}
#'      \item{\code{A.obs}}
#'      \item{\code{A.shp}}
#'      \item{\code{beta.obs}}
#'      \item{\code{beta.shp}}
#'    }
#'  \item{The input \code{beta.tgt} is not required when \code{ci = TRUE}.}
#'  }
#' }
#'
#' @example ./inst/example/chorussell_example.R
#'
#' @export
#'
chorussell <- function(data = NULL, lpmodel, beta.tgt = NULL, n = NULL, R = 100,
                       Rmulti = 1.25, kappa = 0, norm = 2, estimate = TRUE,
                       solver = NULL, ci = NULL, alpha = 0.05, tol = 1e-4,
                       progress = TRUE, remove.const = TRUE) {
  # ---------------- #
  # Step 1: Update call, check and update the arguments; initialize df.error
  # ---------------- #
  # Obtain call information
  call <- match.call()

  # Check the arguments
  chorussell.return <- chorussell.check(data, lpmodel, beta.tgt, R, Rmulti,
                                        kappa, norm, n, estimate, solver,
                                        ci, alpha, tol, progress)

  # Update the arguments
  ci <- chorussell.return$ci
  data <- chorussell.return$data
  solver <- chorussell.return$solver
  solver.name <- chorussell.return$solver.name
  test.logical <- chorussell.return$test.logical
  logical.lb <- chorussell.return$logical.lb
  logical.ub <- chorussell.return$logical.ub
  n <- chorussell.return$n

  # Compute the maximum number of iterations
  maxR <- ceiling(R * Rmulti)

  # Sort kappa
  kappa <- sort(kappa)

  # Sort alpha
  alpha <- sort(alpha)

  ### Case 1: test.logical == 1. Proceed with the calculation because
  ### beta.tgt is inside the logical bounds. test.logical also does not matter
  ### if the user would like to construct a confidence interval.
  if ((test.logical == 1) | (isTRUE(ci))) {
    # ---------------- #
    # Step 2: Obtain estimated bounds
    # ---------------- #
    # Initialization
    ub <- list()
    lb <- list()
    delta <- list()

    # Loop through each elements in the list
    for (i in seq_along(kappa)) {
      if (is.list(lpmodel$beta.obs)) {
        lpm <- lpmodel
        lpm$beta.obs <- lpmodel$beta.obs[[1]]
      } else {
        lpm <- lpmodel
      }
      estb.out <- estbounds(data, lpm, kappa[i], norm, estimate,
                            solver.name)
      ub[[i]] <- estb.out$ub
      lb[[i]] <- estb.out$lb
      delta[[i]] <- ub[[i]] - lb[[i]]
    }

    # ---------------- #
    # Step 3: Obtain estimated bounds from bootstrap data
    # ---------------- #
    cr.bs.return <- chorussell.bs(data, lpmodel, beta.tgt, R, maxR, kappa,
                                  norm, n, estimate, solver.name,
                                  progress)
    ub.bs <- cr.bs.return$ub
    lb.bs <- cr.bs.return$lb
    df.error <- cr.bs.return$df.error
    R.succ <- cr.bs.return$R.succ

    # ---------------- #
    # Step 4: Solve the optimization problem in Cho-Russell
    # ---------------- #
    # Initialization
    lb.can1 <- list()
    lb.can2 <- list()
    ub.can1 <- list()
    ub.can2 <- list()
    cr.lp.return <- list()

    ## Candidates for the upper and lower bounds. If the bound is inifinte, then
    ## the correspond (1 - alpha)-bound is also inifinte so the procedure of
    ## computing the candidates will be skipped.
    for (i in seq_along(kappa)) {
      # Lower bounds
      if (lb[[i]] != -Inf) {
        lb.can1[[i]] <- sqrt(n) * (lb.bs[[i]] - lb[[i]])
        lb.can2[[i]] <- sqrt(n) * (lb.bs[[i]] - lb[[i]] - delta[[i]])
      } else {
        lb.can1[[i]] <- NULL
        lb.can2[[i]] <- lb.can2[[i]]
      }

      # Upper bounds
      if (ub[[i]] != Inf) {
        ub.can1[[i]] <- sqrt(n) * (ub.bs[[i]] - ub[[i]])
        ub.can2[[i]] <- sqrt(n) * (ub.bs[[i]] - ub[[i]] + delta[[i]])
      } else {
        ub.can1[[i]] <- NULL
        ub.can2[[i]] <- ub.can1[[i]]
      }

      # Either find the p-value or the confidence interval
      if ((lb[[i]] == Inf) & (ub[[i]] == Inf)) {
        cr.lp.return[[i]] <- list()
        cr.lp.return[[i]]$pval <- 1
        cr.lp.return[[i]]$c.ub <- Inf
        cr.lp.return[[i]]$c.lb <- -Inf
        cr.lp.return[[i]]$unique <- TRUE
      } else {
        # Simplify the constraints if needed
        if (isTRUE(remove.const)) {
          simp.return <- chorussell.simp(lb.can1[[i]], lb.can2[[i]],
                                         ub.can1[[i]], ub.can2[[i]],
                                         progress)
          df.lb <- simp.return$df.lb
          df.ub <- simp.return$df.ub
        } else {
          df.lb <- NULL
          df.ub <- NULL
        }
        
        cr.lp.return[[i]] <- chorussell.eval(beta.tgt, lb.can1[[i]],
                                             lb.can2[[i]], ub.can1[[i]],
                                             ub.can2[[i]], n, R, ci, alpha,
                                             tol, ub[[i]], lb[[i]],
                                             logical.ub, logical.lb,
                                             remove.const, kappa[i],
                                             progress, df.lb, df.ub)
      }
    }

    # ---------------- #
    # Step 5; Consolidate the list of outputs
    # ---------------- #
    output <- list(ub = ub,
                   lb = lb,
                   ub.bs = ub.bs,
                   lb.bs = lb.bs,
                   test.logical = test.logical,
                   logical.lb = logical.lb,
                   logical.ub = logical.ub,
                   df.error = df.error,
                   ci = ci,
                   call = call,
                   norm = norm,
                   R.succ = R.succ,
                   kappa = kappa,
                   solver = solver.name)

    # Append the p-value or confidence interval based on ci
    if (isFALSE(ci)) {
      # p-values
      pval.df <- data.frame(matrix(vector(), nrow = length(kappa), ncol = 2))
      colnames(pval.df) <- c("kappa", "pvalue")
      for (i in seq_along(kappa)) {
        pval.df[i, 1] <- kappa[i]
        pval.df[i, 2] <- cr.lp.return[[i]]$pval
      }
      output$pval <- pval.df
    } else {
      # Confidence intervals
      ci.df <- data.frame(matrix(vector(),
                                 nrow = length(kappa) * length(alpha),
                                 ncol = 4))
      colnames(ci.df) <- c("alpha", "kappa", "lb", "ub")
      k <- 1
      for (i in seq_along(kappa)) {
        for (j in seq_along(alpha)) {
          ci.df[k, 1] <- alpha[j]
          ci.df[k, 2] <- kappa[i]
          ci.df[k, 3] <- as.numeric(cr.lp.return[[i]][[j]]$bd[1])
          ci.df[k, 4] <- as.numeric(cr.lp.return[[i]][[j]]$bd[2])
          k <- k + 1
        }
      }
      output$ci.df <- ci.df
      output$alpha <- alpha
    }
  } else {
    output <- list(pval = 0,
                   solver = solver.name,
                   call = call,
                   test.logical = test.logical)
  }

  # Turn returned lists for bounds as vectors if parameter is not multivalued
  if (length(kappa) == 1) {
    for (x in c("ub", "lb", "ub.bs", "lb.bs")) {
      output[[x]] <- unlist(output[[x]])
    }
  }

  # ---------------- #
  # Step 6: Issue warning message for estimate = FALSE
  # ---------------- #
  if (isFALSE(estimate)) {
    if (length(nrow(NULL) > 0)) {
      warning(sprintf(paste0("Among the %s bootstrap replications, ",
                             "%s failed when estimate is set as FALSE"),
                      R.succ + nrow(df.error), nrow(df.error)))
    }
  }

  # Assign the class of output
  attr(output, "class") = "chorussell"

  return(output)
}

#' Simplifies the candidates to be considered in
#' \code{\link[lpinfer]{chorussell}}
#'
#' @import future.apply progressr
#'
#' @description This function simplifies the list of candidates to be
#'   considered in the optimization problem in the
#'   \code{\link[lpinfer]{chorussell}} function. In particular, because
#'   \deqn{\mathbf{1}\left[\sqrt{n}\left(\hat{\theta}_{\rm lb}^b -
#'   \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb}\right]
#'   \geq
#'   \mathbf{1}\left[\sqrt{n}\left(\hat{\theta}^b_{\rm lb} -
#'   \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb} \quad \mathrm{ and } \quad
#'   -c_{\rm ub} \leq \sqrt{n} \left(\hat{\theta}^b_{\rm ub} -
#'   \hat{\theta}_{\rm ub} + \Delta\right)
#'   \right],}
#'   the values of \eqn{c_{\rm lb}} that satisfy
#'   \deqn{\frac{1}{B}\sum^B_{b=1} \mathbf{1}
#'   \left[\sqrt{n}\left(\hat{\theta}_{\rm lb}^b -
#'   \hat{\theta}_{\rm lb}\right) \leq c_{\rm lb}\right] < 1 - \alpha}
#'   will be removed from the set of the values under consideration. Similarly,
#'   the list of \eqn{c_{\rm ub}} that satisfy
#'   \deqn{\frac{1}{B}\sum^B_{b=1} \mathbf{1}
#'   \left[- \sqrt{n}\left(\hat{\theta}_{\rm ub}^b -
#'   \hat{\theta}_{\rm ub}\right) \leq c_{\rm lb}\right] < 1 - \alpha}
#'   will be removed from the set of the values under consideration.
#'
#' @param lb.can1 The vector of values that corresponds to
#'   \eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm lb} - \hat{\theta}_{\rm lb}\right)}.
#' @param lb.can2 The vector of values that corresponds to
#'   \eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm lb} - \hat{\theta}_{\rm lb} -
#'   \Delta \right)}.
#' @param ub.can1 The vector of values that corresponds to
#'   \eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm ub} - \hat{\theta}_{\rm ub}\right)}.
#' @param ub.can2 The vector of values that corresponds to
#'   \eqn{\sqrt{n}\left(\hat{\theta}^b_{\rm ub} - \hat{\theta}_{\rm ub} +
#'   \Delta \right)}.
#' @inheritParams chorussell
#'
#' @return Returns the list of updated candidates of lower bounds and upper
#'   bounds.
#'   \item{lb}{The updated list of candidates of lower bounds.}
#'   \item{ub}{The updated list of candidates of upper bounds.}
#'
#' @export
#'
chorussell.simp <- function(lb.can1, lb.can2, ub.can1, ub.can2, progress) {
  # ---------------- #
  # Step 1: Combine the candidates
  # ---------------- #
  lb.can <- c(lb.can1, lb.can2)
  ub.can <- -c(ub.can1, ub.can2)

  # ---------------- #
  # Step 2: Obtain the updated candidates
  # ---------------- #
  # Set the default for progress bar
  progressr::handlers("progress")

  # Lower bound
  progressr::with_progress({
    if (isTRUE(progress)) {
      pbar <- progressr::progressor(along = seq_along(lb.can))
    } else {
      pbar <- NULL
    }

    lb.return <- future.apply::future_lapply(lb.can,
                                             FUN = chorussell.simp.fn,
                                             can1 = lb.can1,
                                             pbar = pbar,
                                             bd = "lower",
                                             progress = progress)
    df.lb <- data.frame(bd = unlist(sapply(lb.return, "[", "bd"),
                                    use.names = FALSE),
                        sum = unlist(sapply(lb.return, "[", "sum"),
                                     use.names = FALSE))
  })

  # Upper bound
  progressr::with_progress({
    if (isTRUE(progress)) {
      pbar <- progressr::progressor(along = seq_along(ub.can))
    } else {
      pbar <- NULL
    }

    ub.return <- future.apply::future_lapply(ub.can,
                                             FUN = chorussell.simp.fn,
                                             can1 = ub.can1,
                                             pbar = pbar,
                                             bd = "upper",
                                             progress = progress)
    df.ub <- data.frame(bd = unlist(sapply(ub.return, "[", "bd"),
                                    use.names = FALSE),
                        sum = unlist(sapply(ub.return, "[", "sum"),
                                     use.names = FALSE))
  })

  return(list(df.lb = df.lb,
              df.ub = df.ub))
}

#' Checks one candidate in \code{\link[lpinfer]{chorussell}}
#'
#' @description This function checks one candidate for the
#'   \eqn{(1-\alpha)}-confidence interval to see whether it should be dropped
#'   from the list of considerations. For details, please see the description
#'   of the \code{\link[lpinfer]{chorussell.simp}} function.
#'
#' @param x One candidate of the upper or lower bound.
#' @param can1 This refers to either \code{lb.can1} or \code{ub.can1}.
#' @param bd A string indicating which set of candidates that the function is
#'   trying to simplify.
#' @inheritParams chorussell
#' @inheritParams dkqs.bs.fn
#'
#' @return Returns the decision to drop or not to drop the candidate. If the
#'   decision is to drop the candidate bound, \code{NULL} is returned.
#'   Otherwise, the bound is returned back via the \code{future_lapply}
#'   function.
#'
#' @export
#'
chorussell.simp.fn <- function(x, can1, pbar, bd, progress) {
  # Message for progress bar
  if (isTRUE(progress)) {
    pbar(sprintf("(Checking the candidates for the %s bound)", bd))
  }

  if (bd == "lower") {
    return(list(bd = x,
                sum = mean(can1 <= x)))
  } else if (bd == "upper") {
    return(list(bd = x,
                sum = mean(-x <= can1)))
  }
}

#' Bootstrap procedure for the \code{\link[lpinfer]{chorussell}} procedure
#'
#' @description This function carries out the bootstrap procedure of the
#'   \code{\link[lpinfer]{chorussell}} procedure. This function supports
#'   parallel programming via the \code{future.apply} package.
#'
#' @import future.apply progressr
#'
#' @inheritParams chorussell
#' @inheritParams dkqs.bs
#'
#' @return Returns a list of output that are obtained from the Cho-Russell
#'   procedure:
#'   \item{ub.bs}{The list of upper bounds from bootstrap data.}
#'   \item{lb.bs}{The list of lower bounds from bootstrap data.}
#'   \item{df.error}{A table showing the id of the bootstrap replication(s)
#'     with error(s) and the corresponding error message(s).}
#'   \item{error.list}{A list of error messages.}
#'   \item{R.eval}{The number of bootstrap replications that have been
#'     conducted.}
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'
#' @export
#'
chorussell.bs <- function(data, lpmodel, beta.tgt, R, maxR, kappa, norm,
                          n, estimate, solver, progress) {
  # ---------------- #
  # Step 1: Initialization
  # ---------------- #
  R.succ <- 0
  R.eval <- 0
  eval.count <- 0
  ub.list <- list()
  lb.list <- list()
  error.list <- list()
  kappa.error <- list()

  # Check if there is any list objects in 'lpmodel'
  any.list <- lpmodel.anylist(lpmodel)

  # If there is some list objects, set maxR as the max length of the list
  if (isTRUE(any.list$list)) {
    maxR <- length(any.list$consol)
  }

  # ---------------- #
  # Step 2: Bootstrap replications
  # ---------------- #
  while ((R.succ < R) & (R.eval != maxR)) {
    # Evaluate the list of indices to be passed to 'future_lapply'
    bs.temp <- bs.assign(R, R.eval, R.succ, maxR, any.list, lpmodel, data,
                         n)
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
      chorussell.return <- future.apply::future_lapply(bs.list,
                                                       FUN = chorussell.bs.fn,
                                                       future.seed = TRUE,
                                                       data = data,
                                                       lpmodel = lpmodel,
                                                       beta.tgt = beta.tgt,
                                                       kappa = kappa,
                                                       norm = norm,
                                                       n = n,
                                                       estimate = estimate,
                                                       solver = solver,
                                                       pbar = pbar,
                                                       eval.count = eval.count,
                                                       n.bs = i1 - i0 + 1,
                                                       progress = progress)
      eval.count <- eval.count + 1
    })

    # Update the list and parameters
    post.return <- post.bs(chorussell.return, i0, i1, R.eval, T.list = NULL,
                           beta.list = NULL, error.list, list.param = kappa,
                           error.param = kappa.error, ub.list = ub.list,
                           lb.list = lb.list)

    ub.list <- post.return$ub.list
    lb.list <- post.return$lb.list
    kappa.error <- post.return$error.param
    error.list <- post.return$error.list
    R.succ <- length(lb.list)
    R.eval <- post.return$R.eval
  }

  # ---------------- #
  # Step 3: Consolidate the upper and lower bounds
  # ---------------- #
  ub.bs <- list()
  lb.bs <- list()

  # Unlist ub.bs and lb.bs
  ubb <- unlist(ub.list, use.names = FALSE)
  lbb <- unlist(lb.list, use.names = FALSE)

  for (i in seq_along(unlist(kappa))) {
    # Make sure that the remainder is 0 if i == length(kappa)
    if (i == length(kappa)) {
      k <- 0
    } else {
      k <- i
    }
    ub.bs[[i]] <- unlist(ubb[seq_along(ubb) %% length(kappa) == k],
                         use.names = FALSE)
    lb.bs[[i]] <- unlist(lbb[seq_along(lbb) %% length(kappa) == k],
                         use.names = FALSE)
  }

  # ---------------- #
  # Step 4: Consolidate the error messages
  # ---------------- #
  if (R.eval != R.succ) {
    # Create data.frame for error messages
    df.error <- data.frame(id = NA,
                           kappa = unlist(kappa.error),
                           message = unlist(error.list))

    # Match the id of the error messages
    df.error <- error.id.match(error.list, df.error)
  } else {
    df.error <- NULL
  }

  return(list(ub = ub.bs,
              lb = lb.bs,
              df.error = df.error,
              error.list = error.list,
              R.succ = R.succ,
              R.eval = R.eval))
}

#' Carries out one bootstrap replication for the Cho-Russell procedure
#'
#' @description This function carries out one bootstrap replication of the
#'   Cho-Russell procedure. This function is used in the
#'   \code{\link[lpinfer]{chorussell.bs}} function via the \code{future_lapply}
#'   command.
#'
#' @inheritParams chorussell
#' @inheritParams chorussell.bs
#' @inheritParams dkqs.bs
#' @inheritParams dkqs.bs.fn
#' @param x This is either the list of indices that represent the bootstrap
#'   replications, or the list of bootstrap components of the \code{lpmodel}
#'   object passed from the user.
#'
#' @return Returns a list of output that are obtained from the Cho-Russell
#'   procedure:
#'   \item{ub}{The bootstrap estimate of the upper bound.}
#'   \item{ub}{The bootstrap estimate of the lower bound.}
#'   \item{kappa.error}{The \code{kappa} parameter that leads to an error (if
#'     applicable).}
#'   \item{msg}{An error message (if applicable).}
#'
#' @export
#'
chorussell.bs.fn <- function(x, data, lpmodel, beta.tgt, kappa, norm, n,
                             estimate, solver, pbar, eval.count, n.bs,
                             progress) {
  # ---------------- #
  # Step 1: Print the progress bar
  # ---------------- #
  if (isTRUE(progress)) {
    if (eval.count == 0) {
      pbar(sprintf("(Computing %s bootstrap estimates)", n.bs))
    } else {
      pbar(sprintf("(Computing %s extra bootstrap estimates)", n.bs))
    }
  }

  # ---------------- #
  # Step 2: Initialization
  # ---------------- #
  # Initialization
  ub.list <- list()
  lb.list <- list()
  kappa.error <- NULL
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
  } else {
    data.bs <- NULL
  }

  for (i in seq_along(kappa)) {
    # Bootstrap estimator
    result <- tryCatch({
      estb.return <- estbounds(data.bs, lpm, kappa[i], norm, estimate, solver)
      estb.return$status <- "NOERROR"
      estb.return
    }, warning = function(w) {
      return(list(status = "warning",
                  msg = w))
    }, error = function(e) {
      return(list(status = "error",
                  msg = e))
    })

    if (!(result$status %in% c("warning", "error"))) {
      ub.list[[i]] <- result$ub
      lb.list[[i]] <- result$lb
      kappa.error <- NULL
      msg <- NULL
    } else {
      ub.list <- NULL
      lb.list <- NULL
      kappa.error <- kappa[i]
      msg <- result$msg$message
      break()
    }
  }

  return(list(ub = ub.list,
              lb = lb.list,
              kappa.error = kappa.error,
              msg = msg))
}

#' Computes the required object in the \code{\link[lpinfer]{chorussell}}
#' procedure
#'
#' @description This function computes the required object in the
#'   \code{\link[lpinfer]{chorussell}} procedure. If \code{ci} is \code{TRUE},
#'   this function returns the \eqn{(1-\alpha)}-confidence interval. Otherwise,
#'   it returns the \eqn{p}-value.
#'
#' @inheritParams chorussell
#' @inheritParams chorussell.simp
#' @param logical.ub The logical upper bound.
#' @param logical.lb The logical lower bound.
#' @param ub The sample upper bound.
#' @param lb The sample lower bound.
#'
#' @return Depending on \code{ci}, this function either returns the
#'   \eqn{p}-value or the \eqn{(1-\alpha)}-confidence interval.
#'
#' @export
#'
chorussell.eval <- function(beta.tgt, lb.can1, lb.can2, ub.can1, ub.can2, n, R,
                            ci, alpha, tol, ub, lb, logical.ub, logical.lb,
                            remove.const, kappa, progress, df.lb = NULL,
                            df.ub = NULL) {
  if (isFALSE(ci)) {
    # ---------------- #
    # Case 1: ci == FALSE, i.e. computes the p-value via bisection method
    # ---------------- #
    # Predefine the two end-points
    a <- 0
    b <- 1

    # Initialization
    k <- 2
    while (abs(b - a) > tol) {
      c <- (a + b)/2
      c.lp <- chorussell.lp(lb.can1, lb.can2, ub.can1, ub.can2, n, R, c, ub,
                            lb, logical.ub, logical.lb, remove.const, ci,
                            kappa, k, progress, df.lb, df.ub)
      c.inout <- chorussell.pt(c.lp, beta.tgt)
      if (isFALSE(c.inout)) {
        b <- c
      } else {
        a <- c
      }
      k <- k + 1
    }
    return(list(pval = c))
  } else {
    # ---------------- #
    # Case 2: ci == TRUE, i.e. computes the (1 - alpha) confidence interval
    # ---------------- #
    # Computes the confidence interval for each alpha
    cr.bd.return <- list()
    for (j in seq_along(alpha)) {
      cr.bd.temp <- chorussell.lp(lb.can1, lb.can2, ub.can1, ub.can2,  n, R,
                                  alpha[j], ub, lb, logical.ub, logical.lb,
                                  remove.const, ci, kappa, 0, progress, df.lb,
                                  df.ub)
      cr.bd.return[[j]] <- cr.bd.temp
    }
    return(cr.bd.return)
  }
}

#' Checks if \code{beta.tgt} is inside the \eqn{(1-\alpha)}-confidence interval
#'
#' @description This function checks whether \code{beta.tgt} is inside the
#'   \eqn{(1-\alpha)}-confidence interval from the
#'   \code{\link[lpinfer]{chorussell}} procedure.
#'
#' @param cr.lp.return List of objects returned from
#'   \code{\link[lpinfer]{chorussell}}.
#' @inheritParams chorussell
#'
#' @return Returns the decision.
#'   \item{decision}{A boolean variable that equals to \code{TRUE} if
#'     \code{beta.tgt} is inside the \eqn{(1-\alpha)}-confidence interval, and
#'     equals to \code{FALSE} otherwise.}
#'
#' @export
#'
chorussell.pt <- function(cr.lp.return, beta.tgt) {
  # Obtain the lower and upper bounds
  lb <- cr.lp.return$bd[1]
  ub <- cr.lp.return$bd[2]

  # Obtain the decision
  if ((beta.tgt <= ub) & (beta.tgt >= lb)) {
    decision <- TRUE
  } else {
    decision <- FALSE
  }

  return(decision)
}


#' Computes the \eqn{(1-\alpha)}-confidence interval in the
#' \code{\link[lpinfer]{chorussell}} procedure
#'
#' @description This function computes the \eqn{(1-\alpha)}-confidence
#' interval in the \code{\link[lpinfer]{chorussell}} procedure by solving the
#' optimization problem.
#'
#' @import future.apply progressr
#'
#' @inheritParams chorussell
#' @inheritParams chorussell.simp
#' @inheritParams chorussell.eval
#' @param k Iteration number.
#'
#' @return Returns the following list of objects:
#'   \item{bd}{A vector that represents the \eqn{(1-\alpha)}-confidence
#'     interval.}
#'   \item{unique}{An indicator variable of whether the solution is unique.}
#'
#' @export
#'
chorussell.lp <- function(lb.can1, lb.can2, ub.can1, ub.can2, n, R, alpha,
                          ub, lb, logical.ub, logical.lb, remove.const, ci,
                          kappa, k = 0, progress, df.lb, df.ub) {
  # ---------------- #
  # Step 1: Select the candidates of the optimization problem
  # ---------------- #
  if (isTRUE(remove.const)) {
    lb.can <- df.lb[df.lb$sum >= 1 - alpha, "bd"]
    ub.can <- df.ub[df.ub$sum >= 1 - alpha, "bd"]
  } else {
    lb.can <- c(lb.can1, lb.can2)
    ub.can <- -c(ub.can1, ub.can2)
  }

  # Include the logical bounds if none of the candidates satisfy the answer
  delta <- ub - lb
  if (is.null(lb.can)) {
    lb.can <- c(sqrt(n) * (logical.lb - lb),
                sqrt(n) * (logical.lb - lb - delta))
  }

  if (is.null(ub.can)) {
    ub.can <- -c(sqrt(n) * (logical.ub - ub),
                 sqrt(n) * (logical.ub - ub + delta))
  }

  # ---------------- #
  # Step 2: Compute the list of bounds that satisfy the constraints
  # of the optimization problem
  # ---------------- #
  if ((lb == -Inf) | (ub == Inf)) {
    # This corresponds to the case where one of the bounds is unbounded based
    # on sample data
    progressr::with_progress({
      if (isTRUE(progress)) {
        pbar <- progressr::progressor(along = seq_along(lb.can))
      } else {
        pbar <- NULL
      }
      cr.lp.return <- future.apply::future_lapply(lb.can,
                                                  FUN = chorussell.lp.fn.unbd,
                                                  lb.can1 = lb.can1,
                                                  lb.can2 = lb.can2,
                                                  ub.can1 = ub.can1,
                                                  ub.can2 = ub.can2,
                                                  ub.can = ub.can,
                                                  lb = lb,
                                                  ub = ub,
                                                  alpha = alpha,
                                                  pbar = pbar,
                                                  ci = ci,
                                                  kappa = kappa,
                                                  k = k,
                                                  progress = progress)
    })
  } else {
    # This corresponds to the case where both of the bounds are bounded based
    # on sample data
    progressr::with_progress({
      if (isTRUE(progress)) {
        pbar <- progressr::progressor(along = seq_along(lb.can))
      } else {
        pbar <- NULL
      }
      cr.lp.return <- future.apply::future_lapply(lb.can,
                                                  FUN = chorussell.lp.fn,
                                                  lb.can1 = lb.can1,
                                                  lb.can2 = lb.can2,
                                                  ub.can1 = ub.can1,
                                                  ub.can2 = ub.can2,
                                                  ub.can = ub.can,
                                                  alpha = alpha,
                                                  pbar = pbar,
                                                  ci = ci,
                                                  kappa = kappa,
                                                  k = k,
                                                  progress = progress)
    })
  }

  # Merge the list of tables from 'cr.lp.return'
  cr.lp.bds <- Reduce(rbind, cr.lp.return)

  # ---------------- #
  # Step 3: Compute the tightest bounds
  # ---------------- #
  # Find the smallest length
  min.len <- min(cr.lp.bds$len)

  # Find the corresponding bounds
  min.bd <- unique(cr.lp.bds[cr.lp.bds$len == min.len, 1:2])

  # Check if its unique
  if (nrow(min.bd) == 1) {
    unique.bd <- 1
  } else {
    unique.bd <- 0
  }

  # Compute the bounds
  bd <- c(lb - min.bd[1,1]/sqrt(n), ub + min.bd[1,2]/sqrt(n))

  return(list(bd = bd,
              unique = unique.bd))
}

#' Computes whether the candidate bounds satisfy the constraints in the
#' \code{\link[lpinfer]{chorussell}} procedure when the sample estimates are
#' non-infinite
#'
#' @description This function carries out one set of computation in the
#'   \code{\link[lpinfer]{chorussell}} procedure to check whether the set of
#'   values satisfy the constraints in the \code{\link[lpinfer]{chorussell}}
#'   procedure. In particular, given one candidate of the lower bound, it
#'   checks all possible pairs with the upper bounds and determine whether they
#'   satisfy the constraints.
#'
#' @param x A candidate lower bound.
#' @param ub.can All candidates of the upper bounds.
#' @inheritParams chorussell
#' @inheritParams chorussell.simp
#' @inheritParams chorussell.simp.fn
#' @inheritParams dkqs.bs.fn
#' @inheritParams chorussell.lp
#'
#' @return Returns a data frame.
#'   \item{df}{A data frame that contains the pairs of feasible lower bounds
#'     and upper bounds, as well as their sum.}
#'
#' @export
#'
chorussell.lp.fn <- function(x, lb.can1, lb.can2, ub.can1, ub.can2, ub.can,
                             alpha, pbar, ci, kappa, k, progress) {
  # ---------------- #
  # Step 1: Print the progress bar
  # ---------------- #
  if (isTRUE(progress)) {
    if (isFALSE(ci)) {
      pbar(sprintf("(Computing p-value for kappa = %s, iteration %s)",
                   kappa, k))
    } else {
      pbar(sprintf("(Constructing confidence interval for kappa = %s)", kappa))
    }
  }

  # ---------------- #
  # Step 2: Initialize data.frame
  # ---------------- #
  df <- data.frame(matrix(vector(), nrow = 0, ncol = 3))
  colnames(df) <- c("lb", "ub", "len")

  # ---------------- #
  # Step 3: Check if the candidate bounds satisfy the inequalities
  # ---------------- #
  for (i in seq_along(ub.can)) {
    ind1 <- (mean((lb.can1 <= x) * (-ub.can[i] <= ub.can2)) >= (1 - alpha))
    ind2 <- (mean((lb.can2 <= x) * (-ub.can[i] <= ub.can1)) >= (1 - alpha))
    if (isTRUE(ind1) & isTRUE(ind2)) {
      j <- nrow(df) + 1
      df[j, 1] <- x
      df[j, 2] <- ub.can[i]
      df[j, 3] <- ub.can[i] + x
    }
  }

  return(df)
}

#' Computes whether the candidate bounds satisfy the constraints in the
#' \code{\link[lpinfer]{chorussell}} procedure when one of the sample
#' estimates of the bounds is infinite
#'
#' @description This function finds the candidate bounds that satisfy
#'   the constraints in the \code{\link[lpinfer]{chorussell}} procedure when
#'   one of the sample estimates of the bounds is infinite. If the upper bound
#'   is infinite, then the third column of the data frame will store the
#'   negative of the lower bound so the largest candidate for the lower bound
#'   will be chosen in order to minimize the interval length. Similarly, if the
#'   lower bound is infinite, then the third column of the data frame will
#'   store the upper bound.
#'
#' @inheritParams chorussell
#' @inheritParams chorussell.simp
#' @inheritParams chorussell.simp.fn
#' @inheritParams chorussell.lp.fn
#'
#' @return Returns a data frame.
#'   \item{df}{A data frame that contains the pairs of feasible lower bounds
#'     and upper bounds, as well as their sum.}
#'
#' @export
#'
chorussell.lp.fn.unbd <- function(x, lb.can1, lb.can2, ub.can1, ub.can2,
                                  lb.can, ub.can, lb, ub, alpha, pbar, ci, k,
                                  progress) {
  # ---------------- #
  # Step 1: Print the progress bar
  # ---------------- #
  if (isTRUE(progress)) {
    if (isFALSE(ci)) {
      pbar(sprintf("(Computing p-value for kappa = %s, iteration %s)",
                   kappa, k))
    } else {
      pbar(sprintf("(Constructing confidence interval for kappa = %s)", kappa))
    }
  }

  # ---------------- #
  # Step 2: Initialize data.frame
  # ---------------- #
  df <- data.frame(matrix(vector(), nrow = 0, ncol = 3))
  colnames(df) <- c("lb", "ub", "len")

  # ---------------- #
  # Step 3: Assign the bounds depending on whether lb or ub is bounded
  # ---------------- #
  if (lb == -Inf) {
    can1 <- -ub.can1
    can2 <- -ub.can2
    x <- -x
  } else {
    can1 <- lb.can1
    can2 <- lb.can2
  }

  # ---------------- #
  # Step 4: Check if the candidate bounds satisfy the inequalities
  # ---------------- #
  for (i in seq_along(ub.can)) {
    ind1 <- (mean((can1 <= x)) >= (1 - alpha))
    ind2 <- (mean((can2 <= x)) >= (1 - alpha))
    if (isTRUE(ind1) & isTRUE(ind2)) {
      j <- nrow(df) + 1
      if (lb == -Inf) {
        df[j, 1] <- -Inf
        df[j, 2] <- x
      } else {
        df[j, 1] <- x
        df[j, 2] <- Inf
      }
      df[j, 3] <- x
    }
  }

  return(df)
}

#' Checks and updates the input in the \code{\link[lpinfer]{chorussell}}
#' procedure
#'
#' @description This function checks and updates the input from the user in the
#'    \code{\link[lpinfer]{chorussell}} function. If there is any invalid input,
#'    the function will be terminated and error messages will be printed.
#'
#' @inheritParams chorussell
#' @inheritParams dkqs
#'
#' @return Returns the updated parameters back to the function
#'   \code{chorussell}. The following information are updated:
#'    \itemize{
#'       \item{\code{ci}}
#'       \item{\code{data}}
#'       \item{\code{solver}}
#'       \item{\code{solver.name}}
#'       \item{\code{test.logical}}
#'       \item{\code{n}}
#'       \item{\code{logical.lb}}
#'       \item{\code{logical.ub}}
#'    }
#'
#' @export
#'
chorussell.check <- function(data, lpmodel, beta.tgt, R, Rmulti, kappa,
                             norm, n, estimate, solver, ci, alpha, tol,
                             progress) {
  # ---------------- #
  # Step 1: Check data
  # ---------------- #
  # Check data. If data is NULL, check if n is a positive integer
  if (!is.null(data)) {
    data <- check.dataframe(data)
    n <- nrow(data)
  } else {
    check.samplesize(n, "n")
  }

  # ---------------- #
  # Step 2: Check lpmodel
  # ---------------- #
  lpmodel <- check.lpmodel(data = data,
                           lpmodel = lpmodel,
                           name.var = "lpmodel",
                           A.tgt.cat = c("matrix", "function_mat", "list"),
                           A.obs.cat = c("matrix", "function_mat", "list"),
                           A.shp.cat = c("matrix", "function_mat", "list"),
                           beta.obs.cat = c("function_mat",
                                            "list",
                                            "function_obs_var"),
                           beta.shp.cat = c("matrix", "function_mat", "list"),
                           R = R)

  # ---------------- #
  # Step 3: Check solver
  # ---------------- #
  solver.return <- check.solver(solver, "solver")
  solver <- solver.return$solver
  solver.name <- solver.return$solver.name

  # ---------------- #
  # Step 4: Assign ci
  # ---------------- #
  if (is.null(beta.tgt) & is.null(ci)) {
    # If beta.tgt and ci are not passed, build a confidence interval
    ci <- TRUE
  } else if ((!is.null(beta.tgt)) & is.null(ci)) {
    # If beta.tgt is given but ci is not passed, compute the p-value
    ci <- FALSE
  }

  # ---------------- #
  # Step 5: Check numerics
  # ---------------- #
  if (isFALSE(ci)) {
    check.numeric(beta.tgt, "beta.tgt") 
  }
  check.positiveinteger(R, "R")
  check.nonnegative(tol, "tol")
  check.norm(norm, "norm")
  check.numrange(Rmulti, "Rmulti", "closed", 1, "open", Inf)

  # Alpha can be a vector
  for (i in seq_along(alpha)) {
    check.numrange(alpha[i], "alpha", "closed", 0, "closed", 1)
  }

  # Kappa can be a vector
  for (i in seq_along(kappa)) {
    check.nonnegative(kappa[i], "kappa")
  }

  # ---------------- #
  # Step 6: Check whether beta.tgt is within the logical bounds
  # ---------------- #
  if (isFALSE(ci)) {
    test.return <- check.betatgt(data, lpmodel, beta.tgt, solver)
    test.logical <- test.return$inout
    logical.lb <- test.return$lb
    logical.ub <- test.return$ub 
  } else {
    test.logical <- 1
    logical.lb <- check.betatgt.lp(data, lpmodel, "min", solver)
    logical.ub <- check.betatgt.lp(data, lpmodel, "max", solver)
  }

  # ---------------- #
  # Step 7: Check Boolean
  # ---------------- #
  check.boolean(estimate, "estimate")
  check.boolean(ci, "ci")
  check.boolean(progress, "progress")

  # ---------------- #
  # Step 8: Return updated items
  # ---------------- #
  return(list(ci = ci,
              data = data,
              solver = solver,
              solver.name = solver.name,
              test.logical = test.logical,
              n = n,
              logical.lb = logical.lb,
              logical.ub = logical.ub))
}

#' Print results from \code{\link[lpinfer]{chorussell}}
#'
#' @description This function either prints the \eqn{p}-values or a
#'   \eqn{(1-\alpha)} confidence interval from
#'   \code{\link[lpinfer]{chorussell}}.
#'
#' @param x The output objects returned from \code{\link[lpinfer]{chorussell}}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.chorussell <- function(x, ...) {
  # Prints confidence interval if ci is TRUE. Otherwise, prints p-value
  digits <- 5
  if ((x$test.logical == 1) | isTRUE(x$ci)) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Prints confidence interval if ci is TRUE. Otherwise, prints p-value
    digits <- 5
    if (isFALSE(x$ci)) {
      # Print the p-values
      df.pval <- x$pval

      # Print the p-values
      if (nrow(df.pval) == 1) {
        cat(sprintf("p-value: %s\n", round(df.pval[1, 2], digits = digits)))
      } else {
        cat("p-values: \n")
        colnames(df.pval) <- c("kappa", "p-value")
        df.pval[, 2] <- round(df.pval[, 2], digits = digits)
        print(df.pval, row.names = FALSE)
      }
    } else {
      # Print the confidence intervals
      ci.df <- x$ci.df
      n.ci.df <- nrow(ci.df)

      if (n.ci.df == 1) {
        cat(sprintf("%s%% confidence interval: [%s, %s]\n",
                    round(100 * (1 - x$alpha), digits = digits),
                    round(x$ci.df[1, 3], digits = digits),
                    round(x$ci.df[1, 4], digits = digits)))
      } else {
        cat("Confidence intervals: \n")
        ci.df.int <- data.frame(matrix(vector(), nrow = n.ci.df, ncol = 3))
        colnames(ci.df.int) <- c("Significance level",
                                 "kappa",
                                 "Confidence intervals")
        ci.df.int[, 1] <- ci.df[, 1]
        ci.df.int[, 2] <- ci.df[, 2]
        for (i in 1:n.ci.df) {
          ci.df.int[i, 3] <- sprintf("[%s, %s]",
                                     round(x$ci.df[i, 3], digits = digits),
                                     round(x$ci.df[i, 4], digits = digits))
        }
        print(ci.df.int, row.names = FALSE)
      }
    }
  } else {
    # Case 2: 'beta.tgt' is outside the logical bound
    infeasible.pval.msg()
  }
}

#' Summary of results from \code{\link[lpinfer]{chorussell}}
#'
#' @description This function prints a summary of the results obtained from
#'   \code{\link[lpinfer]{chorussell}}.
#'
#' @inheritParams print.chorussell
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.chorussell <- function(x, ...) {
  if ((x$test.logical == 1) | isTRUE(x$ci)) {
    # Case 1: 'beta.tgt' is within the logical bound
    # Print the p-value or the corresponding confidence interval
    print(x)

    # Print kappa if x$pval or x$ci.df is one dimensional
    if ((isTRUE(x$ci) & isTRUE(nrow(x$ci.df) == 1)) |
        (isFALSE(x$ci) & isTRUE(nrow(x$pval) == 1))) {
      cat(sprintf("kappa: %s\n", x$kappa))
    }
    
    # Print test statistic, solver used, norm used, and the number of
    # successful bootstrap replications
    cat(sprintf("Norm: %s\n", x$norm))
    cat(sprintf("Solver: %s\n", x$solver))
    cat(sprintf("Number of successful bootstrap replications: %s\n", x$R.succ))

    # Number of failed bootstrap replications
    if (!is.null(x$df.error)) {
      nerr <- nrow(x$df.error)
      errstring <- "Number of failed bootstrap"
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
