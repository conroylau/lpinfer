#' Constructs confidence interval
#'
#' @description This function constructs the confidence interval using the
#'    bisection method.
#'
#' @param f The function that represents a testing procedure.
#' @param farg The list of arguments to be passed to the function of testing
#'    procedure.
#' @param alpha The significance level(s). This can be a vector.
#' @param init.lb The initial brackets to search for the lower bound. This is
#'   not required if the \code{\link[lpinfer]{chorussell}} is used.
#' @param init.ub The initial brackets to search for the upper bound. This is
#'   not required if the \code{\link[lpinfer]{chorussell}} is used.
#' @param tol The tolerance level in the bisection method.
#' @param max.iter The maximum number of iterations in the bisection method.
#' @param pvals The data frame that consists the points and the corresponding
#'    \eqn{p}-values that have been tested in the previous iterations.
#' @param dp The number of decimal places in the output.
#' @param progress The boolean variable for whether the result messages should
#'    be displayed in the procedure of constructing confidence interval. If
#'    it is set as \code{TRUE}, the messages are displayed throughout the
#'    procedure. Otherwise, the messages will not be displayed.
#'
#' @returns Returns the confidence interval and a data frame that contains the
#'    points being tested in the procedure.
#'    \item{pvals}{The data frame that consists of the points and the
#'       corresponding \eqn{p}-values that have been tested in constructing
#'       the confidence intervals.}
#'    \item{df_ub}{The data frame storing the information for the bisection
#'       method in each iteration when evaluating the upper bound of the
#'       confidence interval.}
#'    \item{df_lb}{The data frame storing the information for the bisection
#'       method in each iteration when evaluating the lower bound of the
#'       confidence interval.}
#'    \item{alpha}{The significance levels.}
#'    \item{tol}{The tolerance level in the bisection method.}
#'    \item{iter}{The total number of iterations taken.}
#'    \item{call}{The matched call.}
#'    \item{para.name}{The name of the tuning parameters involved.}
#'    \item{para.vals}{The values of the tuning parameters involved.}
#'    \item{ci}{The confidence intervals constructed.}
#'
#' @details The number of decimal places displayed in the messages (if
#'    \code{progress} is set as \code{TRUE}) is equal to the number of decimal
#'    places in the variable \code{tol}.
#'
#' @example ./inst/example/invertci_example.R
#'
#' @export
#'
invertci <- function(f, farg = list(), alpha = .05, init.lb = NULL,
                     init.ub = NULL, tol = .0001, max.iter = 20,
                     pvals = NULL, dp = 5, progress = TRUE) {
  # ---------------- #
  # Step 1: Update call, check and update the arguments
  # ---------------- #
  # Extract the current RNG state
  rngstate <- .Random.seed

  # Obtain call information
  call <- match.call()

  # Check and update
  invertci.return <- invertci.check(f, farg, alpha, init.lb, init.ub, tol,
                                    max.iter, pvals, progress)

  # Updates the input
  lb0 <- invertci.return$lb0
  lb1 <- invertci.return$lb1
  ub0 <- invertci.return$ub0
  ub1 <- invertci.return$ub1

  # Return pvals
  pvals <- invertci.return$pvals

  # Initialize lists
  df_ub_list <- NULL
  df_lb_list <- NULL
  termination_list <- NULL
  ub_list <- NULL
  lb_list <- NULL
  iter_list <- NULL

  # Initialize message for constructing confidence intervals
  comp.bound <- "\n === Computing %s bound of confidence interval ===\n"

  # ---------------- #
  # Step 2: Compute the left end-point for the first tuning parameter and
  # the first alpha
  # ---------------- #
  # This step is skipped if the chorussell procedure is used
  if (!identical(f, chorussell)) {
    # Append the test point
    farg$beta.tgt <- ub1

    # Run the procedure once and obtain the pval data frame
    init.return <- do.call(f, farg)
    pval <- init.return$pval

    # Parameter names are located in the first (n - 1) columns of pval
    pval.col <- ncol(pval)

    # Parameter name and the list of parameters
    para.name <- colnames(pval)[1:(pval.col - 1)]
    para.vals <- data.frame(pval[, 1:(pval.col - 1)])
    colnames(para.vals) <- para.name
  }

  if (is.null(pvals)) {
    ## Construct pvals if it is null
    pvals <- data.frame(matrix(vector(),
                               nrow = 0,
                               ncol = length(para.name) + 2))
    colnames(pvals) <- c(para.name, "point", "p-value")
  } else {
    # Check pvals if it already exists
    ## Check whether the number of columns match
    if (ncol(pvals) != length(para.name) + 2) {
      stop(paste0("The number of columns in 'pvals' needs to equal to ",
                  "the number of tuning parameters that can be ",
                  "multivalued plus 2."))
    }

    ## Check whether the column names match
    if (!setequal(c(para.name, "point", "p-value"), colnames(pvals))) {
      stop(paste0("The column names in 'pvals' need to contain the names ",
                  "of the tuning parameters that can be multivalued and ",
                  "the two strings 'point' and 'p-value'."))
    }
  }

  # ---------------- #
  # Step 3: Return confidence interval and data frame
  # ---------------- #
  # Sort alpha
  alpha <- sort(alpha)

  # Initialize empty lists
  df_ub_list <- list()
  df_lb_list <- list()
  termination_list <- list()
  lb_list <- list()
  ub_list <- list()
  iter_list <- list()

  # Loop through each alpha and each set of tuning parameters
  for (i in seq_along(alpha)) {
    # Initialize list per alpha
    df_ub_list[[i]] <- list()
    df_lb_list[[i]] <- list()
    termination_list[[i]] <- list()
    lb_list[[i]] <- list()
    ub_list[[i]] <- list()
    iter_list[[i]] <- list()

    if (!identical(f, chorussell)) {
      for (j in 1:nrow(pval)) {
        # Assign the new farg object based on the j-th set of tuning parameters
        for (k in seq_along(para.name)) {
          farg[[para.name[k]]] <- para.vals[j, k]
        }

        # If `f` is not `chorussell`, use the bisection method in `invertci`
        if (i > 1 & (isTRUE(progress))) {
          cat("\n")
        }

        # Append the result that has already been evaluated in step 2
        pvals[nrow(pvals) + 1, "point"] <- ub1
        pvals[nrow(pvals), "p-value"] <- pval[j, 2]
        pvals[nrow(pvals), 1:(ncol(pvals) - 2)] <- para.vals[j, ]

        para.match <- para.vals[j, ]
        if (is.null(nrow(para.match))) {
          para.match <- data.frame(para.match)
          colnames(para.match) <- para.name
        }

        ### Compute upper bound of confidence interval
        termination <- NULL
        if (isTRUE(progress)) {
          # Print the significance level being considered
          cat(sprintf(paste0("********** Constructing %s%% ",
                             "confidence interval **********"),
                      round((1 - alpha[i]) * 100, digits = 5)))

          # Print the parameters being considered
          invertci.show.param(para.name, para.vals, j)

          # Print the interval being considered
          cat(sprintf(comp.bound, "upper"))
        }
        ub_return <- ci.bisection(f, farg, alpha[i], ub1, ub0, tol, max.iter,
                                  pvals, progress, 1, dp, rngstate, para.match)
        # Update data frame
        pvals <- ub_return$pvals
        # Data frame storing all messages in each iteration
        df_ub <- ub_return$df_bis
        # Obtain termination message
        termination$ub <- ub_return$last_iter_msg
        # Obtain upper bound
        ub <- ub_return$pt

        ### Compute lower bound of confidence interval
        if (isTRUE(progress)) {
          cat(sprintf(comp.bound, "lower"))
        }
        lb_return <- ci.bisection(f, farg, alpha[i], lb0, lb1, tol, max.iter,
                                  pvals, progress, -1, dp, rngstate, para.match)
        # Update data frame
        pvals <- lb_return$pvals
        # Data frame storing all messages in each iteration
        df_lb <- lb_return$df_bis
        # Obtain termination message
        termination$lb <- lb_return$last_iter_msg
        # Obtain lower bound
        lb <- lb_return$pt

        # Compute the number of Iterations
        iter <- lb_return$iter + ub_return$iter

        # ---------------- #
        # Step 4: Store the information as a list
        # ---------------- #
        df_ub_list[[i]][[j]] <- df_ub
        df_lb_list[[i]][[j]] <- df_lb
        termination_list[[i]][[j]] <- termination
        lb_list[[i]][[j]] <- lb
        ub_list[[i]][[j]] <- ub
        iter_list[[i]][[j]] <- iter
      }
    } else {
      # If `f` is `chorussell`, use the bisection method in `chorussell`
      farg$tol <- tol
      farg$alpha <- alpha[i]
      farg$ci <- TRUE
      chorussell.return <- do.call(chorussell, farg)
      lb <- chorussell.return$ci.df[, 3]
      ub <- chorussell.return$ci.df[, 4]
      iter <- chorussell.return$iter

      # Consolidate the results
      lb_list[[i]] <- lb
      ub_list[[i]] <- ub
      iter_list[[i]] <- iter
      para.vals <- data.frame(chorussell.return$ci.df$kappa)
      para.name <- "kappa"
    }
  }

  # ---------------- #
  # Step 5: Consolidate the data frame of confidence intervals
  # ---------------- #
  # Initialize data frame
  nr.p <- nrow(para.vals)
  nc.p <- ncol(para.vals)
  ci <- data.frame(matrix(vector(),
                          nrow = nr.p * length(alpha),
                          ncol = ncol(para.vals) + 3))
  colnames(ci) <- c("alpha", para.name, "lb", "ub")

  # Consolidate data
  for (i in seq_along(alpha)) {
    ## Append the alphas
    ci[(nr.p * (i - 1) + 1):(nr.p * i), 1] <- alpha[i]

    ## Append the tuning parameters
    ci[(nr.p * (i - 1) + 1):(nr.p * i), 2:(nc.p + 1)] <- para.vals

    ## Append the confidence intervals
    for (j in 1:nr.p) {
      ci[j + (i - 1) * nr.p, ncol(ci) - 1] <- lb_list[[i]][[j]]
      ci[j + (i - 1) * nr.p, ncol(ci)] <- ub_list[[i]][[j]]
    }
  }

  # ---------------- #
  # Step 6: Assign the return list and return output
  # ---------------- #
  output <- list(pvals = pvals,
                 tol = tol,
                 alpha = alpha,
                 iter = iter_list,
                 call = call,
                 ci = ci)

  # The following information are returned only if `f` is not `chorussell`
  if (!identical(f, chorussell)) {
    output$max.iter <- max.iter
    output$df_ub <- df_ub_list
    output$df_lb <- df_lb_list
    output$termination <- termination_list
  }
  output$para.name <- para.name
  output$para.vals <- para.vals
  attr(output, "class") <- "invertci"

  # Return output
  return(output)
}

#' Bisection method for constructing confidence intervals
#'
#' @description This function constructs the two-sided confidence interval
#'    of a given testing procedure using the bisection method.
#'
#' @param type The type of the confidence interval. Set \code{type} as 1 for the
#'    upper bound of the confidence interval and set \code{type} as -1 for the
#'    lower bound of the confidence interval.
#' @param b0 The lower bound of the initial bracket.
#' @param b1 The upper bound of the initial bracket.
#' @param dp The number of decimal places to be displayed for the \eqn{p}-values
#'    and confidence intervals in the messages if \code{progress} is set
#'    as \code{TRUE}.
#' @param rngstate The current RNG state obtained from \code{.Random.seed}.
#' @param para.match The list of parameters to be matched.
#' @inheritParams invertci
#'
#' @return Return the solution of the bisection method and the updated
#'    data frame.
#'    \item{pt}{The last test point.}
#'    \item{iter}{The number of iterations.}
#'    \item{pvals}{The data frame that consists of the points and the
#'       corresponding \eqn{p}-values that have been tested in constructing the
#'       confidence intervals.}
#'    \item{df_bis}{The data frame storing the information for the bisection
#'       method in each iteration.}
#'    \item{last_iter_msg}{The message for the last iteration. This refers to
#'       whether why the bisection method stops.}
#'
#' @export
#'
ci.bisection <- function(f, farg, alpha, b0, b1, tol, max.iter, pvals,
                         progress, type, dp, rngstate, para.match) {

  # ---------------- #
  # Step 1: Evaluate the end-points and the mid-point of b0 and b1
  # ---------------- #
  ### Initialize data frame to collect the information in the bisection method
  df_bis <- data.frame(matrix(vector(), 0, 6,
                              dimnames = list(c(),
                                              c("iteration",
                                                "left",
                                                "right",
                                                "point",
                                                "p-value",
                                                "decision"))),
                       stringsAsFactors = F)
  # Divide alpha by 2
  alpha_2sided <- alpha/2

  ### Left end-point a
  a <- b0
  fb0_return <- bisec.eval(f, farg, a, pvals, rngstate, para.match)
  pvals <- fb0_return$pvals
  # Print information
  if (isTRUE(progress)) {
    cat(paste0(" Iteration\t Lower bound \t Upper bound \t Test point \t ",
               "p-value\t Reject?\n"))
  }
  df_bis <- bisec.print("left end", alpha_2sided, fb0_return, a, "NA",
                        progress, dp, df_bis)$df_bis

  ### Right end-point b
  b <- b1
  fb1_return <- bisec.eval(f, farg, b, pvals, rngstate, para.match)
  pvals <- fb1_return$pvals
  # Print information
  df_bis <- bisec.print("right end", alpha_2sided, fb1_return, "NA", b,
                        progress, dp, df_bis)$df_bis

  # If fb1 and fb0 are of the same sign, ask user to choose another interval
  # Compute mid-point and evaluate the corresponding p-value
  c <- (b + a)/2
  fc_return <- bisec.eval(f, farg, c, pvals, rngstate, para.match)
  fc <- fc_return$pval
  pvals <- fc_return$pvals

  # ---------------- #
  # Step 2: Bisection method
  # ---------------- #
  for (i in 1:max.iter) {
    # Bisection method is completed if the difference between the two points
    # is below the tolerance level.
    if (abs(b - a) < tol) {
      if (isTRUE(progress)) {
        tol_msg <- paste0(" >>> Length of interval is below tolerance level. ",
                          "Bisection method is completed.\n", sep = "")
      }
      last_iter_msg <- "Length of interval is below tolerance level"
      if (isTRUE(progress)) {
        cat(tol_msg)
      }
      break
    }

    # Update interval based on whether the left or the right segment is chosen
    # Print information
    df_bis <- bisec.print(i , alpha_2sided, fc_return, a, b, progress,
                          dp, df_bis)$df_bis

    if (ci.inout(fc, alpha_2sided, type) == "left") {
      b <- c
    } else {
      a <- c
    }

    # Evaluate new mid-point
    c <- (a + b)/2

    # Update data frame and p-value
    fc_return <- bisec.eval(f, farg, c, pvals, rngstate, para.match)
    fc <- fc_return$pval
    pvals <- fc_return$pvals
  }

  # Only called when the maximum number of iterations is reached
  if (i == max.iter) {
    last_iter_msg <- "Reached maximum number of iterations"
    iter_msg <- paste(" >>> Reached the maximum number of iterations. ",
                      "Bisection method is completed.\n", sep = "")
    if (isTRUE(progress)) {
      cat(iter_msg)
    }
  }

  # ---------------- #
  # Step 3: Return results
  # ---------------- #
  invisible(list(pt = c,
                 iter = i,
                 pvals = pvals,
                 df_bis = df_bis,
                 last_iter_msg = last_iter_msg))
}

#' Evaluation of test statistic and check if the point has been evaluated
#'
#' @import plyr
#'
#' @description This function is used in the \code{\link[lpinfer]{invertci}}
#' procedure. It checks if the \eqn{p}-value for the point considered has
#' already been evaluated in previous iterations or provided by the user. The
#' function will compute the \eqn{p}-value if it is not provided. Otherwise,
#' it will use the previous data.
#'
#' @inheritParams invertci
#' @inheritParams ci.bisection
#' @param pt The point to be evaluated in the bisection method.
#'
#' @return Returns the \eqn{p}-value of the point considered and an updated
#'    data frame that contains the points and the \eqn{p}-values.
#'    \item{pval}{The \eqn{p}-value that corresponds to the test point.}
#'    \item{pvals}{The updated data frame that contains the \eqn{p}-values that
#'      correspond to each set of tuning parameters and each test point.}
#'
#' @export
#'
bisec.eval <- function(f, farg, pt, pvals, rngstate, para.match) {
  # ---------------- #
  # Step 1: Check if the data point has appeared in previous iterations.
  # ---------------- #
  para.match.temp <- para.match
  para.match.temp$point <- pt
  suppressMessages(df_match <- plyr::match_df(pvals, para.match.temp))
  df_n <- nrow(pvals)

  # ---------------- #
  # Step 2: If the points have been evaluated in pvals, it will be taken
  # directly. Otherwise, it will be computed.
  # ---------------- #
  if ((is.null(df_match) == TRUE) | (dim(df_match)[1] == 0)) {
    farg$beta.tgt <- pt
    assign(x = ".Random.seed", value = rngstate, envir = .GlobalEnv)
    test_return <- do.call(f, farg)
    if (is.data.frame(test_return$pval)) {
      pval <- test_return$pval[1, 2]
    } else {
      pval <- test_return$pval
    }
    pvals[df_n + 1, 1:(ncol(pvals) - 2)] <- para.match[1, ]
    pvals[df_n + 1, "point"] <- pt
    pvals[df_n + 1, "p-value"] <- pval
  } else {
    pval <- unique(df_match$`p-value`)
  }

  # ---------------- #
  # Step 3: Return results
  # ---------------- #
  return(list(pval = pval,
              pvals = pvals))
}

#' Determine whether a point is inside the confidence interval or not
#'
#' @description This function determines whether the bisection method is going
#'    to be updated by choosing the left segment or the right segment of the
#'    interval as the updated interval for the next iteration.
#'
#' @param pval The \eqn{p}-value of the test statistic.
#' @inheritParams invertci
#' @inheritParams ci.bisection
#'
#' @return Returns whether the part of the interval to be selected in the
#'    next iteration of the bisection method.
#'    \item{part}{A string indicating whether the left or right segment of the
#'    interval is chosen.}
#'
#' @export
#'
ci.inout <- function(pval, alpha, type) {
  if (type == 1) {
    # Type == 1: Upper bound
    if (pval < alpha) {
      part <- "left"
    } else {
      part <- "right"
    }
  } else if (type == -1) {
    # Type == -1: Lower bound
    if (pval < alpha) {
      part <- "right"
    } else {
      part <- "left"
    }
  }

  return(part)
}

#' Print messages in bisection procedure and store results
#'
#' @description This function prints the information in
#'   \code{\link[lpinfer]{ci.bisection}} and store the result for each
#'   iteration.
#'
#' @param procedure The variable indicating whether the function is evaluating
#'    the end-points or first mid-point, or is iterating through the bisection
#'    procedure.
#' @param alphahalf Half of significance value that is used to evaluate the
#'    confidence interval.
#' @param returnlist The list of information obtained from running
#'    \code{\link[lpinfer]{bisec.eval}}.
#' @param a The lower bound of the current interval. This is \code{NULL} if
#'   the initial end-points are being evaluated.
#' @param b The upper bound of the current interval. This is \code{NULL} if
#'   the initial end-points are being evaluated.
#' @param df_bis A data frame storing the information from the bisection
#'   method.
#' @inheritParams invertci
#'
#' @return Return the updated data frame that stores the information for the
#'    the iteration.
#'    \item{df_bis}{The data frame storing the information for the bisection
#'       method in each iteration.}
#'
#' @export
#'
bisec.print <- function(procedure, alphahalf, returnlist, a, b, progress, dp,
                        df_bis) {
  # ---------------- #
  # Step 1: Obtain information about the current data frame
  # ---------------- #
  df_bis_row <- nrow(df_bis)
  # Update decision
  if (returnlist$pval < alphahalf) {
    decision <- TRUE
  } else {
    decision <- FALSE
  }

  # ---------------- #
  # Step 2: Print information from the iteration
  # ---------------- #
  if (is.numeric(procedure) == FALSE) {
    # Case A: 'procedure' is not numeric if evaluating the initial 3 points
    if (procedure == "left end") {
      df_bis[df_bis_row + 1, 4] <- a
    } else if (procedure == "right end") {
      df_bis[df_bis_row + 1, 4] <- b
    }
  } else {
    # Case B: 'procedure' is numeric if evaluating the bisection method
    df_bis[df_bis_row + 1, 4] <- (a + b)/2
  }

  # ---------------- #
  # Step 3: Update data frame
  # ---------------- #
  # Update column 1, i.e. whether evaluating end-points or iterations
  if (procedure == "left end") {
    df_bis[df_bis_row + 1, 1] <- "Left end pt."
  } else if (procedure == "right end") {
    df_bis[df_bis_row + 1, 1] <- "Right end pt."
  } else if (is.numeric(procedure) == TRUE) {
    df_bis[df_bis_row + 1, 1] <- procedure
  }
  df_bis[df_bis_row + 1, 2] <- a
  df_bis[df_bis_row + 1, 3] <- b
  df_bis[df_bis_row + 1, 5] <- returnlist$pval
  df_bis[df_bis_row + 1, 6] <- decision

  # ---------------- #
  # Step 4: Print information
  # ---------------- #
  if (isTRUE(progress)) {
    summary.bisection.print(df_bis, df_bis_row + 1)
  }

  # ---------------- #
  # Step 5: Return information
  # ---------------- #
  invisible(list(df_bis = df_bis))
}

#' Print results from \code{\link[lpinfer]{invertci}}
#'
#' @description This function prints the results from
#'   \code{\link[lpinfer]{invertci}}.
#'
#' @param x The output object returned from \code{\link[lpinfer]{invertci}}.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned
#'
#' @export
#'
print.invertci <- function(x, ...) {
  if ((length(x$alpha) == 1) & (nrow(x$para.vals) == 1)) {
    print.invertci_single(x)
  } else {
    print.invertci_multiple(x)
  }
}

#' Print results from \code{\link[lpinfer]{invertci}} with a single
#' significance level
#'
#' @description This function prints the results from
#'   \code{\link[lpinfer]{invertci}}.
#'
#' @inheritParams print.invertci
#'
#' @return Nothing is returned
#'
#' @export
#'
print.invertci_single <- function(x, ...) {
  cat(sprintf("Confidence interval: [%s, %s]\n",
              round(x$ci[1, 3], digits = 5),
              round(x$ci[1, 4], digits = 5)))
}

#' Print results from \code{\link[lpinfer]{invertci}} with multiple
#' significance levels
#'
#' @description This function uses the print method on the return list of the
#'    function \code{\link[lpinfer]{invertci}}.
#'
#' @inheritParams print.invertci
#'
#' @return Nothing is returned
#'
#' @export
#'
print.invertci_multiple <- function(x, alphas = NULL, ...) {
  cat("Confidence intervals:\n")
  df <- x$ci
  # Create column of confidence intervals
  colnames(df)[1] <- "Significance level"
  df$ci <- sprintf("[%s, %s]",
                   round(df$lb, digits = 5),
                   round(df$ub, digits = 5))
  colnames(df)[ncol(df)] <- "Confidence interval"

  # Drop the columns of upper and lower bounds
  df$lb <- NULL
  df$ub <- NULL

  # Print object
  if (is.null(alphas)) {
    print(df, row.names = FALSE, digits = 5)
  } else {
    print(df[df[, 1] %in% alphas, ], row.names = FALSE, digits = 5)
  }
}

#' Summary of results from \code{\link[lpinfer]{invertci}}
#'
#' @description This function summarizes the results for
#'   \code{\link[lpinfer]{invertci}}.
#'
#' @param x The output object returned from \code{\link[lpinfer]{invertci}}.
#' @param alphas A list of alphas that the user would like to print.
#' @inheritParams print.invertci
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.invertci <- function(x, alphas = NULL, ...) {
  # String containing the general message to be printed
  msg.bound <- "=== Iterations in constructing %s bound:"

  if ((length(x$alpha) == 1) & (nrow(x$para.vals) == 1)) {
    summary.invertci_single(x, alphas, msg.bound)
  } else {
    summary.invertci_multiple(x, alphas, msg.bound)
  }
}

#' Summary of results from \code{\link[lpinfer]{invertci}} for a single
#' significance level
#'
#' @description This function summarizes the results for
#'   \code{\link[lpinfer]{invertci}}.
#'
#' @inheritParams summary.invertci
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.invertci_single <- function(x, alphas, msg.bound, ...) {
  # ---------------- #
  # Step 1: Summary of results
  # ---------------- #
  cat(sprintf("Significance level: %s\n", round(x$alpha, digits = 5)))
  cat(sprintf("Confidence interval: [%s, %s]\n",
              round(x$ci[1, 3], digits = 5),
              round(x$ci[1, 4], digits = 5)))
  cat(sprintf("\nMaximum number of iterations: %s\n", x$max.iter[[1]][[1]]))
  cat(sprintf("Tolerance level: %s\n", x$tol))

  # Display the details if `f` is not `chorussell` (where x$df_ub and
  # x$df_lb will be NULL)
  if (!is.null(x$df_ub)) {
    cat("\n")
    # ---------------- #
    # Step 2: Messages in constructing the upper bound
    # ---------------- #
    cat("Details:\n\n")
    cat(sprintf(msg.bound, "upper"))
    consolidate.invertci(x$df_ub[[1]][[1]], x$termination$ub)
    cat("\n")

    # ---------------- #
    # Step 3: Messages in constructing the lower bound
    # ---------------- #
    cat(sprintf(msg.bound, "lower"))
    consolidate.invertci(x$df_lb[[1]][[1]], x$termination$lb)
  }
}

#' Summary of results from \code{\link[lpinfer]{invertci}} for multiple
#' significance levels
#'
#' @description This function summarizes the results for
#'   \code{\link[lpinfer]{invertci}}.
#'
#' @inheritParams summary.invertci_single
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.invertci_multiple <- function(x, alphas, msg.bound, ...) {
  # ---------------- #
  # Step 1: Extract the list of alphas that the user would like to check
  # ---------------- #
  alphas <- sort(alphas)
  # If null, then simply assign it as the alphas considered in the procedure
  if (is.null(alphas)) {
    alphas <- x$alpha
  } else {
    alphas <- alphas[which(x$alpha %in% alphas)]
  }

  # Do not display results if `alphas` provided do not match x$alphas
  if (length(which(x$alpha %in% alphas)) != 0) {
    # ---------------- #
    # Step 2: Print the confidence intervals and relevant parameters
    # ---------------- #
    # Print confidence interval
    print.invertci_multiple(x, alphas)
    cat("\n")

    # Print the relevant parameters
    cat(sprintf("Maximum number of iterations: %s\n", x$max.iter))
    cat(sprintf("Tolerance level: %s\n", x$tol))

    # ---------------- #
    # Step 3: Print the details of each iteration
    # ---------------- #
    # Display the details if `f` is not `chorussell` (where x$df_ub and
    # x$df_lb will be NULL)
    if (!is.null(x$df_ub)) {
      cat("\n")
      cat("Details:\n\n")
      for (i in seq_along(x$alpha)) {
        # Only print the result if alpha appears in alphas
        if (x$alpha[i] %in% alphas) {
          # Print the significance level
          cat(sprintf(paste0("********** Confidence interval for significance",
                             "level = %s **********"),
                      x$alpha[i]))

          for (j in 1:nrow(x$para.vals)) {
            # Print the parameters
            invertci.show.param(x$para.name, x$para.vals, j)

            # Upper bound
            cat(sprintf(msg.bound, "upper"))
            consolidate.invertci(x$df_ub[[i]][[j]],
                                 x$termination[[i]][[j]]$ub)

            # Lower bound
            cat(sprintf(msg.bound, "lower"))
            consolidate.invertci(x$df_lb[[i]][[j]],
                                 x$termination[[i]][[j]]$lb)
          }
          cat("\n")
        }
      }
    }
  }
}

#' Print results in constructing bounds in bisection method
#'
#' @description This function is used to display the message when constructing
#'    the bounds and used in \code{\link[lpinfer]{summary.invertci}} to print
#'    the results in each step of the bisection method.
#'
#' @inheritParams bisec.print
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.bisection.print <- function(df_bis, i) {
  # ---------------- #
  # Step 1: Data cleaning
  # ---------------- #
  print.iter1 <- df_bis[i, 1]
  if (print.iter1 == "Left end pt." | print.iter1 == "Right end pt.") {
    print.iter1 <- paste("", print.iter1, "\t")
  } else {
    print.iter1 <- paste("\r", as.character(print.iter1), "\t\t")
  }

  print.iter2 <- df_bis[i, 2]
  if (print.iter2 != "NA") {
    print.iter2 <- format(round(as.numeric(print.iter2), digits = 5),
                          nsmall = 5)
  } else {
    print.iter2 <- "NA\t"
  }

  print.iter3 <- df_bis[i, 3]
  if (print.iter3 != "NA") {
    print.iter3 <- format(round(as.numeric(print.iter3), digits = 5),
                          nsmall = 5)
  } else {
    print.iter3 <- "NA\t"
  }

  print.iter4 <- df_bis[i, 4]

  print.iter5 <- format(round(as.numeric(df_bis[i, 5]), digits = 5),
                        nsmall = 5)

  print.iter6 <- df_bis[i, 6]

  # ---------------- #
  # Step 2: Print results
  # ---------------- #
  cat(paste(print.iter1,
            print.iter2, "\t",
            print.iter3, "\t",
            format(round(print.iter4, digits = 5), nsmall = 5), "\t",
            print.iter5, "\t",
            print.iter6, "\n"))
}

#' Checks and updates the input of the function \code{\link[lpinfer]{invertci}}
#'
#' @description This function checks and updates the input from the user in the
#'    \code{\link[lpinfer]{invertci}} function. If there is any invalid input,
#'    the function will be terminated and error messages will be printed.
#'
#' @inheritParams invertci
#'
#' @return Returns the updated parameters back to the function
#' \code{subsample}. The following information are updated:
#'    \itemize{
#'       \item{\code{pvals}}
#'       \item{\code{lb0}}
#'       \item{\code{lb1}}
#'       \item{\code{ub0}}
#'       \item{\code{ub1}}
#'    }
#'
#' @export
#'
invertci.check <- function(f, farg, alpha, init.lb, init.ub, tol, max.iter,
                           pvals, progress) {
  # ---------------- #
  # Step 1: Conduct the checks
  # ---------------- #
  # Part 1. Check f
  if (class(f) != "function") {
    stop("The class of function ('f') has to be function.", call. = FALSE)
  }

  # Part 2. Check farg
  if (class(farg) != "list") {
    stop("The argument of the function ('farg') has to be a list.",
         call. = FALSE)
  }

  # Part 3. Check alpha
  for (i in seq_along(alpha)) {
    check.numrange(alpha[i], "alpha", "closed", 0, "closed", 1)
  }

  # Part 4. Check init.lb and init.ub
  # Note: lb and ub are not required if f is chorussell
  if (!identical(f, chorussell)) {
    ## Check init.lb and init.ub
    init.lb.return <- check.initb(init.lb, "init.lb", "lb")
    lb0 <- init.lb.return$lb
    lb1 <- init.lb.return$ub

    ## Check init.ub
    init.ub.return <- check.initb(init.ub, "init.ub", "ub")
    ub1 <- init.ub.return$lb
    ub0 <- init.ub.return$ub

    ## Retrieve the logical bounds if either init.lb or init.ub is empty
    if (is.null(init.ub) | is.null(init.lb)) {
      farg$beta.tgt <- 0
      freturn <- do.call(f, farg)
      logical.lb <- freturn$logical.lb
      logical.ub <- freturn$logical.ub

      ## Return an error if the logical bounds are infinite
      bd.infinite.msg <- paste0("The logical %s bound%s %s infinite. Please ",
                                "provide the corresponding initial bracket%s.")
      if (is.infinite(logical.lb) & is.infinite(logical.ub)) {
        stop(sprintf(bd.infinite.msg, "upper and lower", "s", "are", "s"))
      } else if (!is.infinite(logical.lb) & is.infinite(logical.ub)) {
        if (is.null(init.ub)) {
          stop(sprintf(bd.infinite.msg, "upper", "", "is", ""))
        }
      } else if (is.infinite(logical.lb) & !is.infinite(logical.ub)) {
        if (is.null(init.lb)) {
          stop(sprintf(bd.infinite.msg, "lower", "", "is", ""))
        }
      }
    }

    ## Assign the bounds if necessary
    if (is.null(lb0)) {
      lb0 <- logical.lb
    }

    if (is.null(ub0)) {
      ub0 <- logical.ub
    }

    if (is.null(lb1)) {
      lb1 <- ub0
    }

    if (is.null(ub1)) {
      ub1 <- lb0
    }
  } else {
    lb0 <- NULL
    lb1 <- NULL
    ub0 <- NULL
    ub1 <- NULL
  }

  # Part 7: Check tol
  check.positive(tol, "tol")

  # Part 8: Check max.iter
  check.positiveinteger(max.iter, "max.iter")

  # Part 9: Check pvals
  if (is.null(pvals) == TRUE) {
    # Part A: If pvals is null
    pvals <- NULL
  } else {
    # Part B: If pvals is non-null
    if (class(pvals) %in% c("data.frame", "matrix") == TRUE) {
      # Set pvals as a data frame
      pvals = as.data.frame(pvals)
      # Check the column names
      if (("point" %in% colnames(pvals)) & ("p-value" %in% colnames(pvals))) {
        stop("The column names of the data frame 'pvals' need to contain 'point'
             and 'p-value'.", call. = FALSE)
      }
      # Check if the values are numeric
      if (is.numeric(unlist(pvals)) == FALSE) {
        stop("The data frame 'pvals' has to be numeric.")
      }
      # Check if the p-values are bounded between [0, 1]
      if ((sum(pvals$`p-value` <= 1) != nrow(pvals)) |
          (sum(pvals$`p-value` >= 0) != nrow(pvals))) {
        stop("The p-values have to be in the interval [0,1].")
      }
    } else {
      stop(gsub("\\s+", " ",
                "The data provided 'pvals' must either be a data.frame,
                a data.table, or a matrix."), call. = FALSE)
    }
  }

  # Part 10: Check progress
  check.boolean(progress, "progress")

  # Step 11: Return the updated information
  return(list(pvals = pvals,
              lb0 = lb0,
              lb1 = lb1,
              ub0 = ub0,
              ub1 = ub1))
}

#' Consolidates and prints the \code{summary} table in
#' \code{\link[lpinfer]{invertci}}
#'
#' @description This function is used to print and consolidate the summary
#'   table for display via the \code{summary} command.
#'
#' @param df Data frame of the details for the bisection method.
#' @param msg Message that indicates the reason for the bisection method to
#'   terminate.
#'
#' @return Nothing is returned.
#'
#' @export
#'
consolidate.invertci <- function(df, msg) {
  df.temp <- df[, 2:6]
  df.temp[1, 2] <- NA
  df.temp[2, 1] <- NA
  for (i in 1:4) {
    df.temp[,i] <- formatC(as.numeric(df.temp[,i]), digits = 5, format = "f")
  }
  df.consol <- rbind(c("  Lower bound",
                       "  Upper bound",
                       "  Test point",
                       "  p-value",
                       "  Reject?"), df.temp)
  rownames(df.consol) <- c("Iteration", df[, 1])
  colnames(df.consol) <- NULL
  print(df.consol)
  cat(sprintf("Reason for termination: %s\n", msg))
}

#' Print the parameters used in the bisection method
#'
#' @description This function is used to print the tuning parameters that are
#'   relevant to the current set of iterations in the
#'   \code{\link[lpinfer]{invertci}} function.
#'
#' @param para.name The name(s) of the parameter(s).
#' @param para.vals The value(s) of the parameter(s).
#' @param j The row number for the parameters being used.
#'
#' @return Nothing is returned.
#'
#' @export
#'
invertci.show.param <- function(para.name, para.vals, j) {
  # Initialize the data.frame
  para.temp <- data.frame(matrix(vector(),
                                 ncol = length(para.name),
                                 nrow = 2))

  # Assign the parameters and the row names
  para.temp[1, ] <- para.name
  para.temp[2, ] <- round(para.vals[j, ], digits = 5)
  colnames(para.temp) <- NULL
  rownames(para.temp) <- c("Parameters:   ", "")

  # Print the parameters
  print(para.temp, digits = 5)
}
