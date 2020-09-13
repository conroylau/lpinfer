#' Constructs confidence interval
#'
#' @description This function constructs the confidence interval using the
#'    bisection method.
#'
#' @param f Function that represents a testing procedure.
#' @param farg List of arguments to be passed to the function of testing
#'    procedure.
#' @param alpha Significance level of the test.
#' @param lb0 Logical lower bound for the confidence interval. This is not
#'    required if the function \code{chorussell} is used.
#' @param lb1 Maximum possible lower bound for the confidence interval. This is
#'    not required if the function \code{chorussell} is used.
#' @param ub0 Logical upper bound for the confidence interval. This is not
#'    required if the function \code{chorussell} is used.
#' @param ub1 Minimum possible upper bound for the confidence interval. This
#'    is not required if the function \code{chorussell} is used.
#' @param tol Tolerance level in the bisection method.
#' @param max.iter Maximum number of iterations in the bisection method.
#' @param df_ci Data frame that consists the points and the corresponding
#'    \eqn{p}-values that have been tested in the previous iterations.
#' @param dp Number of decimal places in the output.
#' @param progress The boolean variable for whether the result messages should
#'    be displayed in the procedure of constructing confidence interval. If
#'    it is set as \code{TRUE}, the messages are displayed throughout the
#'    procedure. Otherwise, the messages will not be displayed.
#'
#' @returns Returns the confidence interval and a data frame that contains the
#'    points being tested in the procedure.
#'    \item{ub}{Upper bound of the confidence interval.}
#'    \item{lb}{Lower bound of the confidence interval.}
#'    \item{df_ci}{Data frame that consists of the points and the corresponding
#'       \eqn{p}-values that have been tested in constructing the confidence
#'       intervals.}
#'    \item{df_ub}{Data frame storing the information for the bisection
#'       method in each iteration when evaluating the upper bound of the
#'       confidence interval.}
#'    \item{df_lb}{Data frame storing the information for the bisection
#'       method in each iteration when evaluating the lower bound of the
#'       confidence interval.}
#'    \item{tol}{Tolerance level in the bisection method.}
#'    \item{iter}{Total number of iterations taken.}
#'    \item{call}{The function that has been called.}
#'
#' @details The number of decimal places displayed in the messages (if
#'    \code{progress} is set as \code{TRUE}) is equal to the number of decimal
#'    places in the variable \code{tol}.
#'
#' @example ./example/invertci_example.R
#'
#' @export
#'
invertci <- function(f, farg = list(), alpha = .05, lb0 = NULL, lb1 = NULL,
                     ub0 = NULL, ub1 = NULL, tol = .0001, max.iter = 20,
                     df_ci = NULL, dp = 5, progress = TRUE) {
  # ---------------- #
  # Step 1: Update call, check and update the arguments
  # ---------------- #
  # Extract the current RNG state
  rngstate <- .Random.seed

  # Obtain call information
  call <- match.call()

  # Check and update
  invertci.return <- invertci.check(f, farg, alpha, lb0, lb1, ub0, ub1, tol,
                                    max.iter, df_ci, progress)

  # Updates the input
  lb0 <- invertci.return$lb0
  lb1 <- invertci.return$lb1
  ub0 <- invertci.return$ub0
  ub1 <- invertci.return$ub1

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

    for (j in 1:nrow(pval)) {
      if (!identical(f, chorussell)) {
        # Assign the new farg object based on the j-th set of tuning parameters
        for (k in seq_along(para.name)) {
          farg[[para.name[k]]] <- para.vals[j, k]
        }

        # If `f` is not `chorussell`, use the bisection method in `invertci`
        if (i > 1 & (isTRUE(progress))) {
          cat("\n")
        }

        # Append the result that has already been evaluated in step 2
        df_ci <- invertci.return$df_ci
        df_ci[nrow(df_ci) + 1, 1] <- ub1
        df_ci[nrow(df_ci), 2] <- pval[j ,2]

        termination <- NULL

        ### Compute upper bound of confidence interval
        if (isTRUE(progress)) {
          # Print the significance level being considered
          cat(sprintf(paste0("********** Constructing confidence interval ",
                             "for significance level = %s **********"),
                      alpha[i]))

          # Print the parameters being considered
          invertci.show.param(para.name, para.vals, j)

          # Print the interval being considered
          cat(sprintf(comp.bound, "upper"))
        }
        ub_return <- ci.bisection(f, farg, alpha[i], ub1, ub0, tol, max.iter,
                                  df_ci, progress, 1, dp, rngstate)
        # Update data frame
        df_ci <- ub_return$df_ci
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
                                  df_ci, progress, -1, dp, rngstate)
        # Update data frame
        df_ci <- lb_return$df_ci
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
      } else {
        # If `f` is `chorussell`, use the bisection method in `chorussell`
        farg$tol <- tol
        farg$alpha <- alpha[i]
        farg$ci <- TRUE
        chorussell.return <- do.call(chorussell, farg)
        lb <- chorussell.return$ci.df[1, 2]
        ub <- chorussell.return$ci.df[1, 3]
        iter <- chorussell.return$iter

        # Consolidate the results
        lb_list[[i]][[j]] <- lb
        ub_list[[i]][[j]] <- ub
        iter_list[[i]][[j]] <- iter
      }
    }
  }

  # ---------------- #
  # Step 5: Consolidate the data frame of confidence intervals
  # ---------------- #
  # Initialize data frame
  nr.p <- nrow(para.vals)
  nc.p <- ncol(para.vals)
  consol.ci <- data.frame(matrix(vector(),
                                 nrow = nr.p * length(alpha),
                                 ncol = ncol(para.vals) + 2))
  colnames(consol.ci) <- c("alpha", para.name, "ci")

  # Consolidate data
  for (i in seq_along(alpha)) {
    ## Append the alphas
    consol.ci[(nr.p * (i - 1) + 1):(nr.p * i), 1] <- alpha[i]

    ## Append the tuning parameters
    consol.ci[(nr.p * (i - 1) + 1):(nr.p * i), 2:(nc.p + 1)] <- para.vals

    ## Append the confidence intervals
    for (j in 1:nr.p) {
      consol.ci[j + (i - 1) * nr.p, ncol(consol.ci)] <-
        sprintf("[%s, %s]",
                lb_list[[i]][[j]],
                ub_list[[i]][[j]])
    }
  }

  # ---------------- #
  # Step 6: Assign the return list and return output
  # ---------------- #
  output <- list(ub = ub_list,
                 lb = lb_list,
                 df_ci = df_ci,
                 tol = tol,
                 alpha = alpha,
                 iter = iter_list,
                 call = call,
                 consol.ci = consol.ci)

  # The following information are returned only if `f` is not `chorussell`
  if (!identical(f, chorussell)) {
    output$max.iter <- max.iter
    output$df_ub <- df_ub_list
    output$df_lb <- df_lb_list
    output$termination <- termination_list
    output$para.name <- para.name
    output$para.vals <- para.vals
  }
  attr(output, "class") <- "invertci"

  # Return output
  return(output)
}

#' Bisection method for constructing confidence intervals
#'
#' @description This function constructs the two-sided confidence interval
#'    of a given testing procedure using the bisection method.
#'
#' @param type Type of the confidence interval. Set \code{type} as 1 for the
#'    upper bound of the confidence interval and set \code{type} as -1 for the
#'    lower bound of the confidence interval.
#' @param b0 Logical lower or upper bound for the confidence interval.
#' @param b1 Maximum possible lower bound or minimum possible upper bound for
#'    the confidence interval.
#' @param dp Number of decimal places to be displayed for the \eqn{p}-values
#'    and confidence intervals in the messages if \code{progress} is set
#'    as \code{TRUE}.
#' @param rngstate The current RNG state obtained from \code{.Random.seed}.
#' @inheritParams invertci
#'
#' @return Return the solution of the bisection method and the updated
#'    data frame.
#'    \item{soln}{Solution to the bisection method.}
#'    \item{df_ci}{Data frame that consists of the points and the
#'       corresponding \eqn{p}-values that have been tested in constructing the
#'       confidence intervals.}
#'    \item{df_bis}{Data frame storing the information for the bisection
#'       method in each iteration.}
#'
#' @export
#'
ci.bisection <- function(f, farg, alpha, b0, b1, tol, max.iter, df_ci,
                         progress, type, dp, rngstate) {

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
  fb0_return <- bisec.eval(f, farg, a, df_ci, rngstate)
  df_ci <- fb0_return$df_ci
  # Print information
  if (isTRUE(progress)) {
    cat(paste0(" Iteration\t Lower bound \t Upper bound \t Test point \t ",
               "p-value\t Reject?\n"))
  }
  df_bis <- bisec.print("left end", alpha_2sided, fb0_return, a, "NA",
                        progress, dp, df_bis)$df_bis

  ### Right end-point b
  b <- b1
  fb1_return <- bisec.eval(f, farg, b, df_ci, rngstate)
  df_ci <- fb1_return$df_ci
  # Print information
  df_bis <- bisec.print("right end", alpha_2sided, fb1_return, "NA", b,
                        progress, dp, df_bis)$df_bis

  # If fb1 and fb0 are of the same sign, ask user to choose another interval
  # Compute mid-point and evaluate the corresponding p-value
  c <- (b + a)/2
  fc_return <- bisec.eval(f, farg, c, df_ci, rngstate)
  fc <- fc_return$pval
  df_ci <- fc_return$df_ci

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
    fc_return <- bisec.eval(f, farg, c, df_ci, rngstate)
    fc <- fc_return$pval
    df_ci <- fc_return$df_ci
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
                 df_ci = df_ci,
                 df_bis = df_bis,
                 last_iter_msg = last_iter_msg))
}

#' Evaluation of test statistic and check if the point has been evaluated
#'
#' @description This function checks if the \eqn{p}-value for the point
#'    considered has already been evaluated in previous iterations or provided
#'    by the user. The function will compute the \eqn{p}-value if it has been
#'    evaluated. Otherwise, it will use the previous data.
#'
#' @inheritParams invertci
#' @inheritParams ci.bisection
#' @param pt Point to be evaluated in the bisection method.
#'
#' @return Returns the \eqn{p}-value of the point considered and an updated
#'    data frame that contains the points and the \eqn{p}-values.
#'    \item{pval}{\eqn{p}-value of the point.}
#'    \item{df_ci}{Updated data frame that consists of the points that have
#'       been tested in constructing the confidence intervals.}
#'
#' @export
#'
bisec.eval <- function(f, farg, pt, df_ci, rngstate) {
  # ---------------- #
  # Step 1: Check if the data point has appeared in previous iterations.
  # ---------------- #
  df_match <- df_ci[df_ci[, "point"] == pt,]
  df_n <- nrow(df_ci)

  # ---------------- #
  # Step 2: If the points have been evaluated in df_ci, it will be taken
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
    df_ci[df_n + 1, 1] <- pt
    df_ci[df_n + 1, 2] <- pval
  } else {
    pval <- df_match[2]
  }

  # ---------------- #
  # Step 3: Return results
  # ---------------- #
  return(list(pval = pval,
              df_ci = df_ci))
}

#' Determine whether a point is inside the confidence interval or not
#'
#' @description This function determines whether the bisection method is going
#'    to be updated by choosing the left segment or the right segment of the
#'    interval as the updated interval for the next iteration in the bisection
#'    method.
#'
#' @param pval \eqn{p}-value of the test statistic.
#' @inheritParams invertci
#' @inheritParams ci.bisection
#'
#' @return Returns whether the part of the interval to be selected in the
#'    next iteration of the bisection method.
#'    \item{part}{Left or right segment of the interval.}
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
#' @description This function prints the information in \code{ci\_bisection}
#'    and store the result for each iteration.
#'
#' @param procedure Variable indicating whether the function is evaluating
#'    the end-points or first mid-point, or is iterating through the bisection
#'    procedure.
#' @param alphahalf Half of significance value that is used to evaluate the
#'    confidence interval.
#' @param returnlist List of information obtained from running
#'    \code{bisec.eval}.
#' @param a Lower bound of the current interval. This is \code{NULL} if the
#'    initial end-points are being evaluated.
#' @param b Upper bound of the current interval. This is \code{NULL} if the
#'    initial end-points are being evaluated.
#' @param df_bis Data frame storing the information from the bisection method.
#' @inheritParams invertci
#'
#' @return Return the updated data frame that stores the information for the
#'    the iteration.
#'    \item{df_bis}{Data frame storing the information for the bisection
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

#' Print results from \code{invertci}
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#'
#' @return Print the basic set of results from \code{invertci}.
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

#' Print results from \code{invertci} with a single alpha
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#'
#' @return Print the basic set of results from \code{invertci}.
#'
#' @export
#'
print.invertci_single <- function(x, ...) {
  cat(sprintf("Confidence interval: [%s, %s]\n",
              round(x$lb[[1]][[1]], digits = 5),
              round(x$ub[[1]][[1]], digits = 5)))
}

#' Print results from \code{invertci} with vector-valued alpha
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#'
#' @return Print the basic set of results from \code{invertci}.
#'
#' @export
#'
print.invertci_multiple <- function(x, alphas = NULL, ...) {
  cat("Confidence intervals:\n")
  df <- x$consol.ci
  colnames(df)[1] <- "Significance level"
  colnames(df)[ncol(df)] <- "Confidence interval"
  if (is.null(alphas)) {
    print(df, row.names = FALSE, digits = 5)
  } else {
    print(df[df[, 1] %in% alphas, ], row.names = FALSE, digits = 5)
  }
}

#' Summary of results from \code{invertci}
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param alphas List of alphas that the user would like to print.
#' @param ... Additional arguments.
#'
#' @return Print the summary of the basic set of results from \code{invertci}.
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

#' Summary of results from \code{invertci} for a single alpha
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param msg.bound String containing the general message to be printed.
#' @param ... Additional arguments.
#' @inheritParams invertci
#'
#' @return Print the summary of the basic set of results from \code{invertci}.
#'
#' @export
#'
summary.invertci_single <- function(x, alphas, msg.bound, ...) {
  # ---------------- #
  # Step 1: Summary of results
  # ---------------- #
  cat(sprintf("Significance level: %s\n", round(x$alpha, digits = 5)))
  cat(sprintf("Confidence interval: [%s, %s]\n",
              round(x$lb[[1]][[1]], digits = 5),
              round(x$ub[[1]][[1]], digits = 5)))
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

#' Summary of results from \code{invertci} with multiple alphas
#'
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#' @inheritParams summary.invertci
#' @inheritParams summary.invertci_single
#'
#' @return Print the summary of the basic set of results from \code{invertci}.
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
#'    the bounds and used in \code{summary.invertci} to print the results in
#'    each step of the bisection method.
#'
#' @inheritParams bisec.print
#'
#' @return Nothing is returned in this function.
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

#' Checks and updates the input of the function \code{invertci}
#'
#' @description This function checks and updates the input from the user for
#'    the function \code{invertci}. If there is any invalid input, this
#'    function ill terminate the procedure and generate appropriate error
#'    messages.
#'
#' @inheritParams invertci
#'
#' @return Returns the updated parameters back to the function
#' \code{subsample}. The following information are updated:
#'    \itemize{
#'       \item{\code{df_ci}}
#'       \item{\code{lb0}}
#'       \item{\code{lb1}}
#'       \item{\code{ub0}}
#'       \item{\code{ub1}}
#'    }
#'
#' @export
#'
invertci.check <- function(f, farg, alpha, lb0, lb1, ub0, ub1, tol, max.iter,
                           df_ci, progress) {
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

  # lb and ub are not required if f is chorussell
  if (!identical(f, chorussell)) {
    # Part 4: Check lb0 and ub0
    farg$beta.tgt <- 0
    freturn <- do.call(f, farg)
    ## Lower bound
    if (is.null(lb0)) {
      lb0 <- freturn$logical.lb
    } else {
      if (!(is.numeric(lb0) == TRUE & length(lb0) == 1)) {
        stop("The argument 'lb0' must be a scalar.", call. = FALSE)
      }
    }
    ## Upper bound
    if (is.null(ub0)) {
      ub0 <- freturn$logical.ub
    } else {
      if (!(is.numeric(ub0) == TRUE & length(ub0) == 1)) {
        stop("The argument 'ub0' must be a scalar.", call. = FALSE)
      }
    }

    # Part 5: Check lb1 and ub1
    # Part A: Check lb1
    if (is.null(lb1)) {
      # If lb1 is null, assign lb1 as ub0
      lb1 = ub0
    } else {
      # If lb1 is nonnull, check whether its numeric
      check.numeric(lb1, "lb1")
    }
    # Part B: Check ub1
    if (is.null(ub1)) {
      # If ub1 is null, assign ub1 as lb0
      ub1 = lb0
    } else {
      # If ub1 is nonnull, check whether its numeric
      check.numeric(lb0, "lb0")
    }

    # Part 6: Check the difference between lb0 vs lb1, and ub0 vs ub1
    if (lb0 > lb1) {
      stop("The logical lower bound 'lb0' cannot be larger than the maximum
         possible lower bound 'lb1'.")
    }
    if (ub0 < ub1) {
      stop("The logical upper bound 'ub0' cannot be smaller than the minimum
         possible upper bound 'ub1'.")
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

  # Part 9: Check df_ci
  if (is.null(df_ci) == TRUE) {
    # Part A: If df_ci is null
    df_ci <- data.frame(matrix(vector(), 0, 2,
                               dimnames = list(c(), c("point", "value"))),
                        stringsAsFactors = F)
  } else {
    # Part B: If df_ci is non-null
    if (class(df_ci) %in% c("data.frame", "matrix") == TRUE) {
      # Set df_ci as a data frame
      df_ci = as.data.frame(df_ci)
      # Check the column names
      if (sum(colnames(df_ci) == c("point", "value")) != 2) {
        stop("The column names of the data frame 'df_ci' have to be 'point'
             and 'value'.", call. = FALSE)
      }
      # Check if the values are numeric
      if (is.numeric(unlist(df_ci)) == FALSE) {
        stop("The data frame 'df_ci' has to be numeric.")
      }
      # Check if the p-values are bounded between [0, 1]
      if ((sum(df_ci[, 2] <= 1) != nrow(df_ci)) |
          (sum(df_ci[, 2] >= 0) != nrow(df_ci))) {
        stop("The p-values have to be in the interval [0,1].")
      }
    } else {
      stop(gsub("\\s+", " ",
                "The data provided 'df_ci' must either be a data.frame,
                a data.table, or a matrix."), call. = FALSE)
    }
  }

  # Part 10: Check progress
  check.boolean(progress, "progress")

  # Step 11: Return the upated information
  return(list(df_ci = df_ci,
              lb0 = lb0,
              lb1 = lb1,
              ub0 = ub0,
              ub1 = ub1))
}

#' Function to consolidate the print the summary table for the interactions
#'
#' @description This function is used to consolidate the summary table for
#'   display via the \code{summary} command.
#'
#' @inheritParams bisec.print
#'
#' @return Returns the consolidate table.
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
#'   relevant to the current set of iterations in the \code{invertci} function.
#'
#' @param para.name Name(s) of the parameter(s).
#' @param para.vals Value(s) of the parameter(s).
#' @param j Row number for the parameters being used.
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
