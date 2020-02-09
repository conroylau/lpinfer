#' Constructs confidence interval
#' 
#' @description This function constructs the confidence interval using the
#'    bisection method. 
#'    
#' @import parallel foreach doMC
#' 
#' @param f Function that represents a testing procedure.
#' @param farg List of arguments to be passed to the function of testing 
#'    procedure.
#' @param alpha Significance level of the test.
#' @param lb0 Logical lower bound for the confidence interval.
#' @param lb1 Maximum possible lower bound for the confidence interval.
#' @param ub0 Logical upper bound for the confidence interval.
#' @param ub1 Minimum possible upper bound for the confidence interval.
#' @param tol Tolerance level in the bisection method.
#' @param max_iter Maximum number of iterations in the bisection method.
#' @param df_ci Data frame that consists the points and the corresponding
#'    \eqn{p}-values that have been tested in the previous iterations.
#' @param progress The boolean variable for whether the result messages should
#'    be displayed in the procedure of constructing confidence interval. If 
#'    it is set as \code{TRUE}, the messages are displayed throughout the 
#'    procedure. Otherwise, the messages will not be displayed.
#' 
#' @returns Returns the confidence interval and a data frame that contains the 
#'    points being tested in the procedure.
#'    \item{up}{Upper bound of the confidence interval.}
#'    \item{low}{Lower bound of the confidence interval.}
#'    \item{df_ci}{data frame that consists of the points and the corresponding
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
#' @export 
#' 
invertci <- function(f, farg, alpha = .05, lb0 = NULL, lb1 = NULL, ub0 = NULL, 
                 ub1 = NULL, tol = .0001, max_iter = 20, df_ci = NULL,
                 progress = FALSE){
  
  #### Step 1: Obtain call, check and update the input
  # Obtain call information
  call = match.call()
  # Check and update
  check_return = invertci_check(f, farg, alpha, lb0, lb1, ub0, ub1, tol, 
                                max_iter, df_ci, progress)
  # Updates the input
  lb0 = check_return$lb0
  lb1 = check_return$lb1
  ub0 = check_return$ub0
  ub1 = check_return$ub1
  # Compute the number of decimal places in tol
  dp = decimal_places(tol)
  # Initialize lists
  df_ub_list = NULL
  df_lb_list = NULL
  termination_list = NULL
  ub_list = NULL
  lb_list = NULL
  iter_list = NULL
  
  #### Step 2: Return confidence interval and data frame
  
  for (i in 1:length(alpha)){
    df_ci = check_return$df_ci
    termination = NULL
    
    cat(sprintf("< Constructing confidence interval for alpha = %s >\n", 
                alpha[i]))
    ### Compute upper bound of confidence interval
    if (progress == TRUE){
      cat("\n=== Computing upper bound of confidence interval ===\n")
    }
    ub_return = ci_bisection(f, farg, alpha[i], ub1, ub0, tol, max_iter, df_ci, 
                             progress, 1, dp)
    # Update data frame
    df_ci = ub_return$df_ci
    # Data frame storing all messages in each iteration
    df_ub = ub_return$df_bis
    # Obtain termination message
    termination$ub = ub_return$last_iter_msg
    # Obtain upper bound
    ub = ub_return$pt
    
    ### Compute lower bound of confidence interval
    if (progress == TRUE){
      cat("\n=== Computing lower bound of confidence interval ===\n")
    }
    lb_return = ci_bisection(f, farg, alpha[i], lb0, lb1, tol, max_iter, df_ci, 
                             progress, -1, dp)
    # Update data frame
    df_ci = lb_return$df_ci
    # Data frame storing all messages in each iteration
    df_lb = lb_return$df_bis
    # Obtain termination message
    termination$lb = lb_return$last_iter_msg
    # Obtain lower bound
    lb = lb_return$pt
    
    # Iterations 
    iter = lb_return$iter + ub_return$iter
    
    #### Step 3: Print confidence interval and return results
    if (progress == TRUE){
      cat("-----------------------------------\n")
      cat(paste("< Significance level: ", alpha[i], " >\n", sep = ""))
      cat(paste("Tolerance level: ", tol, ".\n", sep = ""))
      cat(paste("Total number of iterations: ", iter, ".\n", sep = ""))
      cat(paste("Confidence interval: [", 
                round(lb, digits = dp), ", ", 
                round(ub, digits = dp), "].\n\n", sep = "")) 
    }
    
    #### Step 4: Store information as list if alphas is a vector
    if (length(alpha) > 1){
      df_ub_list = c(df_ub_list, list(df_ub))
      df_lb_list = c(df_lb_list, list(df_lb))
      termination_list = c(termination_list, list(termination))
      lb_list = c(lb_list, list(lb))
      ub_list = c(ub_list, list(ub))
      iter_list = c(iter_list, list(iter))
    }
  }
  
  #### Step 5: Rename variables for output
  if (length(alpha) > 1){
    df_lb = df_lb_list
    df_ub = df_ub_list
    termination = termination_list
    lb = lb_list
    ub = ub_list
    iter = iter_list
  }
  
  # Assign the return list
  output = list(ub = ub,
                lb = lb,
                df_ci = df_ci,
                df_ub = df_ub,
                df_lb = df_lb,
                tol = tol,
                alpha = alpha,
                iter = iter,
                call = call,
                termination = termination)
  attr(output, "class") = "invertci"
  
  invisible(output)
}

#' Bisection method
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
ci_bisection <- function(f, farg, alpha, b0, b1, tol, max_iter, df_ci,
                         progress, type, dp){
  
  #### Step 1: Evaluate the end-points and the mid-point of b0 and b1
  # Initialize data frame to collect the information in the bisection method
  df_bis = data.frame(matrix(vector(), 0, 6,
                             dimnames=list(c(), 
                                           c("iteration", 
                                             "left",
                                             "right",
                                             "point",
                                             "p-value",
                                             "decision"))),
                      stringsAsFactors=F)

  # Divide alpha by 2
  alpha_2sided = alpha/2
  ## Left end-point a
  a = b0
  fb0_return = bisec_eval(f, farg, a, df_ci)
  fb0 = fb0_return$pval
  df_ci = fb0_return$df_ci
  # Print information
  cat("Iteration\t Lower bound \t Upper bound \t Test point \t p-value\t Reject?\n")
  df_bis = bisec_print("left end", alpha_2sided, fb0_return, a, "NA", progress, 
                       dp, df_bis)$df_bis
  
  ## Right end-point b
  b = b1
  fb1_return = bisec_eval(f, farg, b, df_ci)
  fb1 = fb1_return$pval
  df_ci = fb1_return$df_ci
  # Print information
  df_bis = bisec_print("right end", alpha_2sided, fb1_return, "NA", b, 
                       progress, dp, df_bis)$df_bis

  # If fb1 and fb0 are of the same sign, ask user to choose another interval
  if ((fb1 - alpha_2sided) * (fb0 - alpha_2sided) > 0){
    stop("Please choose another interval.")
  }
  # Compute mid-point and evaluate the corresponding p-value
  c = (b+a)/2
  fc_return = bisec_eval(f, farg, c, df_ci)
  fc = fc_return$pval
  df_ci = fc_return$df_ci  

  #### Step 2: Bisection method
  for (i in 1:max_iter){
    # Bisection method complete if the difference between the two points is
    # below the tolereance level.
    if (abs(b-a) < tol){
      tol_msg = paste(">>> Length of interval is below tolerance level. ",
                      "Bisection method is completed.\n", sep = "")
      last_iter_msg = ">>> Length of interval is below tolerance level"
      if (progress == TRUE){
        cat(tol_msg) 
      }
      break
    }
    # Update interval based on whether the left section or the section of the 
    # interval is chosen
    # Print information
    df_bis = bisec_print(i , alpha_2sided, fc_return, a, b, progress,
                         dp, df_bis)$df_bis
    
    if (ci_inout(fc, alpha_2sided, type) == "left"){
      b = c
    } else {
      a = c
    }
    # Evaluate new mid-point
    c = (a+b)/2
    # Update data frame and p-value
    fc_return = bisec_eval(f, farg, c, df_ci)
    fc = fc_return$pval
    df_ci = fc_return$df_ci  
  }
  
  # Only called when the maximum number of iterations is reached
  iter_msg = paste(">>> Reached the maximum number of iterations. ",
                   "Bisection method is completed.\n", sep = "") 
  last_iter_msg = ">>> Reached maximum number of iterations"
  if (progress == TRUE & i == max_iter){
    cat(iter_msg) 
    
  }
  
  #### Step 3: Return results
  # Return the results
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
#' @param pt Point to be evaluated in the bisection method.
#' 
#' @return Returns the \eqn{p}-value of the point considered and an updated
#'    data frame that contains the points and the \eqn{p}-values.
#'    \item{pval}{\eqn{p}-value of the point.}
#'    \item{df_ci}{Updated data frame that consists of the points that have been 
#'       tested in constructing the confidence intervals.}
#' 
#' @export
#' 
bisec_eval <- function(f, farg, pt, df_ci){
  #### Step 1: Check if the data point has appeared in previous iterations.
  df_match = df_ci[df_ci[, "point"] == pt,]
  df_n = nrow(df_ci)
  
  #### Step 2: If such point has been evaluated in the past, then it is 
  #### obatined from df_ci. Otherwise, it will be evaluated.
  if ((is.null(df_match) == TRUE) | (dim(df_match)[1] == 0)){
    farg$beta_tgt = pt
    test_return = do.call(f, farg)
    pval = test_return$p_val
    df_ci[df_n+1, 1] = pt
    df_ci[df_n+1, 2] = pval
  } else {
    pval = df_match[2]
  }

  #### Step 3: Return results
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
#' @inheritParams ci_bisection
#' 
#' @return Returns whether the part of the interval to be selected in the 
#'    next iteration of the bisection method.
#'    \item{part}{Left or right segment of the interval.}
#'    
#' @export
#' 
ci_inout <- function(pval, alpha, type){
  if (type == 1){
    ### Type == 1: Upper bound
    if (pval < alpha){
      part = "left"
    } else {
      part = "right"
    }
  } else if (type == -1){
    ### Type == -1: Lower bound
    if (pval < alpha){
      part = "right"      
    } else {
      part = "left"
    }
  }
  return(part)
}

#' Checks and updates the input of the function \code{invertci}
#' 
#' @description This function checks and updates the input from the user for 
#'    the function \code{invertci}. If there is any invalid input, this function 
#'    will terminate the procedure and generate appropriate error messages.
#'    
#' @inheritParams invertci
#' 
#' @return Returns the updated value of the parameters back to the function 
#'    \code{invertci} in the correct format.
#'    \item{df_ci}{Data frame that stores the points that has been evaluated
#'       and the corresponding \eqn{p}-values.}
#'    \item{lb0}{Logical lower bound for the confidence interval.}
#'    \item{lb1}{Maximum possible lower bound for the confidence interval.}
#'    \item{ub0}{Logical upper bound for the confidence interval.}
#'    \item{ub1}{Minimum possible upper bound for the confidence interval.}
#' 
#' @export
#' 
invertci_check <- function(f, farg, alpha, lb0, lb1, ub0, ub1, tol, max_iter, 
                       df_ci, progress){

  #### Part 1. Check f
  if (class(f) != "function"){
    stop("The class of function ('f') has to be function.", call. = FALSE)
  }    

  #### Part 2. Check farg
  if (class(farg) != "list"){
    stop("The arguemnt of the function ('farg') has to be a list.", 
         call. = FALSE)
  }   
  
  #### Part 3. Check alpha
  if ((is.numeric(alpha) == TRUE & sum(alpha >= 0) == length(alpha)
       & sum(alpha <= 1) == length(alpha)) == FALSE){
    stop("The significance level ('alpha') has to be in the interval [0,1].", 
         call. = FALSE)
  }    
  
  #### Part 4: Check lb0 and ub0
  if (is.null(lb0) | is.null(ub0)){
    ### Part A: If lb0 and ub0 are null, check whether if the function is dkqs
    ### If yes, compute lb0 and ub0. Otherwise, return terminate the function.
    if (as.character(substitute(f)) == "dkqs"){
      logicalb_return = dkqs_logicalb(f, farg)
      if (is.null(lb0)){
        lb0 = logicalb_return$lb0
      }
      if (is.null(ub0)){
        ub0 = logicalb_return$ub0
      }
    } else{
      if (is.null(ub0) | !is.null(lb0)){
        stop("Please provide the logical upper bound 'ub0'.", call. = FALSE) 
      } else if (is.null(lb0) | !is.null(ub0)){
        stop("Please provide the logical lower bound 'lb0'.", call. = FALSE) 
      } else {
        stop("Please provide the logical lower bound 'lb0' and the logical
             upper bound 'ub0'.", call. = FALSE) 
      }
    }
  } else{
    ### Part B: If lb0 and ub0 are both nonnull, check whether they are numeric.
    if (!(is.numeric(lb0) == TRUE & length(lb0) == 1)) {
      stop("The argument 'lb0' must be a scalar.", call. = FALSE)
    }  
    if (!(is.numeric(ub0) == TRUE & length(ub0) == 1)) {
      stop("The argument 'ub0' must be a scalar.", call. = FALSE)
    }  
  }

  #### Part 5: Check lb1 and ub1
  ### Part A: Check lb1
  if (is.null(lb1)){
    # If lb1 is null, assign lb1 as ub0
    lb1 = ub0
  } else {
    # If lb1 is nonnull, check whether its numeric
    if (!(is.numeric(lb1) == TRUE & length(lb1) == 1)) {
      stop("The argument 'lb1' must be a scalar.", call. = FALSE)
    } 
  }
  ### Part B: Check ub1
  if (is.null(ub1)){
    # If ub1 is null, assign ub1 as lb0
    ub1 = lb0
  } else {
    # If ub1 is nonnull, check whether its numeric
    if (!(is.numeric(ub1) == TRUE & length(ub1) == 1)) {
      stop("The argument 'ub1' must be a scalar.", call. = FALSE)
    } 
  }
  
  #### Part 6: Check the difference between lb0 vs lb1, and ub0 vs ub1
  if (lb0 > lb1){
    stop("The logical lower bound 'lb0' cannot be larger than the maximum 
         possible lower bound 'lb1'.")
  }
  if (ub0 < ub1){
    stop("The logical upper bound 'ub0' cannot be smaller than the minimum 
         possible upper bound 'ub1'.")
  }
  
  #### Part 7: Check tol
  if ((is.numeric(tol) == TRUE & length(tol) == 1 & tol > 0) == FALSE){
    stop("The tolerance level ('tol') must be a positive scalar.",
         call. = FALSE)
  }  
  
  #### Part 8: Check max_iter
  if ((is.numeric(max_iter) == TRUE & length(max_iter) == 1 & max_iter >= 0 &
       max_iter%%1 == 0) == FALSE){
    stop("The number of iterations ('max_iter') must be a positive integer.",
         call. = FALSE)
  }
  
  #### Part 9: Check df_ci
  if (is.null(df_ci) == TRUE){
    ### Part A: If df_ci is null
    df_ci = data.frame(matrix(vector(), 0, 2,
                              dimnames=list(c(), c("point", "value"))),
                       stringsAsFactors=F)
  } else {
    ### Part B: If df_ci is non-null
    if (class(df_ci) %in% c("data.frame", "matrix") == TRUE){
      # Set df_ci as a data frame
      df_ci = as.data.frame(df_ci)  
      # Check the column names
      if (sum(colnames(df_ci) == c("point", "value")) != 2){
        stop("The column names of the data frame 'df_ci' have to be 'point' 
             and 'value'.", call. = FALSE)
      }
      # Check if the values are numeric
      if (is.numeric(unlist(df_ci)) == FALSE){
        stop("The data frame 'df_ci' has to be numeric.")
      }
      # Check if the p-values are bounded between [0, 1]
      if ((sum(df_ci[,2] <= 1) != nrow(df_ci)) | 
          (sum(df_ci[,2] >= 0) != nrow(df_ci))){
        stop("The p-values have to be in the interval [0,1].")
      }      
      
    } else {
      stop(gsub("\\s+", " ",
                "The data povided 'df_ci' must either be a data.frame, 
                a data.table, or a matrix."), call. = FALSE)    
    }
  }

  #### Part 10: Check progress
  if (!(progress == TRUE | progress == FALSE)){
    stop("The argument 'progress' has to be boolean.")
  }
  
  ### Step 11: Return the upated information
  return(list(df_ci = df_ci,
              lb0 = lb0,
              lb1 = lb1,
              ub0 = ub0,
              ub1 = ub1))
}

#' Compute the logical upper and lower bounds for dkqs
#' 
#' @description This function computes the logical upper and lower bounds for
#'    the test \code{dkqs}.
#'    
#' @inheritParams invertci
#' 
#' @return Returns the logical upper and lower bounds for dkqs.
#'    \item{lb0}{Logical lower bound for \code{dkqs}.}
#'    \item{ub0}{Logical upper bound for \code{dkqs}.}
#' 
#' @export
#' 
dkqs_logicalb <- function(f, farg){
  
  #### Step 1: Assign a value to beta_tgt for returning the results
  farg$beta_tgt = 0
  
  #### Step 2: Run dkqs to obtain the logical bounds
  dkqs_return = do.call(f, farg)
  lb0 = dkqs_return$lb0
  ub0 = dkqs_return$ub0

  #### Step 3: Return results
  invisible(list(lb0 = lb0, 
                 ub0 = ub0))
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
#'    \code{bisec_eval}.
#' @param a Lower bound of the current interval. This is \code{NULL} if the
#'    initial end-points are being evaluated.
#' @param b Upper bound of teh current interval. This is \code{NULL} if the
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
bisec_print <- function(procedure, alphahalf, returnlist, a, b, progress, dp,
                        df_bis){
  
  #### Step 1: Obtain information about the current data frame
  df_bis_row = nrow(df_bis)
  space6 = "      "
  # Update decision
  if (returnlist$pval < alphahalf){
    decision = TRUE
  } else {
    decision = FALSE
  }
  
  
  # The messages are being displayed only if 'progress' is set to TRUE
  if (progress == TRUE){
    #### Step 2: Print the iteration
    if (is.numeric(procedure) == FALSE){
      # Case A: 'procedure' is not numeric if evaluating the initial 3 points
      if (procedure == "left end"){
        # cat(paste0(space6, "* Point being evaluated: ", 
        #            round(a, digits = dp), "\n"))  
        df_bis[df_bis_row + 1, 4] = a
      } else if (procedure == "right end"){
        df_bis[df_bis_row + 1, 4] = b
      }
    } else {
      # Case B: 'procedure' is numeric if evaluating the bisection method
      df_bis[df_bis_row + 1, 4] = (a+b)/2
    }
    
  }
  
  #### Step 5: Update data frame
  # Update column 1, i.e. whether evaluating end-points or iterations
  if (procedure == "left end"){
    df_bis[df_bis_row + 1,1] = "Left end pt."
  } else if (procedure == "right end"){
    df_bis[df_bis_row + 1,1] = "Right end pt."
  } else if (is.numeric(procedure) == TRUE){
    df_bis[df_bis_row + 1,1] = procedure
  }
  
  df_bis[df_bis_row + 1, 2] = a
  df_bis[df_bis_row + 1, 3] = b
  df_bis[df_bis_row + 1, 5] = returnlist$pval
  df_bis[df_bis_row + 1, 6] = decision

  #### Step 6: Print information
  summary_bisection_print(df_bis, df_bis_row + 1)
  
  #### Step 7: Return information
  invisible(list(df_bis = df_bis))
}

#' Auxiliary function: Count decimal places
#' 
#' @description This function returns the number of decimal places of a real 
#'    snumber.
#' 
#' @param x Real number.
#' 
#' @return Returns the number of decimal places of the number.
#'    \item{dp}{Number of decimal places of the number \code{x}.}
#'    
#' @export
#' 
decimal_places <- function(x){
  if ((x %% 1) == 0){
    #### Case 1: x is an integer
    invisible(0)
  } else {
    #### Case 2: x is not an integer
    x = as.character(x)
    xsplit = strsplit(x, ".", fixed = TRUE)
    # Count the number of digits after the dot
    dp = nchar(xsplit[[1]][2])
    invisible(dp) 
  }
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
print.invertci <- function(x, ...){
  if (length(x$alpha) == 1){
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
print.invertci_single <- function(x, ...){
  cat(sprintf("< Significance level: %s >\n", round(x$alpha, digits = 5)))
  cat(sprintf("Total number of iterations: %s.\n", round(x$iter, digits = 5)))  
  cat(sprintf("Confidence interval: [%s, %s].\n", 
              round(x$lb, digits = 5),
              round(x$ub, digits = 5)))
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
print.invertci_multiple <- function(x, ...){
  for (i in 1:length(x$alpha)){
    if (i > 1){
      cat("--------------------------------------\n")
    }
    cat(sprintf("< Significance level: %s >\n", round(x$alpha[i], digits = 5)))
    cat(sprintf("Total number of iterations: %s.\n", round(x$iter[[i]], 
                                                           digits = 5)))  
    cat(sprintf("Confidence interval: [%s, %s].\n", 
                round(x$lb[[i]], digits = 5),
                round(x$ub[[i]], digits = 5)))
  }
}

#' Summary of results from \code{invertci}
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'    
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{invertci}.
#' 
#' @export
#' 
summary.invertci <- function(x, ...){
  if (length(x$alpha) == 1){
    summary.invertci_single(x)
  } else {
    summary.invertci_multiple(x)
  }
}


#' Summary of results from \code{invertci} for a single alpha
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'    
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{invertci}.
#' 
#' @export
#' 
summary.invertci_single <- function(x, ...){
  #### Part 1: Display what has been the function
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  
  #### Part 2: Basic results
  cat(sprintf("Significance level: %s.\n", round(x$alpha, digits = 5)))  
  
  cat("\n") 
  
  #### Part 3: Messages in constructing the upper bound
  cat("=== Iterations in constructing upper bound:\n")
  cat(paste0("Iteration\t Lower bound \t Upper bound \t ",
             "Test point \t p-value\t Reject?\n"))
  for(j in 1:nrow(x$df_ub)){
    summary_bisection_print(x$df_ub, j)
  }
  cat(x$termination$ub)
  cat("\n\n")  
  
  #### Part 4: Messages in constructing the lower bound
  cat("=== Iterations in constructing lower bound:\n")
  cat(paste0("Iteration\t Lower bound \t Upper bound \t ",
             "Test point \t p-value\t Reject?\n"))
  for(j in 1:nrow(x$df_lb)){
    summary_bisection_print(x$df_lb, j)
  }
  cat(x$termination$lb)
  cat("\n\n")   
    
}

#' Summary of results from \code{invertci} with multiple alphas
#' 
#' @description This function uses the print method on the return list of the
#'    function \code{invertci}.
#'    
#' @param x Object returned from \code{invertci}.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{invertci}.
#' 
#' @export
#' 
summary.invertci_multiple <- function(x, ...){
  #### Part 1: Display what has been the function
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  
  #### Part 2: Basic results
  cat("Significance levels considered: ")
  for (i in 1:length(x$alpha)){
    cat(paste(round(unlist(x$alpha[i]), digits = 5), " "))
  }
  
  cat("\n\n") 
  
  for (i in 1:length(x$alpha)){
    cat(sprintf("< Confidence interval for alpha = %s >\n", x$alpha[i]))
    
    #### Part 3: Messages in constructing the upper bound
    cat("=== Iterations in constructing upper bound:\n")
    cat(paste0("Iteration\t Lower bound \t Upper bound \t ",
               "Test point \t p-value\t Reject?\n"))
    for(j in 1:nrow(x$df_ub[[i]])){
      summary_bisection_print(x$df_ub[[i]], j)
    }
    cat(x$termination[[i]]$ub)
    cat("\n\n")  
    
    #### Part 4: Messages in constructing the lower bound
    cat("=== Iterations in constructing lower bound:\n")
    cat(paste0("Iteration\t Lower bound \t Upper bound \t ",
               "Test point \t p-value\t Reject?\n"))
    for(j in 1:nrow(x$df_lb[[i]])){
      summary_bisection_print(x$df_lb[[i]], j)
    }
    cat(x$termination[[i]]$lb)
    cat("\n\n")   
  }
}



#' Print results in constructing bounds in bisection method
#' 
#' @description This function is used to display the message when constructing
#'    the bounds and used in \code{summary.invertci} to print the results in 
#'    each step of the bisection method.
#' 
#' @inheritParams bisec_print
#' 
#' @return Nothing is returned in this function.
#' 
#' @export
#' 
summary_bisection_print <- function(df_bis, i){
  
  ### Data cleaning
  print_iter1 = df_bis[i, 1]
  if (print_iter1 == "Left end pt." | print_iter1 == "Right end pt."){
    print_iter1 = paste(print_iter1, "\t")
  } else {
    print_iter1 = paste(print_iter1, "\t\t")
  }
  
  print_iter2 = df_bis[i, 2]
  if (print_iter2 != "NA"){
    print_iter2 = format(round(as.numeric(print_iter2), digits = 5), nsmall = 5)
  } else {
    print_iter2 = "NA\t"
  }
  
  print_iter3 = df_bis[i, 3]
  if (print_iter3 != "NA"){
    print_iter3 = format(round(as.numeric(print_iter3), digits = 5), nsmall = 5)
  } else {
    print_iter3 = "NA\t"
  }
  
  print_iter4 = df_bis[i, 4]
  
  print_iter5 = format(round(as.numeric(df_bis[i, 5]), digits = 5), nsmall = 5)
  
  print_iter6 = df_bis[i, 6]
  
  ### Print results
  cat(paste(print_iter1, 
            print_iter2,"\t",
            print_iter3,"\t",
            format(round(print_iter4, digits = 5), nsmall = 5),"\t",
            print_iter5,"\t",
            print_iter6,"\t\n"))
  
}

