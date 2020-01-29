#' Constructs confidence interval
#' 
#' @description This function constructs the confidence interval using the
#'    bisection method. 
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
#' @param noisy The boolean variable for whether the result messages should
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
#' 
#' @export 
#' 
qpci <- function(f, farg, alpha = .05, lb0 = NULL, lb1 = NULL, ub0 = NULL, 
                 ub1 = NULL, tol = .0001, max_iter = 20, df_ci = NULL,
                 noisy = TRUE){
  
  #### Step 1: Check and update the dependencies
  check_return = qpci_check(f, farg, alpha, lb0, lb1, ub0, ub1, tol, max_iter, 
                            df_ci, noisy)
  # Updates the dependencies
  df_ci = check_return$df_ci
  lb0 = check_return$lb0
  lb1 = check_return$lb1
  ub0 = check_return$ub0
  ub1 = check_return$ub1
  
  #### Step 2: Return confidence interval and data frame
  ### Compute upper bound of confidence interval
  if (noisy == TRUE){
    cat(">>> Computing upper bound of confidence interval")
  }
  up_return = ci_bisection(f, farg, alpha, ub1, ub0, tol, max_iter, df_ci, 
                           noisy, 1)
  # Update data frame
  df_ci = up_return$df_ci
  ### Compute lower bound of confidence interval
  if (noisy == TRUE){
    cat("\n>>> Computing lower bound of confidence interval")
  }
  down_return = ci_bisection(f, farg, alpha, lb0, lb1, tol, max_iter, df_ci, 
                             noisy, -1)
  # Update data frame
  df_ci = down_return$df_ci

  #### Step 3: Print confidence interval and return results
  if (noisy == TRUE){
    cat("-----------------------------------\n")
    cat(paste("Total number of iterations: ", 
              down_return$iter + up_return$iter, ".\n", sep = ""))
    cat(paste("Tolerance level: ", tol, ".\n", sep = ""))
    cat(paste("Confidence interval: [", round(down_return$pt, digits = 5), 
              ", ", round(up_return$pt, digits = 5), "].\n", sep = "")) 
  }
  
  invisible(list(up = up_return$pt,
                 down = down_return$pt,
                 df_ci = df_ci))
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
#' @inheritParams qpci
#' 
#' @return Return the solution of the bisection method and the updated 
#'    data frame.
#'    \item{soln}{Solution to the bisection method.}
#'    \item{df_ci}{data frame that consists of the points and the corresponding
#'       \eqn{p}-values that have been tested in constructing the confidence 
#'       intervals.}   
#' 
#' @export 
#'
ci_bisection <- function(f, farg, alpha, b0, b1, tol, max_iter, df_ci, noisy,
                         type){
  
  #### Step 1: Evaluate the end-points and the mid-point of b0 and b1
  # Divide alpha by 2
  alpha_2sided = alpha/2
  # Left end-point a
  a = b0
  fb0_return = bisec_eval(f, farg, a, df_ci)
  fb0 = fb0_return$pval
  df_ci = fb0_return$df_ci
  # Right end-point b
  b = b1
  fb1_return = bisec_eval(f, farg, b, df_ci)
  fb1 = fb1_return$pval
  df_ci = fb1_return$df_ci

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
    if ((b-a) < tol){
      if (noisy == TRUE){
        cat(paste("\n       Length of interval is below tolerance level. ",
                  "Bisection method is completed.\n", sep = "")) 
      }
      break
    }
    # Display a dot to represent an iteration is completed
    if (noisy == TRUE){
      cat(".") 
    }
    # Update interval based on whether the left section or the section of the 
    # interval is chosen
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
  if (noisy == TRUE & i == max_iter){
    cat(paste("\n       Reached the maximum number of iterations.\n", 
              sep = "")) 
  }
  
  #### Step 3: Return results
  # Return the results
  invisible(list(pt = c,
                 iter = i,
                 df_ci = df_ci))
}

#' Evaluation of test statistic and check if the point has been evaluated
#' 
#' @description This function checks if the \eqn{p}-value for the point 
#'    considered has already been evaluated in previous iterations or provided 
#'    by the user. The function will compute the \eqn{p}-value if it has been 
#'    evaluated. Otherwise, it will use the previous data.
#' 
#' @inheritParams qpci
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
#' @inheritParams qpci
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

#' Checks and updates the input of the function \code{qpci}
#' 
#' @description This function checks and updates the input from the user for 
#'    the function \code{qpci}. If there is any invalid input, this function 
#'    will terminate the procedure and generate appropriate error messages.
#'    
#' @inheritParams qpci
#' 
#' @return Returns the updated value of the parameters back to the function 
#'    \code{qpci} in the correct format.
#'    \item{df_ci}{Data frame that stores the points that has been evaluated
#'       and the corresponding \eqn{p}-values.}
#'    \item{lb0}{Logical lower bound for the confidence interval.}
#'    \item{lb1}{Maximum possible lower bound for the confidence interval.}
#'    \item{ub0}{Logical upper bound for the confidence interval.}
#'    \item{ub1}{Minimum possible upper bound for the confidence interval.}
#' 
#' @export
#' 
qpci_check <- function(f, farg, alpha, lb0, lb1, ub0, ub1, tol, max_iter, 
                       df_ci, noisy){

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
  if ((is.numeric(alpha) == TRUE & length(alpha) == 1 & alpha >= 0 & 
       alpha <= 1) == FALSE){
    stop("The significance level ('alpha') has to be in the interval [0,1].", 
         call. = FALSE)
  }    
  
  #### Part 4: Check lb0 and ub0
  if (is.null(lb0) | is.null(ub0)){
    ### Part A: If lb0 and ub0 are null, check whether if the function is dkqs
    ### If yes, compute lb0 and ub0. Otherwise, return terminate the function.
    if (as.character(substitute(f)) == "dkqs_cone"){
      logicalb_return = dkqs_cone_logicalb(f, farg)
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
      if (colnames(df_ci) != c("point", "value")){
        stop("The column names of the data frame 'df_ci' have to be 'point' 
             and 'value'.", call. = FALSE)
      }
      # Check if the values are numeric
      if (is.numeric(df_ci) == FALSE){
        print(is.numeric(df_ci))
        print(df_ci)
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

  #### Part 10: Check noisy
  if (!(noisy == TRUE | noisy == FALSE)){
    stop("The argument 'noisy' has to be boolean.")
  }
  
  ### Step 11: Return the upated information
  return(list(df_ci = df_ci,
              lb0 = lb0,
              lb1 = lb1,
              ub0 = ub0,
              ub1 = ub1))
}

#' Compute the logical upper and lower bounds for dkqs_cone
#' 
#' @description This function computes the logical upper and lower bounds for
#'    the test \code{dkqs_cone}.
#'    
#' @inheritParams qpci
#' 
#' @return Returns the logical upper and lower bounds for dkqs_cone.
#'    \item{lb0}{Logical lower bound for \code{dkqs_cone}.}
#'    \item{ub0}{Logical upper bound for \code{dkqs_cone}.}
#' 
#' @export
#' 
dkqs_cone_logicalb <- function(f, farg){
  
  #### Step 1: Assign a value to beta_tgt for returning the results
  farg$beta_tgt = 0
  
  #### Step 2: Run dkqs_cone to obtain the logical bounds
  dkqs_return = do.call(f, farg)
  lb0 = dkqs_return$lb0
  ub0 = dkqs_return$ub0

  #### Step 3: Return results
  invisible(list(lb0 = lb0, 
                 ub0 = ub0))
}


#' Wrapper function for \code{qpci}
#' 
#' @description This function serves a wrapper function for computing 
#'    confidence interval with many significance levels.
#' 
#' @param alphas The list of significance levels to be used in constructing
#'    the confidence intervals.
#' @param noisy_one The boolean variable for whether the result messages should
#'    be displayed in running the function \code{qpci}. If it is set as 
#'    \code{TRUE}, the messages are displayed throughout the procedure. 
#'    Otherwise, the messages will not be displayed.
#' @param noisy_many The boolean variable for whether the result messages 
#'    should be displayed in running the function \qpci{many_qpci}, i.e. the 
#'    current function. If it is set as \code{TRUE}, the messages are displayed 
#'    throughout the procedure. Otherwise, the essages will not be displayed.
#' @inheritParams qpci 
#' 
#' @return Returns a list of confidence intervals.
#' 
#' @export
#' 
many_qpci <- function(f, farg, alphas = c(.05), lb0 = NULL, lb1 = NULL, 
                      ub0 = NULL, ub1 = NULL, tol = NULL, max_iter = 10, 
                      df_ci = NULL, noisy_one = FALSE, noisy_many = TRUE){
  
  #### Step 1: Check input of alphas
  if (!(is.numeric(alphas) == TRUE)) {
    stop("The argument 'alphas' must be in numeric format.", call. = FALSE)
  }  
  
  #### Step 2: Initialization and hard-coded information
  # Count number of alphas
  alpha_n = length(alphas)
  # Empty data frame to collect the confidence intervals
  df_many_ci = data.frame(matrix(vector(), 0, 3,
                                 dimnames=list(c(), c("alpha", "down", "up"))),
                          stringsAsFactors=F)
  
  #### Step 3: For-loop to compute the confidence intervals
  for (i in 1:alpha_n){
    # Call qpci to compute the confidence interval for a particular value of 
    # alpha
    qpci_return = qpci(f, farg, alphas[i], lb0, lb1, ub0, ub1, tol, max_iter, 
                       df_ci, noisy_one)
    # Obtain result
    df_many_ci[i, "alpha"] = alphas[i]
    df_many_ci[i, "up"] = qpci_return$up
    df_many_ci[i, "down"] = qpci_return$down
    # Update data frame
    df_ci = qpci_return$df_ci
    # Print result if noisy_many == TRUE
    if (noisy_many == TRUE){
      cat(paste("Confidence interval for significance level ", alphas[i], 
                ": [", round(qpci_return$down, digits = 5), ", ", 
                round(qpci_return$up, digits = 5),  "].\n", sep = ""))
    }
  }
  
  #### Step 4: Return result
  invisible(list(df_many_ci = df_many_ci))
}

