## ========================================================================= ##
##
##  Example for the missing data problem
##
##  This is an example code of generating an `lpmodel` object and drawing a
##  DGP based on the missing data problem.
##
## ========================================================================= ##

#' Creates the \code{lpmodel} object for the missing data problem
#' 
#' @description This function creates the \code{lpmodel} object for the 
#'   missing data problem.
#'
#' @param J The number of distinct outcomes in Y.
#' @param info The string that indicates whether the full-information or the 
#'   two-moments approach is used.
#' @param data The data set that represents the missing data problem.
#' 
#' @return An \code{lpmodel} object.
#' 
missingdata_lpm <- function(J = 5, info = "mean", data = NULL) {
  # Initialize the parameters
  J1 <- J + 1

  # If data is not provided, draw the data here
  if (is.null(data)) {
    data <- missingdata_draw(J = J)
  }

  # Construct the components in the lpmodel object
  # The 'A.tgt' component are common in both approaches
  yp <- seq(0, 1, 1/J)
  A.tgt <- matrix(c(yp, yp), nrow = 1)
  A.shp <- rep(1, ncol(A.tgt))
  beta.shp <- 1
  if (info == "mean") {
    A.obs <- matrix(c(rep(0, J1), yp, rep(0, J1), rep(1, J1)),
                    nrow = 2,
                    byrow = TRUE)

    beta.obs <- function(data) {
      # Initialize beta
      beta <- matrix(c(0,0), nrow = 2)

      # Computes the two moments E[YD] and E[D]
      n <- nrow(data)
      beta[1, 1] <- sum(data$Y * data$D)/n
      beta[2, 1] <- sum(data$D)/n

      # Construct the variance matrix
      var <- matrix(0L, nrow = 2, ncol = 2)
      var[1, 1] <- sum((data$Y)^2 * (data$D)^2)/n - beta[1, 1]^2
      var[2, 2] <- sum((data$D)^2)/n - beta[2, 1]^2

      return(list(beta = beta,
                  var = var))
    }
  } else if (info == "full") {
    A.obs <- cbind(matrix(rep(0, J1 * J1), nrow = J1), diag(1, J1))

    beta.obs <- function(data) {
      # Count total number of rows
      n <- nrow(data)

      # Construct the beta vector and variance-covariance matrix
      beta <- NULL
      var <- matrix(0L, nrow = J1, ncol = J1)
      for (i in seq_along(yp)) {
        beta <- c(beta, sum((data$Y == yp[i]) * (data$D == 1))/n)
        var[i, i] <- sum(((data$Y == yp[i]) * (data$D == 1))^2)/n - beta[i]^2
      }
      beta <- as.matrix(beta)

      return(list(beta = beta,
                  var = var))
    }
  } else {
    stop(paste0("Please only use either 'mean' to represent the ",
                "two-moments approach or 'full' to represent the ",
                "full-information approach."))
  }

  # Construct the lpmodel object
  lpm <- lpmodel(A.obs = A.obs,
                 A.shp = A.shp,
                 A.tgt = A.tgt,
                 beta.obs = beta.obs,
                 beta.shp = beta.shp)
  
  return(lpm)
}

#' Simulate the missing data for the missing data problem
#' 
#' @description This function generates the simulated data for the missing
#'   data problem. There are two variables in this DGP (D, Y), where Y is the
#'   outcome variable and D is the variable that determines whether Y is 
#'   observed.
#'
#' @param J The number of distinct outcomes in Y.
#' @param n The number of observations.
#' @param seed The seed.
#' @param prob.obs The probability that Y is observed.
#' 
#' @return A data frame that contains the simulated data.
#' 
missingdata_draw <- function(J = 5, n = 500, seed = 1, prob.obs = .5) {
  set.seed(seed)

  # Initialize the data frame
  data <- data.frame(matrix(vector(), nrow = n, ncol = 2))
  colnames(data) <- c("D", "Y")

  # Draw indicator variable of whether Y is observed
  data$D <- e1071::rdiscrete(n,
                             probs = c(1 - prob.obs, prob.obs),
                             values = c(0, 1))

  # Draw outcome variable
  J1 <- J + 1
  data$Y <- e1071::rdiscrete(n, probs = rep(1/J1, J1), values = seq(0, 1, 1/J))

  return(data)
}
