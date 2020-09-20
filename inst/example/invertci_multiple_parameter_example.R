## ========================================================================= ##
##
##  Example for the invertci function for multiple tuning parameters
##
##  This followings illustrate how the invertci function can construct
##  confidence intervals when there are multiple tuning paramters in the
##  testing procedure.
##
## ========================================================================= ##
rm(list = ls())

# ---------------- #
# Part 1: Load required packages
# ---------------- #
library(lpinfer)

# ---------------- #
# Part 2: Create a procedure that accepts two tuning parameters
# ---------------- #
eg_procedure <- function(beta.tgt, beta, gamma) {
  # Data frame that contains the p-values
  df <- data.frame(matrix(vector(),
                          nrow = length(beta) * length(gamma),
                          ncol = 3))
  colnames(df) <- c("beta", "gamma", "p-value")

  # Assign the p-values
  k <- 1
  for (i in seq_along(beta)) {
    for (j in seq_along(gamma)) {
      df[k, 1] <- beta[i]
      df[k, 2] <- gamma[j]
      df[k, 3] <- max(abs(beta.tgt - beta * 5 + gamma * 2), 1)
      k <- k + 1
    }
  }
  return(list(pval = df,
              logical.lb = 0,
              logical.ub = 1))
}

# ---------------- #
# Part 3: Construct confidence intervals when the two tuning paramters are
# vectors (i.e. multiple values for beta and gamma). There are also multiple
# significance levels.
# ---------------- #
set.seed(1)
r <- invertci(f = eg_procedure,
              farg = list(beta = c(.1, .3),
                          gamma = c(.1, .5)),
              init.lb = c(0, .4),
              init.ub = c(.6, 1),
              max.iter = 5,
              alpha = c(.9, .95),
              progress = FALSE)
print(r)
