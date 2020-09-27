#' General function to create a table of critical values
#'
#' @import scales
#'
#' @description This function generates a table of critical values for the
#'   testing procedures By default, the table includes the sample test
#'   statistic, and critical values at the 90\%, 95\% and 99\% level. The table
#'   can be retrieved from \code{cv.table} in the output of the functions.
#'
#' @param param The tuning parameter of the testing procedure.
#' @param param.name The name of the tuning parameter.
#' @param ts.sample.list A vector or sample test statistics.
#' @param ts.bs.list A vector or a list of bootstrap test statistics.
#' @param levels The confidence levels.
#' @param digits The number of digits in the test statistics.
#'
#' @return Returns a table of critical values.
#'   \item{cv.table}{A table of critical values.}
#'
#' @export
#'
construct.cv.table <- function(param, param.name, ts.sample.list,
                               ts.bs.list, levels = c(.99, .95, .9),
                               digits = 5) {
  # ---------------- #
  # Step 1: Compute the required the parameters and the data frame
  # ---------------- #
  n.param <- length(param)
  n.levels <- length(levels)
  cv.table <- data.frame(matrix(vector(),
                                nrow = n.levels + 1,
                                ncol = n.param + 1))

  # ---------------- #
  # Step 2: Construct the cv.table
  # ---------------- #
  # Name the columns
  colnames(cv.table) <- c(param.name, round(param, digits = digits))
  percents <- scales::label_percent(accuracy = 1)(levels)
  percents <- paste0("Bootstrap ", percents, " CV")
  cv.table[, 1] <- c("Sample", percents)

  # Construct the table
  for (i in 1:n.param) {
    # Construct the quantiles
    if (is.list(ts.bs.list)) {
      test.level.list <- quan.stat(ts.bs.list[[i]], levels)
    } else {
      test.level.list <- quan.stat(ts.bs.list, levels)
    }
    cv.table[1, i + 1] <- round(ts.sample.list[i], digits = digits)
    for (j in 1:n.levels)
      cv.table[j + 1, i + 1] <- round(test.level.list[j], digits = digits)
  }

  # ---------------- #
  # Step 3: Return the cv.table
  # ---------------- #
  return(cv.table)
}

#' Wrapper for the \code{\link[lpinfer]{cv.table}} function for the
#' \code{\link[lpinfer]{fsst}} procedure
#'
#' @description This is a wrapper of the \code{\link[lpinfer]{cv.table}}
#'   function in order to produce a table that consists of the critical values
#'   for the test statistics, cone component and the range component for the
#'   \code{\link[lpinfer]{fsst}} procedure.
#'
#' @param cone.n.list A list of sample cone test statistics.
#' @param range.n.list A list of sample range test statistics.
#' @param cone.bs.list A vector or a list of bootstrap cone test statistics.
#' @param range.bs.list A vector or a list of bootstrap range test statistics.
#' @param T.bs.list A vector or a list of bootstrap test statistics.
#' @inheritParams construct.cv.table
#'
#' @return Returns a table of critical values.
#'   \item{cv.table}{A table of critical values.}
#'
#' @export
#'
fsst.cv.table <- function(param, param.name, cone.n.list, range.n.list,
                          cone.bs.list, range.bs.list, T.bs.list,
                          levels = c(.99, .95, .9), digits = 5) {
  # ---------------- #
  # Step 1: Initialization
  # ---------------- #
  cv.range <- data.frame(matrix(vector(),
                                nrow = length(levels) + 1,
                                ncol = length(param) + 1))

  # ---------------- #
  # Step 2: Construct the CV table for the cone and range component
  # ---------------- #
  # Create table for the sample test statistics
  cv.test <- construct.cv.table(param, param.name,
                                rep(max(cone.n.list, range.n.list),
                                    length(T.bs.list)),
                                T.bs.list, levels, digits)
  
  # Create table for the cone component
  cv.cone <- construct.cv.table(param, param.name, cone.n.list,
                                cone.bs.list, levels, digits)

  # Create table for the range component
  cv.range.temp <- construct.cv.table(param[1], param.name, range.n.list,
                                      range.bs.list, levels, digits)
  colnames(cv.range) <- colnames(cv.cone)
  cv.range[, 1:2] <- cv.range.temp[, 1:2]

  # ---------------- #
  # Step 3: Construct the full CV table
  # ---------------- #
  # Construct the two parts
  cv.table.temp <- rbind(cv.test, cv.cone, cv.range)

  # Generate the column names
  n.levels <- length(levels)
  cv.table.name <- c("Test statistic", rep("", n.levels),
                     "Cone", rep("", n.levels),
                     "Range", rep("", n.levels))
  cv.table.name <- as.data.frame(cv.table.name)

  # Merge the two tables
  cv.table <- cbind(cv.table.name, cv.table.temp)
  colnames(cv.table)[1] <- ""

  # ---------------- #
  # Step 4: Return the cv.table
  # ---------------- #
  return(cv.table)
}

#' Function that computes the basic quantiles
#'
#' @description This function is used to evaluate the test statistics at
#'   different standard quantiles. By default, it evaluates the test
#'   statistics at the quantiles --- 90\%, 95\% and 99\%.
#'
#' @importFrom pracma ceil
#'
#' @param stat The test statistics.
#' @param quan The quantiles.
#'
#' @return Return the quantile of the test statistics in the order of the
#'   \code{quan} vector.
#'   \item{stat.quan}{The quantile of the test statistics in the order of the
#'   \code{quan} vector.}
#'
#' @export
#'
quan.stat <- function(stat, quan = c(.99, .95, .9)) {
  # ---------------- #
  # Step 1: Compute the basic parameters and initialize
  # ---------------- #
  n.stat <- length(stat)
  stat.order <- sort(stat)
  n.quan <- length(quan)
  stat.quan <- c()

  # ---------------- #
  # Step 2: Compute the quantiles via a for-loop
  # ---------------- #
  for (i in 1:n.quan) {
    temp <- stat.order[pracma::ceil(quan[i]*n.stat)]
    stat.quan <- c(stat.quan, temp)
  }

  return(stat.quan)
}
