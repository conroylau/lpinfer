#' Auxiliary function in the bootstrap replications
#'
#' @description This function is used in the bootstrap replications of the
#'   tests. It returns the indices to be considered in the
#'   bootstrap replications and the list to be passed to the
#'   \code{lpmodel.anylist} function.
#'
#' @inheritParams dkqs
#' @inheritParams subsample.bs
#' @param R.eval Number of bootstrap replications that has been evaluated/
#' @param R.succ Number of successful bootstrap replications.
#' @param maxR Maximum number of bootstrap replications.
#' @param any.list An object returned from the \code{lpmodel.anylist}
#'   function, which indicates whether at least one of the components in
#'   \code{lpmodel} is a \code{list}.
#'
#' @return Returns the following three objects.
#'   \item{i0}{Starting index to be evaluated.}
#'   \item{i1}{Last index to be evaluated.}
#'   \item{bs.list}{List of bootstrap replications in the \code{lpmodel}
#'     object passed from the user.}
#'
#' @export
#'
bs.assign <- function(R, R.eval, R.succ, maxR, any.list, lpmodel, data = NULL,
                      m, replace = TRUE) {
  # Evaluate the indices
  i0 <- min(R.eval + 1, maxR)
  i1 <- min(R - R.succ + R.eval, maxR)

  # Assign the list of objects or indices to be passed to 'future_lapply'
  if (isTRUE(any.list$list)) {
    bs.list <- any.list$consol[i0:i1]
  } else {
    bs.list <- as.list(i0:i1)
  }
  
  # Temporary solution to the missing globals problem
  names <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")
  for (i in names) {
    if (is.function(lpmodel[[i]])) {
      for (j in seq_along(bs.list)) {
        if (!is.list(bs.list[[j]])) {
          bs.list[[j]] <- list()
        }
        data.bs <- as.data.frame(data[sample(1:nrow(data), m, replace),])
        if (i == "beta.obs") {
          bs.list[[j]][[i]] <- lpmodel.beta.eval(data.bs, lpmodel[[i]], 1)
        } else {
          bs.list[[j]][[i]] <- lpmodel.eval(data.bs, lpmodel[[i]], 1)
        }
      }
    }
  }

  return(list(i0 = i0,
              i1 = i1,
              bs.list = bs.list))
}

#' Auxiliary function for the post-bootstrap procedure
#'
#' @description This function is used after the bootstrap replications to
#'   extract the relevant and updated information.
#'
#' @inheritParams dkqs
#' @param test.return List of information returned from the test function.
#' @param i0 Starting index to be evaluated.
#' @param i1 Last index to be evaluated.
#' @param R.eval Number of bootstrap replications that has been evaluated.
#' @param T.list List of bootstrap test statistics.
#' @param beta.list List of the bootstrap replications of the \code{beta.obs}
#'     component.
#' @param error.list List of error messages.
#' @param param.list List of parameters (if applicable).
#' @param error.param List of parameters that lead to errors (if applicable).
#' @param ub.list List of upper bounds (if applicable).
#' @param lb.list List of lower bounds (if applicable).
#'
#' @return Returns the following three objects.
#'   \item{R.succ}{Number of successful bootstrap replications.}
#'   \item{R.eval}{Number of bootstrap replications that has been evaluated.}
#'   \item{T.list}{List of test statistics.}
#'   \item{beta.list}{List of the bootstrap replications of the
#'     \code{beta.obs} component.}
#'   \item{error.list}{List of error messages (if applicable).}
#'   \item{error.param}{List of parameters that lead to errors (if
#'     applicable).}
#'   \item{ub.list}{List of upper messages (if applicable).}
#'   \item{lb.list}{List of lower messages (if applicable).}
#'
#' @export
#'
post.bs <- function(test.return, i0, i1, R.eval, T.list = NULL,
                    beta.list = NULL, error.list = NULL, list.param = NULL,
                    error.param = NULL, ub.list = NULL, lb.list = NULL) {
  # Update the lists
  if (!is.null(T.list)) {
    T.list <- c(unlist(T.list), unname(unlist(sapply(test.return, "[", "Ts"))))
  }
  if (!is.null(beta.list)) {
    beta.list <- c(beta.list, sapply(test.return, "[", "beta"))
    beta.list <- beta.list[lengths(beta.list) != 0]
  }
  if (!is.null(error.list)) {
    error.list <- c(error.list, sapply(test.return, "[", "msg"))
  }
  if (!is.null(error.param)) {
    error.param <- c(error.param, sapply(test.return, "[", "param"))
  }
  if (!is.null(ub.list)) {
    ub.list <- c(ub.list, sapply(test.return, "[", "ub"))
  }
  if (!is.null(lb.list)) {
    lb.list <- c(lb.list, sapply(test.return, "[", "lb"))
  }

  # Update the number of successful bootstrap replications
  R.succ <- length(beta.list)
  R.eval <- (i1 - i0 + 1) + R.eval

  return(list(R.succ = R.succ,
              R.eval = R.eval,
              T.list = T.list,
              beta.list = beta.list,
              error.list = error.list,
              error.param = error.param,
              ub.list = ub.list,
              lb.list = lb.list))
}

#' Matches the id of the error messages
#'
#' @description This function is used after the bootstrap replications. This
#'   function is used to match the id of the error messages, given the list
#'   of error messages returned from the \code{future_lapply} function.
#'
#' @param error.list List of error messages.
#' @param df.error Consolidated data frame of error messages.
#'
#' @return Returns an updated list of error messages.
#'   \item{df.error}{Updated data frame of error messages.}
#'
#' @export
#'
error.id.match <- function(error.list, df.error) {
  # Remove 'NULL' from the list before passing to future_lapply
  df.error.nonnull <- error.list[-which(sapply(error.list, is.null))]

  # Match the id
  for (i in seq_along(df.error.nonnull)) {
    df.error$id[i] <- match(df.error.nonnull[i], error.list)
    error.list[df.error$id[i]] <- NA
  }

  return(df.error)
}
