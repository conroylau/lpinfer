#' Auxiliary function in the bootstrap replications
#'
#' @description This function is used in the bootstrap replications of the
#'   procedures. It returns the indices to be considered in the
#'   bootstrap replications and the list to be passed to the
#'   \code{\link[lpinfer]{lpmodel.anylist}} function.
#'
#' @inheritParams dkqs
#' @inheritParams subsample.bs
#' @param R.eval The number of bootstrap replications that has been evaluated.
#' @param R.succ The number of successful bootstrap replications.
#' @param maxR The maximum number of bootstrap replications.
#' @param any.list An object returned from the
#'   \code{\link[lpinfer]{lpmodel.anylist}} function, which indicates
#'   whether at least one of the components in \code{lpmodel} is a \code{list}.
#'
#' @return Returns the following three objects.
#'   \item{i0}{The starting index to be evaluated.}
#'   \item{i1}{The last index to be evaluated.}
#'   \item{bs.list}{The list of bootstrap replications in the \code{lpmodel}
#'     object passed from the user.}
#'
#' @export
#'
bs.assign <- function(R, R.eval, R.succ, maxR, any.list, lpmodel, data = NULL,
                      m = NULL, replace = TRUE) {
  # Evaluate the indices
  i0 <- min(R.eval + 1, maxR)
  i1 <- min(R - R.succ + R.eval, maxR)

  # Assign the list of objects or indices to be passed to 'future_lapply'
  if (isTRUE(any.list$list)) {
    bs.list <- any.list$consol[i0:i1]
  } else {
    bs.list <- as.list(i0:i1)
  }

  # Assign the subsample size
  if (is.null(m)) {
    m <- nrow(data)
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
          bs.temp <- lpmodel.beta.eval(data.bs, lpmodel[[i]], 1)
          if (is.null(bs.temp$omega)) {
            bs.list[[j]][[i]] <- bs.temp$beta.obs
          } else {
            bs.list[[j]][[i]] <- bs.temp
          }
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

#' Auxiliary function to return the indices for bootstrap
#' replications
#' 
#' @description This function is used to return the starting and ending
#'   indices for bootstrap replications.
#'   
#' @return Returns the following three objects.
#'   \item{i0}{The starting index to be evaluated.}
#'   \item{i1}{The last index to be evaluated.}
#'   
#' @inheritParams dkqs
#' @inheritParams bs.assign
#' 
bs.index <- function(R, R.eval, R.succ, maxR) {
  # Evaluate the indices
  i0 <- min(R.eval + 1, maxR)
  i1 <- min(R - R.succ + R.eval, maxR)
  return(list(i0 = i0,
              i1 = i1))
}

#' Auxiliary function for the post-bootstrap procedure
#'
#' @description This function is used after the bootstrap replications to
#'   extract the relevant and updated information.
#'
#' @inheritParams dkqs
#' @param test.return The list of information returned from the test function.
#' @param i0 The starting index to be evaluated.
#' @param i1 The last index to be evaluated.
#' @param R.eval The number of bootstrap replications that has been evaluated.
#' @param T.list The list of bootstrap test statistics.
#' @param beta.list The list of the bootstrap replications of the
#'   \code{beta.obs} component.
#' @param error.list The list of error messages.
#' @param param.list The list of parameters (if applicable).
#' @param error.param The list of parameters that lead to errors (if applicable).
#' @param ub.list The list of upper bounds (if applicable).
#' @param lb.list The list of lower bounds (if applicable).
#'
#' @return Returns the following three objects.
#'   \item{R.succ}{The number of successful bootstrap replications.}
#'   \item{R.eval}{The number of bootstrap replications that has been evaluated.}
#'   \item{T.list}{The list of test statistics.}
#'   \item{beta.list}{The list of the bootstrap replications of the
#'     \code{beta.obs} component.}
#'   \item{error.list}{The list of error messages (if applicable).}
#'   \item{error.param}{The list of parameters that lead to errors (if
#'     applicable).}
#'   \item{ub.list}{The list of upper messages (if applicable).}
#'   \item{lb.list}{The list of lower messages (if applicable).}
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
#' @param error.list The list of error messages.
#' @param df.error The consolidated data frame of error messages.
#'
#' @return Returns an updated list of error messages.
#'   \item{df.error}{An updated data frame of error messages.}
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

#' Coerces a \code{sparseMatrix} as a \code{matrix}
#' 
#' @description This function coerces a \code{sparseMatrix} as a 
#'   \code{matrix}. This function is used specifically in the 
#'   \code{\link[lpinfer]{gurobi.optim}} function to ensure that the 
#'   arguments are accepted by the \code{gurobi} function.
#' 
#' @param mat The matrix that is used in the \code{model} component.
#' 
#' @return Returns a matrix in the \code{matrix} form.
#' 
#' @export
#' 
smatrixconvert <- function(mat) {
  if (is(mat, "sparseMatrix") | is(mat, "Matrix")) {
    return(as.matrix(mat))
  } else {
    return(mat)
  }
}

