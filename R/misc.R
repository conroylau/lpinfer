#' Auxiliary function in the bootstrap replications
#'
#' @description This function is used in the bootstrap replications of the
#'   tests. It returns the indices to be considered in the
#'   bootstrap replications and the list to be passed to the
#'   \code{lpmodel.anylist} function.
#'
#' @inheritParams dkqs
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
bs.assign <- function(R, R.eval, R.succ, maxR, any.list) {
  # Evaluate the indices
  i0 <- min(R.eval + 1, maxR)
  i1 <- min(R - R.succ + R.eval, maxR)

  # Assign the list of objects or indices to be passed to 'future_lapply'
  if (isTRUE(any.list$list)) {
    bs.list <- any.list$consol[i0:i1]
  } else {
    bs.list <- i0:i1
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
#'
#' @export
#'
post.bs <- function(test.return, i0, i1, R.eval, T.list = NULL,
                    beta.list = NULL, error.list = NULL, list.param = NULL,
                    error.param = NULL) {
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

  # Update the number of successful bootstrap replications
  R.succ <- length(beta.list)
  R.eval <- (i1 - i0 + 1) + R.eval

  return(list(R.succ = R.succ,
              R.eval = R.eval,
              T.list = T.list,
              beta.list = beta.list,
              error.list = error.list,
              error.param = error.param))
}

#' Gets the current seed with RNG kind "L'Ecuyer-CMRG"
#'
#' @description This function gets the current seed with RNG kind
#'   "L'Ecuyer-CMRG".
#'
#' @return Returns a \code{seed}.
#'   \item{seed}{Seed with RNG kind "L'Ecuyer-CMRG".}
#'
lpinfer.seed <- function() {
  RNGkind("L'Ecuyer-CMRG")
  seed <- .Random.seed

  return(seed)
}
