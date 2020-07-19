#' Combine lists in parallel loops
#'
#' @description This function is used to combine the list of objects in the
#'   parallelized bootstrap replications.
#'
#' @param x List of objects returned.
#' @param ... Additional arguments.
#'
#' @details This is used in the parallelized bootstrap replications.
#' @export
#'
para.comb <- function(x, ...) {
  lapply(seq_along(x), function(i) c(x[[i]],
                                     lapply(list(...),
                                            function(y) y[[i]])))
}
