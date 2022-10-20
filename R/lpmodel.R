#' Evaluates an object inside \code{lpmodel}
#'
#' @importFrom methods is
#'
#' @description This function returns the matrix or vector depending on the
#'   class of the variable in the \code{lpmodel} object. In the design of
#'   the \code{lpinfer} module, objects in \code{lpmodel} can have three
#'   different types of classes:
#'   \itemize{
#'      \item{\code{function} --- If the object is a function, then the
#'         object that is evaluated by the function with the data will be
#'         returned.}
#'      \item{\code{list} --- If the object is a list, then the \eqn{i}-th
#'         element will be returned.}
#'      \item{\code{matrix} or \code{numeric} --- If the object is a matrix
#'         or in a vector, then it will be directly returned.}
#'   }
#'
#' @param data A data frame.
#' @param obj An object in lpmodel.
#' @param i An index.
#'
#' @return Returns an object at iteration \code{i}.
#'
#' @export
#'
lpmodel.eval <- function(data, obj, i) {
  if (inherits(obj, "function")) {
    obj.eval <- obj(data)
  } else if (inherits(obj, "list")) {
    obj.eval <- obj[[i]]
  } else if (!is.matrix(obj) & !is.data.frame(obj) &
             !methods::is(obj, "sparseMatrix")) {
    obj.eval <- matrix(obj, nrow = 1)
  } else if (inherits(obj, "data.frame")) {
    obj.eval <- as.matrix(obj)
  } else {
    obj.eval <- obj
  }

  return(obj.eval)
}

#' Evaluates the point estimate and asymptotic variance of \code{beta.obs}
#'
#' @description This function returns the matrix or vector depending on the
#'   class of the variable in the \code{lpmodel} object. In the design of
#'   the \code{lpinfer} module, objects in \code{lpmodel} can have three
#'   different types of classes:
#'   \itemize{
#'      \item{\code{function} --- If the object is a function, then the
#'         object that is evaluated by the function with the data will be
#'         returned.}
#'      \item{\code{list} --- If the object is a list, then the \eqn{i}-th
#'         element will be returned.}
#'      \item{\code{matrix} or \code{numeric} --- If the object is a matrix
#'         or a vector, then it will be directly returned.}
#'   }
#'
#' @param data A data frame.
#' @param obj An object in lpmodel.
#' @param i An index.
#'
#' @return Returns the point estimate and the asymptotic variance of the
#'    \code{beta.obs} object.
#'    \item{beta.obs}{The point estimate of
#'      \eqn{\widehat{\bm{\beta}}_{\rm obs}}.}
#'    \item{omega}{The estimator of the asymptotic variance for
#'    \eqn{\widehat{\bm{\beta}}_{\rm obs}}.}
#'
#' @export
#'
lpmodel.beta.eval <- function(data, obj, i) {
  if (inherits(obj, "function")) {
    beta.return <- obj(data)
    if (inherits(beta.return, "list")) {
      if (length(beta.return) == 2) {
        if (is.null(nrow(beta.return[[1]]))) {
          beta.obs.hat <- beta.return[[1]]
          omega.hat <- beta.return[[2]]
        } else if (nrow(beta.return[[1]] == 1 | ncol(beta.return[[1]]) == 1)) {
          beta.obs.hat <- beta.return[[1]]
          omega.hat <- beta.return[[2]]
        } else {
          beta.obs.hat <- beta.return[[2]]
          omega.hat <- beta.return[[1]]
        }
      } else if (length(beta.return) == 1) {
        beta.obs.hat <- beta.return[[1]]
        omega.hat <- NULL
      } else {
        stop(paste0("If the output when the data is passed to the function ",
                    "is a list, it can have at most two objects in the ",
                    "list."))
      }
    } else {
      beta.obs.hat <- beta.return
      omega.hat <- NULL
    }
  } else if (inherits(obj, "list")) {
    if (inherits(obj[[i]], "list")) {
      if (length(obj[[i]]) != 2) {
        stop(paste0("When the first object of the list 'beta.obs' is a ",
                    "list, it needs to have two objects"))
      } else if (length(obj[[i]]) == 1) {
        beta.obs.hat <- obj[[i]][[1]]
        omega.hat <- NULL
      } else {
        if (is.null(nrow(obj[[i]][[1]]))) {
          beta.obs.hat <- obj[[i]][[1]]
          omega.hat <- obj[[i]][[2]]
        } else if (nrow(obj[[i]][[1]]) == 1 | ncol(obj[[i]][[1]]) == 1) {
          beta.obs.hat <- obj[[i]][[1]]
          omega.hat <- obj[[i]][[2]]
        } else {
          beta.obs.hat <- obj[[i]][[2]]
          omega.hat <- obj[[i]][[1]]
        }
      }
    } else {
      if (length(obj) == 2) {
        if (is.null(nrow(obj[[1]]))) {
          beta.obs.hat <- obj[[1]]
          omega.hat <- obj[[2]]
        } else if (nrow(obj[[1]]) == 1 | ncol(obj[[1]]) == 1) {
          beta.obs.hat <- obj[[1]]
          omega.hat <- obj[[2]]
        } else {
          beta.obs.hat <- obj[[2]]
          omega.hat <- obj[[1]]
        }
      } else {
        beta.obs.hat <- obj[[i]]
        omega.hat <- NULL
      }
    }
  } else if (inherits(obj, "numeric") | inherits(obj, "matrix")) {
    beta.obs.hat <- obj
    omega.hat <- NULL
  }

  if (inherits(beta.obs.hat, "data.frame")) {
    beta.obs.hat <- as.matrix(beta.obs.hat)
  } else if (inherits(beta.obs.hat, "numeric")) {
    beta.obs.hat <- matrix(beta.obs.hat, ncol = 1)
  }

  return(list(beta.obs = beta.obs.hat,
              omega = omega.hat))
}

#' Defines a \code{lpmodel} object
#'
#' @description This function defines the objects required in the
#'    \code{lpinfer} package in the \code{lpmodel} class.
#'
#' @param A.obs A matrix, list or function.
#' @param A.shp A matrix, list or function.
#' @param A.tgt A matrix, list or function.
#' @param beta.obs A vector, list or function.
#' @param beta.shp A vector, list or function.
#'
#' @return Returns a \code{lpmodel} object.
#'
#' @export
#'
lpmodel <- function(A.obs = NULL, A.shp = NULL, A.tgt = NULL, beta.obs = NULL,
                    beta.shp = NULL) {
  # ---------------- #
  # Step 1: Define the lpmodel objects
  # ---------------- #
  lpm <- list()
  lpm$A.obs <- A.obs
  lpm$A.shp <- A.shp
  lpm$A.tgt <- A.tgt
  lpm$beta.obs <- beta.obs
  lpm$beta.shp <- beta.shp

  # ---------------- #
  # Step 2: Define the class of the model
  # ---------------- #
  class(lpm) <- "lpmodel"

  return(lpm)
}

#' Define an \code{lpmodel.natural} form object
#'
#' @description This function defines the objects required in the
#'    \code{lpinfer} package in the \code{lpmodel.natural} class that allows
#'    both equality and inequality constraints.
#'
#' @param sense.shp The sense vector for the shape constraints.
#' @param x.lb The lower bound for the \eqn{\bm{x}} variable.
#' @param x.ub The upper bound for the \eqn{\bm{x}} variable.
#' @inheritParams lpmodel
#'
#' @return Returns a list of \code{lpmodel} objects in the \code{lpmodel} class.
#'
#' @export
#'
lpmodel.natural <- function(A.obs = NULL, A.shp = NULL, A.tgt = NULL,
                            beta.obs = NULL, beta.shp = NULL,
                            sense.shp = NULL, x.lb = NULL, x.ub = NULL) {
  # ---------------- #
  # Step 1: Define the lpmodel objects
  # ---------------- #
  lpm.natural <- list()
  lpm.natural$A.obs <- A.obs
  lpm.natural$A.shp <- A.shp
  lpm.natural$A.tgt <- A.tgt
  lpm.natural$beta.obs <- beta.obs
  lpm.natural$beta.shp <- beta.shp
  lpm.natural$sense.shp <- sense.shp
  lpm.natural$x.lb <- x.lb
  lpm.natural$x.ub <- x.ub

  # ---------------- #
  # Step 2: Define the class of the model
  # ---------------- #
  class(lpm.natural) <- "lpmodel.natural"

  return(lpm.natural)
}

#' Print the \code{lpmodel} or \code{lpmodel.natural} object
#'
#' @importFrom methods is
#'
#' @description This function prints the details of the components that are
#'    contained in the \code{lpmodel} or \code{lpmodel.natural} object.
#'
#' @inheritParams dkqs
#' @param x The \code{lpmodel} object or \code{lpmodel.natural} object.
#' @param lpm.string The string that contains the name of the variables
#'    available in the \code{lpmodel} object.
#' @param ... Some additional arguments passed to the \code{lpm.print}
#'    function.
#'
#' @return Nothing is returned.
#'
#' @export
#'
lpm.print <- function(x, lpm.string, data = NULL, ...) {
  # List of variables
  lpmodel.string <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")
  lpmodel.ind <- NULL
  for (i in 1:length(lpm.string)) {
    if (!is.null(x[[lpm.string[i]]])) {
      lpmodel.ind <- c(lpmodel.ind, i)
    }
  }

  if (length(lpmodel.ind) == 0) {
    cat("'lpmodel' object does not contain the required objects.")
  } else {
    cat("Object     Class \t\tDimension \tLength \n")
    for (i in 1:length(lpm.string)) {
      if (i %in% lpmodel.ind) {
        obj <- x[[lpm.string[i]]]

        # Check class of object
        class.tmp <- class(obj)
        # Concatenate the class names if length(class.tmp) > 1
        # E.g., if class = c("matrix", "array"), then it is "matrix, array"
        class.tmp <- paste(class.tmp, collapse = ", ")

        if (inherits(obj, "list")) {
          class.tmp <- "list  "
          length.tmp <- length(obj)
          dimension.str <- dim(as.matrix(obj[[1]]))
          dimension.tmp <- paste0(dimension.str[1], "x", dimension.str[2])
        } else if (inherits(obj, "function")) {
          # If data is not passed, print "N/A" for dimensions. Otherwise,
          # compute the output for the output object
          if (is.null(data)) {
            length.tmp <- "N/A"
            dimension.tmp <- "N/A"
          } else {
            tmp.obj <- obj(data)
            if (inherits(tmp.obj, "list")) {
              length.tmp <- length(tmp.obj)
            } else {
              length.tmp <- 1
            }
            dimension.str <- dim(as.matrix(tmp.obj))
            dimension.tmp <- paste0(dimension.str[1], "x", dimension.str[2])
          }
        } else if (inherits(obj, "data.frame") |
                   inherits(obj, "matrix") |
                   inherits(obj, "numeric") |
                   methods::is(obj, "sparseMatrix")) {
          dim.obj <- dim(obj)
          if (is.null(dim.obj)) {
            dimension.tmp <- paste0("1x", length(obj))
          } else {
            dimension.tmp <- paste0(dim.obj[1], "x", dim.obj[2])
          }
          length.tmp <- 1
        } else {
          length.tmp <- length(obj)
          dimension.tmp <- "  "
        }
        if (nchar(class.tmp) > 12) {
          st <- "%s\t%s\t\t%s\n"
        } else {
          st <- "%s\t \t%s\t\t%s\n"
        }
        cat(sprintf(paste0("%s",
                           paste(rep(" ", 11 - nchar(lpm.string[i])),
                                 collapse = ""),
                           st),
                           #"%s\t \t%s\t\t%s\n"),
                    lpm.string[i], class.tmp, dimension.tmp, length.tmp))
      } else {
        cat(sprintf(paste0("%s",
                           paste(rep(" ", 11 - nchar(lpm.string[i])),
                                 collapse = ""),
                           "-empty-\t\t-empty-\t\t-empty-\n"),
                    lpm.string[i]))
      }
    }
  }
}

#' Print the \code{lpmodel} object
#'
#' @description This function is a wrapper of the
#'   \code{\link[lpinfer]{lpm.print}} function and prints the details of the
#'   components in the \code{lpmodel} object.
#'
#' @param x An \code{lpmodel} object.
#' @inheritParams dkqs
#' @inheritParams lpm.print
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.lpmodel <- function(x, data = NULL, ...) {
  # List of variables
  lpmodel.string <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")
  lpm.print(x, lpmodel.string, data)
}

#' Summary of the \code{lpmodel} object
#'
#' @description This function is a wrapper of the \code{print.lpmodel}
#'    function and prints the same information for the object \code{lpmodel}.
#'
#' @param x The \code{lpmodel} object.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned
#'
#' @export
#'
summary.lpmodel <- function(x, ...) {
  print(x)
}

#' Print the \code{lpmodel.natural} object
#'
#' @description This function is a wrapper of the
#'   \code{\link[lpinfer]{lpm.print}} function and prints the details of the
#'   components in the \code{lpmodel.natural} object.
#'
#' @param x An \code{lpmodel.natural} object.
#' @inheritParams dkqs
#' @inheritParams lpm.print
#'
#' @return Nothing is returned.
#'
#' @export
#'
print.lpmodel.natural <- function(x, data = NULL, ...) {
  # List of variables
  lpmodel.natural.string <- c("A.obs", "A.shp", "A.tgt", "beta.obs",
                              "beta.shp", "sense.shp", "x.lb", "x.ub")
  lpm.print(x, lpmodel.natural.string, data)
}

#' Summary of the \code{lpmodel.natural} object
#'
#' @description This function is a wrapper of the \code{print.lpmodel.natural}
#'    function and prints the same information as the object
#'    \code{lpmodel.natural}.
#'
#' @param x The \code{lpmodel.natural} object.
#' @param ... Additional arguments.
#'
#' @return Nothing is returned.
#'
#' @export
#'
summary.lpmodel.natural <- function(x, ...) {
  print(x)
}

#' Check if there is any list in the \code{lpmodel} object
#'
#' @description This function checks if there is any components in the
#'   \code{lpmodel} object with class \code{list}.
#'
#' @param lpmodel A \code{lpmodel} object.
#'
#' @return Returns the following objects:
#'   \item{list}{A boolean variable that indicates whether there is any object
#'      with the class \code{list} in the \code{lpmodel} object.}
#'   \item{name}{The names of the components with class \code{list} in the
#'     \code{lpmodel} object.}
#'   \item{consol}{A consolidated and updated \code{lpmodel} object that lists
#'     the \code{lpmodel} by observation instead of by component (if
#'     applicable).}
#'   \item{len}{The length of the component in the \code{lpmodel} object that
#'     is a list (if applicable).}
#'
#' @export
#'
lpmodel.anylist <- function(lpmodel) {
  # Initialize the variables to be returned
  any.list <- FALSE
  name.list <- NULL
  len <- NULL

  # Check the objects in 'lpmodel' one-by-one
  names <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")
  for (i in names) {
    if (inherits(lpmodel[[i]], "list")) {
      any.list <- TRUE
      if (is.null(len)) {
        len <- length(lpmodel[[i]])
      } else {
        len <- min(len, length(lpmodel[[i]]))
      }
      name.list <- c(name.list, i)
    }
  }

  # Assign the component as the list of the bootstrap replications if they
  # are passed by the user. Otherwise, set the component as NULL.
  if (isTRUE(any.list)) {
    consol <- Map(list,
                  A.obs = lpmodel.extractlist(lpmodel$A.obs, len),
                  A.shp = lpmodel.extractlist(lpmodel$A.shp, len),
                  A.tgt = lpmodel.extractlist(lpmodel$A.tgt, len),
                  beta.obs = lpmodel.extractlist(lpmodel$beta.obs, len),
                  beta.shp = lpmodel.extractlist(lpmodel$beta.shp, len))
  } else {
    consol <- NULL
    len <- NULL
  }

  return(list(list = any.list,
              name = name.list,
              consol = consol,
              len = len))
}

#' Extracts the bootstrap replications of the \code{lpmodel} object
#'
#' @description This function extracts the bootstrap replications in the
#'   \code{lpmodel} object if the object is a \code{list}. Otherwise, this
#'   function returns \code{NULL}.
#'
#' @param obj A component in the \code{lpmodel} object.
#' @param len The length of the list.
#'
#' @return Returns one object.
#'   \item{result}{This is either a \code{list} or \code{NULL}.}
#'
#' @export
#'
lpmodel.extractlist <- function(obj, len) {
  if (inherits(obj, "list")) {
    result <- obj[-1]
  } else {
    result <- rep(list(NULL), len - 1)
  }

  return(result)
}

#' Combines deterministic components and one bootstrap estimate in
#' \code{lpmodel}
#'
#' @description This function is used in the bootstrap replications to combine
#'   the deterministic components in \code{lpmodel} to the stochastic component.
#'
#' @param lpm.de The deterministic components of the \code{lpmodel} object.
#' @param lpm.st An \code{lpmodel} object that only contains one bootstrap
#'   replication of the stochastic component(s). The deterministic component
#'   is set as \code{NULL}.
#'
#' @return Returns an \code{lpmodel} object that combines the deterministic
#'   and stochastic component.
#'   \item{lpm.de}{An updated \code{lpmodel} object.}
#'
#' @export
#'
lpmodel.update <- function(lpm.de, lpm.st) {
  # List of the names of the components in lpmodel
  lpmodel.namelist <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")

  # Combine the bootstrap estimate with the deterministic component
  for (i in seq_along(lpmodel.namelist)) {
    if (!is.null(lpm.st[[lpmodel.namelist[i]]])) {
      lpm.de[[lpmodel.namelist[i]]] <- lpm.st[[lpmodel.namelist[i]]]
    }
  }

  return(lpm.de)
}
