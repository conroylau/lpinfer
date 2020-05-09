#' Evaluates an object inside \code{lpmodel}
#' 
#' @description This function returns the matrix or vector depending on the 
#'   class of the variable in the \code{lpmodel} object. In the design of
#'   the \code{lpinfer} module, objects in \code{lpmodel} can have three
#'   different types of classes:
#'   \itemize{
#'      \item{\code{function} --- If the object is a function, then the
#'         object that is evaluated by the function with the data will be
#'         returned.}
#'      \item{\code{list} --- If the object is a list, then the \eqn{i}th
#'         element will be returned.}
#'      \item{\code{matrix} or \code{numeric} --- If the object is a matrix
#'         or in a vector, then it will be directly returned.}
#'   }
#' 
#' @param data Data frame
#' @param obj Object in lpmodel
#' @param i Index 
#' 
#' @return Returns an object at iteration \code{i}.
#' @export
#' 
lpmodel.eval <- function(data, obj, i){
  if (class(obj) == "function"){
    obj.eval <- obj(data)
  } else if (class(obj) == "list"){
    obj.eval <- obj[[i]]
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
#'      \item{\code{list} --- If the object is a list, then the \eqn{i}th
#'         element will be returned.}
#'   }
#' 
#' @param data Data frame
#' @param obj Object in lpmodel
#' @param i Index 
#' 
#' @return Returns the point estimate and the asymptotic variane of the 
#'    \code{beta.obs} object.
#'    \item{beta.obs}{Point estimate of \eqn{\widehat{\bm{\beta}}_{\rm obs}}}
#'    \item{beta.obs}{Asymptotic variance of 
#'       \eqn{\widehat{\bm{\beta}}_{\rm obs}}}
#'    
#' @export
#' 
lpmodel.beta.eval <- function(data, obj, i){
  if (class(obj) == "function"){
    beta.return <- obj(data)
    if (class(beta.return) == "list"){
      if (length(beta.return) == 2){
        if (is.null(nrow(beta.return[[1]]))){
          beta.obs.hat <- beta.return[[1]]
          omega.hat <- beta.return[[2]]
        } else if (nrow(beta.return[[1]] == 1)){
          beta.obs.hat <- beta.return[[1]]
          omega.hat <- beta.return[[2]]      
        } else {
          beta.obs.hat <- beta.return[[2]]
          omega.hat <- beta.return[[1]] 
        }
      } else if (length(beta.return) == 1){
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
  } else if (class(obj) == "list"){
    if (class(obj[[i]]) == "list"){
      if (length(obj[[i]]) != 2){
        stop(paste0("When the first object of the list 'beta.obs' is a ",
                    "list, it needs to have two objects"))
      } else if (length(obj[[i]]) == 1){
        beta.obs.hat <- obj[[i]][[1]]
        omega.hat <- NULL
      } else {
        if (is.null(nrow(obj[[i]][[1]]))){
          beta.obs.hat <- obj[[i]][[1]]
          omega.hat <- obj[[i]][[2]]
        } else if (nrow(obj[[i]][[1]]) == 1 | ncol(obj[[i]][[1]]) == 1){
          beta.obs.hat <- obj[[i]][[1]]
          omega.hat <- obj[[i]][[2]]      
        } else {
          beta.obs.hat <- obj[[i]][[2]]
          omega.hat <- obj[[i]][[1]] 
        }
      }
    } else {
      beta.obs.hat <- obj[[i]]
      omega.hat <- NULL
    }
  } else if (class(obj) == "numeric"){
    beta.obs.hat <- obj
    omega.hat <- NULL
  }
  
  return(list(beta.obs = beta.obs.hat,
              omega = omega.hat))
}

#' Print the \code{lpmodel} object
#' 
#' @description This function prints objects that are contained in the 
#'    list of \code{lpmodel}.
#'    
#' @param x The \code{lpmodel} object.
#' 
#' @return Print the summary of the objects in \code{lpmodel}.
#' 
#' @export
#' 
print.lpmodel <- function(x, data = NULL, ...){
  # List of variables 
  lpmodel.string <- c("A.obs", "A.shp", "A.tgt", "beta.obs", "beta.shp")
  lpmodel.ind <- NULL
  for (i in 1:length(lpmodel.string)){
    if (!is.null(x[[lpmodel.string[i]]])){
      lpmodel.ind <- c(lpmodel.ind, i)
    }
  }
  if (length(lpmodel.ind) == 0){
    cat("'lpmodel' object does not contain the required objects.")
  } else {
    cat("Object     Class \tDimension \tLength \n")
    for (i in 1:length(lpmodel.string)){
      if (i %in% lpmodel.ind){
        obj <- x[[lpmodel.string[i]]]
        # Check class of object 
        class.tmp <- class(obj)
        # Check length of object 
        if (class.tmp == "list"){
          class.tmp <- "list  "
          length.tmp <- length(obj)
          dimension.str <- dim(as.matrix(obj[[1]]))
          dimension.tmp <- paste0(dimension.str[1], "x", dimension.str[2])
        } else if (class.tmp == "function"){
          # If data is not passed, print "N/A" for dimensions. Otherwise,
          # compute the output for the output object
          if (is.null(data)){
            length.tmp <- "N/A"
            dimension.tmp <- "N/A"           
          } else {
            tmp.obj <- obj(data)
            if (class(tmp.obj) == "list"){
              length.tmp <- length(tmp.obj)
            } else {
              length.tmp <- 1
            }
            dimension.str <- dim(as.matrix(tmp.obj))
            dimension.tmp <- paste0(dimension.str[1], "x", dimension.str[2])
          }
        } else if (class.tmp %in% c("data.frame", "matrix", "numeric")){
          length.tmp <- 1
          dimension.str <- dim(as.matrix(obj))
          dimension.tmp <- paste0(dimension.str[1], "x", dimension.str[2])
        } else {
          length.tmp <- length(obj)
          dimension.tmp <- "  "
        }
        # Check 
        if (i <= 3){
          cat(sprintf("%s      %s\t%s\t\t%s\n", 
                      lpmodel.string[i], class.tmp, dimension.tmp, length.tmp))
        } else {
          cat(sprintf("%s   %s\t%s\t\t%s\n", 
                      lpmodel.string[i], class.tmp, dimension.tmp, length.tmp)) 
        }
      } else {
        if (i <= 3){
          cat(sprintf("%s      -empty-\t-empty-\t\t-empty-\n", 
                      lpmodel.string[i]))
        } else {
          cat(sprintf("%s   -empty-\t-empty-\t\t-empty-\n", 
                      lpmodel.string[i])) 
        }
      }
    }
  }
}

#' Define a \code{lpmodel} object
#' 
#' @description This function defines the objects required in the 
#'    \code{lpinfer} module in the \code{lpmodel} class.
#'    
#' @param A.obs A matrix, list or function
#' @param A.shp A matrix, list or function
#' @param A.tgt A matrix, list or function
#' @param beta.obs A vector, list or function
#' @param beta.shp A vector, list or function
#' 
#' @return Returns a list of \code{lpmodel} objects in the \code{lpmodel} 
#'    class.
#'    
#' @export
#' 
lpmodel <- function(A.obs = NULL, A.shp = NULL, A.tgt = NULL, beta.obs = NULL, 
                    beta.shp = NULL){
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

#' Summary of results from \code{lpmodel}
#' 
#' @description This function uses the summary method on the return list of the
#'    function \code{lpmodel}. This is a wrapper of the \code{print} command.
#'    
#' @param x The \code{lpmodel} object.
#' @param ... Additional arguments.
#' 
#' @return Print the summary of the basic set of results from \code{lpmodel}.
#' 
#' @export
#' 
summary.lpmodel <- function(x, ...){
  print(x)
}
