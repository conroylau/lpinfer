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
  } else if (class(obj) == "list"){
    if (is.null(nrow(beta.return[[i]][[1]]))){
      beta.obs.hat <- beta.return[[i]][[1]]
      omega.hat <- beta.return[[i]][[2]]
    } else if (nrow(beta.return[[i]][[1]] == 1)){
      beta.obs.hat <- beta.return[[i]][[1]]
      omega.hat <- beta.return[[i]][[2]]      
    } else {
      beta.obs.hat <- beta.return[[i]][[2]]
      omega.hat <- beta.return[[i]][[1]] 
    }
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
            length.tmp <- length(tmp.obj)
            dimension.str <- dim(as.matrix(tmp.obj[[1]]))
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
# attr(lpmodel, "class") <- "lpmodel"



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
