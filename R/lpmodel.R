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
