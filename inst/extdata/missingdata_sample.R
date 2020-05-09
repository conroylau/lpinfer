#' Simulated data
#'
#' @format A data frame with 1,000 rows and 2 columns.
#' \describe{
#'   \item{Y}{multivariate discrete outcome variable taht takes value of 0 to 1
#'   with step size 0.1}
#'   \item{D}{binary treatment where \eqn{D_i=1} means \eqn{Y_i} is observed}
#' }
#'
#' @source Simulated
#'
sampledata_generation <- function(){
  wd <- getwd()
  sampledata = read.csv(paste0(wd, "sampledata.csv"))
  colnames(sampledata) <- c("D", "Y")
  library("devtools")
  use_data(sampledata, overwrite = TRUE)
  sampledata = sampledata
  return(sampledata)
}

wd <- getwd()
sampledata = read.csv(paste0(wd, "sampledata.csv"))
colnames(sampledata) <- c("D", "Y")
library("devtools")
use_data(sampledata, overwrite = TRUE)
