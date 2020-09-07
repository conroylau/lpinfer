#' Create a mixed logit data generating process
#'
#' @description The model is:
#'  Y = 1(V1 + V2*W2 >= U), where
#'      W2 is an observable random variable
#'      V1 and V2 are unobservables with known support
#'      U is an unobservable with a standard logistic distribution
#' (V1, V2, U) are assumed independent of W.
#'
#' @param dimw Number of points in the support of W2.
#' @param dimv Number of points in the support of (V1, V2).
#' It should be the square of some integer, e.g. 4, 16, 100, 400, etc.
#' @param wbds The lower and upper bound of the support of W2.
#' @param vbds A list with two elements. The first element is a vector
#' indicating the lower and upper bound of the support of V1. The second is the
#' same for V2.
#'
#' @return dgp A list containing two matrices, vdist and wdist. Each matrix has
#' one row for each support point. The first element of the row is the
#' probability that this point is draw. The remaining elements are the actual
#' point of support.
#'
#' @export
mixedlogit_dgp <- function(dimw = 4, dimv = 16,
                           wbds = c(0, 3.0),
                           vbds = list(c(0, .5), c(-3, 0))) {

    dimvmarg <- sqrt(dimv)
    stopifnot(dimvmarg == round(dimvmarg)) # sqrt of dimv should be an integer

    v1supp <- seq(from = vbds[[1]][1], to = vbds[[1]][2], length.out = dimvmarg)
    v2supp <- seq(from = vbds[[2]][1], to = vbds[[2]][2], length.out = dimvmarg)
    vsupp <- expand.grid(v1 = v1supp, v2 = v2supp)
    vsupp <- as.matrix(vsupp)
    vdist <- cbind(pr = rep(1/nrow(vsupp), nrow(vsupp)), vsupp)

    w2supp <- seq(from = wbds[1], to = wbds[2], length.out = dimw)
    wsupp <- cbind(w1 = rep(1, length(w2supp)), w2 = w2supp)
    wdist <- cbind(pr = rep(1/nrow(wsupp), nrow(wsupp)), wsupp)

    return(list(vdist = vdist, wdist = wdist))
}

#' Generate the "Aobs" matrix for a mixedlogit dgp
#'
#' @param dgp Output from mixedlogit_dgp.
#'
#' @return Aobs a matrix with (i,j) element equal to the probability that Y = 1
#' if W is its ith support point and V is its jth support points.
#'
#' @export'
mixedlogit_Aobs <- function(dgp) {
    plogis(dgp$wdist[,-1] %*% t(dgp$vdist[,-1]))
}

#' Generate an elasticity distribution "Atgt" matrix for a mixedlogit dgp
#'
#' @description Atgt matrix corresponding to the probability that the
#' elasticity of choosing Y = 1 with respect to W2 is less than some point
#' eeval when W2 = w2eval, where both eeval and w2eval are parameters that can
#' be adjusted.
#'
#' @param dgp Output from mixedlogit_dgp.
#' @param w2eval The elasticity is calculated at w2eval.
#' @param eeval The probability is that the elasticity is less than eeval.
#'
#' @return A matrix (row vector) of dimension 1 x dimv.
#'
#' @export
mixedlogit_Atgt_dfelast <- function(dgp, w2eval = 1, eeval = -1) {
    vsupp <- dgp$vdist[,-1]
    pr <- plogis(vsupp %*% c(1, w2eval)) # choice probability by type
    dpr <- pr * (1 - pr) * vsupp[, 2] # deriv choice pr, e.g. pg. 58 Train 2009
    elasts <- dpr * w2eval / pr # elasticity
    return(as.numeric(elasts <= eeval))
}

#' Generate the "betaobs" function for the mixedlogit dgp
#'
#' @description This function takes a dataframe and the dgp object and returns
#' a vector P[Y = 1 | W2 = w2] for w2 in the support of W2.
#'
#' @param data Output from mixedlogit_draw.
#' @param dgp Output from mixedlogit_dgp.
#'
#' @return A vector of length dimw.
#'
#' @export
mixedlogit_betaobs <- function(data, dgp) {
    a <- aggregate(data$y, by = list(data$w), FUN = mean)
    print(a)
    bobs <- rep(0, nrow(dgp$wdist))
    bobs[a[,1]] <- a[,2]
    return(bobs)
}

#' Draw a sample of data from a mixedlogit dgp
#'
#' @param dgp Output from mixedlogit_dgp.
#' @param n Sample size.
#' @param seed The seed.
#'
#' @return A dataframe with two columns:
#'  \item{y}{Binary 0,1.}
#'  \item{w}{Integer from 1 up to dimw.}
#'
#' @export'
mixedlogit_draw <- function(dgp, n = 2000, seed = 1) {
    set.seed(seed)

    vint <- sample(1:nrow(dgp$vdist), n, replace = TRUE,
                   prob = dgp$vdist[,1])
    v <- dgp$vdist[vint, -1]

    wint <- sample(1:nrow(dgp$wdist), n, replace = TRUE,
                   prob = dgp$wdist[,1])

    idx <- rowSums(dgp$wdist[wint, -1] * v)
    u <- rlogis(n)
    y <- as.integer(u <= idx)

    data <- data.frame(y, w = wint)
}
