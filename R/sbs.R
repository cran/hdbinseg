#' Sparsified Binary Segmentation
#'
#' Perform the Sparsified Binary Segmentation algorithm detecting change-points in the mean or second-order structure of the data.
#' @param x input data matrix, with each row representing the component time series
#' @param cp.type \code{cp.type=1} specifies change-points in the mean, \code{cp.type=2} specifies change-points in the second-order structure
#' @param thr pre-defined threshold values; when \code{thr = NULL}, bootstrap procedure is employed for the threshold selection; when \code{thr != NULL} and \code{cp.type = 1}, \code{length(thr)} should match \code{nrow(x)}, if \code{cp.type = 2}, \code{length(thr)} should match \code{nrow(x)*(nrow(x)+1)/2*length(scales)}
#' @param trim length of the intervals trimmed off around the change-point candidates; \code{trim = NULL} activates the default choice (\code{trim = round(log(dim(x)[2]))})
#' @param height maximum height of the binary tree; \code{height = NULL} activates the default choice (\code{height = floor(log(dim(x)[2], 2)/2)})
#' @param temporal used when \code{cp.type = 1}; if \code{temporal = FALSE}, rows of \code{x} are scaled by \link{mad} estimates, if \code{temporal = TRUE}, their long-run variance estimates are used
#' @param scales used when \code{cp.type = 2}; negative integers representing Haar wavelet scales to be used for computing \code{nrow(x)*(nrow(x)+1)/2} dimensional wavelet transformation of \code{x}; a small negative integer represents a fine scale
#' @param diag used when \code{cp.type = 2}; if \code{diag = TRUE}, only changes in the diagonal elements of the autocovariance matrices are searched for
#' @param B used when \code{is.null(thr)}; number of bootstrap samples for threshold selection
#' @param q used when \code{is.null(thr)}; quantile of bootstrap test statistics to be used for threshold selection
#' @param do.parallel used when \code{is.null(thr)}; number of copies of R running in parallel, if \code{do.parallel = 0}, \%do\% operator is used, see also \link{foreach}
#' @return S3 \code{bin.tree} object, which contains the following fields:
#'    \item{tree}{a \link{list} object containing information about the nodes at which change-points are detected}
#'    \item{mat}{matrix concatenation of the nodes of \code{tree}}
#'    \item{ecp}{estimated change-points}
#'    \item{thr}{threshold used to construct the tree}
#' @references
#' H. Cho and P. Fryzlewicz (2014) Multiple-change-point detection for high dimensional time series via sparsified binary segmentation. \emph{JRSSB}, vol. 77, pp. 475--507.
#' @examples
#' x <- matrix(rnorm(20*300), nrow=20)
#' sbs.alg(x, cp.type=2, scales=-1, diag=TRUE, do.parallel=0)$ecp
#' \donttest{
#' x <- matrix(rnorm(100*300), nrow=100)
#' x[1:10, 151:300] <- x[1:10, 151:300]*sqrt(2)
#' sbs.alg(x, cp.type=2, scales=-1, diag=TRUE, do.parallel=0)$ecp
#' }
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats ar rgeom mad quantile acf cor
#' @useDynLib hdbinseg, .registration = TRUE
#' @export
sbs.alg <- function(x, cp.type=c(1, 2)[1], thr=NULL, trim=NULL, height=NULL,
                    temporal=TRUE, scales=NULL, diag=FALSE,
                    B=1000, q=.01, do.parallel=4){

  n <- dim(x)[1]; T <- dim(x)[2]
  stopifnot(n >= 1 && T >=1)

  stopifnot(cp.type==1 || cp.type==2)
  stopifnot(is.null(trim) || (trim > 0 && 2*trim+1 < T))
  stopifnot(is.null(height) || 2^height < T)
  if(cp.type==2){
    if(!is.null(scales)){
      stopifnot(sum(scales<0)==length(scales))
      scales <- sort(scales, decreasing=TRUE)
    } else scales <- -1
  }

  if(is.null(thr)){
    stopifnot(B>0)
    stopifnot(q >= 0 && q <= 1)
  } else{
    if(cp.type==1) stopifnot(length(thr)==n)
    if(cp.type==2) stopifnot(length(thr)==(length(scales)*n*(n+1)/2*(!diag) + length(scales)*n*diag))
  }

  if(is.null(trim)) trim <- round(log(T)^1.5)
  if(is.null(height)) height <- floor(log(T, 2)/2)

  if(cp.type==1){
    ccp <- clean.cp(x, type='sbs', trim=trim, height=height)
    cx <- ccp$x
    if(!temporal) tau <- apply(cx, 1, mad) else tau <- apply(cx, 1, mad)/apply(cx, 1, function(z){1-acf(z, lag.max=1, type='correlation', plot=FALSE)$acf[2,,1]})
    if(is.null(thr)) thr <- sbs.thr(cx, interval = c(1, T), cp.type=1, B=B, q=q/n, do.parallel=do.parallel)
    mt <- sbs.make.tree(x, tau, thr, trim, height)
  }

  if(cp.type==2){
    sgn <- sign(cor(t(x)))
    input <- gen.input(x, scales, TRUE, diag, sgn)
    if(is.null(thr)) thr <- sbs.thr(x, interval = c(1, T), cp.type=2, scales=scales, diag=diag, sgn=sgn, B=B, q=q, do.parallel=do.parallel)
    mt <- sbs.make.tree(input, rep(1, dim(input)[1]), thr, trim, height)
  }

  return(mt)

}

#' Growing a binary tree for SBS algorithm
#'
#' Grow a binary tree via Sparsified Binary Segmentation
#' @param input input data matrix, with each row representing the component time series or their transformation
#' @param tau scaling terms the rows of \code{input}
#' @param thr,trim,height see \code{\link{sbs.alg}}
#' @return S3 \code{bin.tree} object, which contains the following fields:
#'    \item{tree}{a \link{list} object containing information about the nodes at which change-points are detected}
#'    \item{mat}{matrix concatenation of the nodes of \code{tree}}
#'    \item{ecp}{estimated change-points}
#'    \item{thr}{threshold used to construct the tree}
#' @import Rcpp
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
sbs.make.tree <- function(input, tau=rep(1, nrow(input)), thr, trim, height){

  N <- dim(input)[1]; len <- dim(input)[2]
  if(is.null(height)) height <- floor(log(len, 2)/2)
  if(is.null(trim)) trim <- round(log(len)^1.5)

  tree <- list(matrix(0, 5, 1))
  mat <- matrix(NA, ncol=0, nrow=6)

  acs <- func_cusum(input)$acs
  stat <- apply(acs*(acs>thr)/tau, 2, sum)
  stat[c((1:trim), (len-trim):len)] <- 0
  sb <- search.b(stat, trim)

  if(sb$FLAG){
    tree[[1]][1, 1] <- 1
    tree[[1]][2, 1] <- 1
    tree[[1]][3, 1] <- sb$b
    tree[[1]][4, 1] <- len
    tree[[1]][5, 1] <- sb$test.stat
    mat <- cbind(mat, c(1, tree[[1]][, ]))

    j <- 1
    while(length(tree)==j & j < height){
      npc <- dim(tree[[j]])[2]
      ncc <- 0; i <- 1
      while(i <= npc){
        if(tree[[j]][3, i]-tree[[j]][2, i]+1 > 2*trim+1){
          s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
          acs <- func_cusum(input[, s:e])$acs
          stat <- apply(acs*(acs>thr)/tau, 2, sum)
          stat[c((1:trim), (e-s-trim+1):(e-s))] <- 0
          sb <- search.b(stat, trim)

          if(sb$FLAG){
            if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]-1
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- s+sb$b-1
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- sb$test.stat
            mat <- cbind(mat, c(j+1, tree[[j+1]][, ncc]))
          }
        }
        if(tree[[j]][4, i]-tree[[j]][3, i] > 2*trim+1){
          s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
          acs <- func_cusum(input[, s:e])$acs
          stat <- apply(acs*(acs>thr)/tau, 2, sum)
          stat[c((1:trim), (e-s-trim+1):(e-s))] <- 0
          sb <- search.b(stat, trim)

          if(sb$FLAG){
            if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- s+sb$b-1
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- sb$test.stat
            mat <- cbind(mat, c(j+1, tree[[j+1]][, ncc]))
          }
        }
        i <- i+1
      }
      j <- j+1
    }
  }
  rownames(mat) <- c('level_index', 'child_index', 's', 'b', 'e', 'test_stat')

  return(structure(list(tree=tree, mat=mat, ecp=sort(mat[4, ], decreasing=FALSE), thr=thr), class='bin.tree'))

}

#' Searching for a change-point in a branch
#'
#' Searches for a change-point in a branch of a binary tree grown by the Sparsified Binary Segmentation
#' @param stat aggregated CUSUM statistics
#' @param trim see \code{\link{sbs.alg}}
#' @return a list containing
#'    \item{b}{a possible location of the change-point}
#'    \item{test.stat}{test statistic}
#'    \item{FLAG}{should the branch be grown?}
#' @import Rcpp
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
search.b <- function(stat, trim){
  FLAG <- FALSE; b <- test.stat <- NA
  while(!FLAG && sum(stat) > 0){
    test.stat <- max(stat)
    b <- min(which(stat==test.stat))
    int <- max(b-trim+1, 1):min(b+trim, length(stat))
    if(sum(stat[int]==0) > 0){
      stat[int] <- 0
    } else FLAG <- TRUE
  }
  return(list(b=b, test.stat=test.stat, FLAG=FLAG))
}

#' Bootstrapping for threshold selection in SBS algorithm
#'
#' Generate thresholds for SBS algorithm via bootstrapping
#' @param z input data matrix, with each row representing the component time series
#' @param interval a vector of two containing the start and the end points of the interval from which the bootstrap test statistics are to be calculated
#' @param cp.type,scales,diag,B,q,do.parallel see \code{\link{sbs.alg}}
#' @param do.clean.cp if \code{do.clean.cp = TRUE} pre-change-point cleaning is performed
#' @param sgn if \code{diag = FALSE}, wavelet transformations of the cross-covariances are computed with the matching signs
#' @return
#' if \code{cp.type = 1}, a vector of length \code{nrow(z)}, each containing the threshold applied to the CUSUM statistics from the corresponding coordinate of \code{z}
#' if \code{cp.type = 2}, a vector of length \code{length(scales)*nrow(z)} (when \code{diag = TRUE}) or \code{length(scales)*nrow(z)*(nrow(z)+1)/2} (when \code{diag = FALSE}), each containing the threshold applied to the CUSUM statistics of the corresponding wavelet transformation of \code{z}
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats cor median ar rgeom mad
#' @useDynLib hdbinseg, .registration = TRUE
#' @export
sbs.thr <- function(z, interval = c(1, dim(z)[2]), cp.type=1, do.clean.cp=TRUE, scales=NULL, diag=FALSE, sgn=NULL, B=1000, q=.01, do.parallel=4){

  if(do.parallel > 0){
    cl <- parallel::makeCluster(do.parallel); doParallel::registerDoParallel(cl)
  }
  `%mydo%` <- ifelse(do.parallel > 0, `%dopar%`, `%do%`)

  n <- dim(z)[1]; len <- dim(z)[2]
  s <- interval[1]; e <- interval[2]
  if(e - s + 1 < ceiling(log(len)^1.5)) stop("Error: the interval too short!")
  if(cp.type==1) d <- n
  if(cp.type==2){
    if(is.null(sgn)) sgn <- sign(cor(t(z[, s:e])))
    d <- n*diag+n*(n+1)/2*(!diag)
  }

  if(cp.type==1 && do.clean.cp) z <- clean.cp(z, type='sbs', trim=round(log(len)^1.5), height=round(log(len, 2)/2))$x

  burnin <- 100
  p <- floor(min(sqrt(e-s+1), 10*log(e-s+1, 10)))
  ar.coef <- matrix(0, nrow=d, ncol=p)
  residuals <- matrix(0, nrow=d, ncol=e-s+1-p)
  k <- 1
  for(i in 1:n){
    ar.fit <- stats::ar(z[i, s:e], order.max=p, method='yw')
    if(length(ar.fit$ar) > 0) ar.coef[k, 1:length(ar.fit$ar)] <- ar.fit$ar
    residuals[k, ] <- ar.fit$resid[-(1:p)]
    k <- k+1
    if(cp.type==2 && i < n && !diag){
      for(j in (i+1):n){
        ar.fit <- stats::ar(z[i, s:e]-sgn[i, j]*z[j, s:e], order.max=p, method='yw')
        if(length(ar.fit$ar) > 0) ar.coef[k, 1:length(ar.fit$ar)] <- ar.fit$ar
        residuals[k, ] <- ar.fit$resid[-(1:p)]
        k <- k+1
      }
    }
  }
  if(sum(apply(abs(ar.coef), 2, max) > 0) > 0){
    ar.coef <- ar.coef[, 1:max(which(apply(abs(ar.coef), 2, max) > 0)), drop=FALSE]; p <- ncol(ar.coef)
  } else p <- 0

  ns <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo','hdbinseg')) %mydo% {
    bx <- residuals[, sample(dim(residuals)[2], len+p+burnin, replace=TRUE)]
    if(p > 0) bx <- func_mvt_ar(ar.coef, bx)
    bx <- bx[, -(1:(p+burnin))]

    if(cp.type==1){
      tmp <- apply(func_cusum(bx)$acs, 1, max)
    }
    if(cp.type==2){
      input <- gen.input(bx, scales, TRUE, TRUE, sgn)
      tmp <- apply(func_cusum(input)$acs, 1, max)
    }
    tmp
  }
  if(do.parallel > 0)  parallel::stopCluster(cl)

  apply(ns, 1, function(z){quantile(z, 1-q)})
}
