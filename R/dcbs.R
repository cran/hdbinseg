#' Double CUSUM Binary Segmentation
#'
#' Perform the Double CUSUM Binary Segmentation algorithm detecting change-points in the mean or second-order structure of the data.
#' @param x input data matrix, with each row representing the component time series
#' @param cp.type \code{cp.type=1} specifies change-points in the mean, \code{cp.type=2} specifies change-points in the second-order structure
#' @param phi choice of parameter for weights in Double CUSUM statistic; 0 <= phi <= 1 or phi = -1 allowed with the latter leading to the DC statistic combining phi = 0 and phi = 1/2, see Section 4.1 of Cho (2016) for further details
#' @param thr pre-defined threshold values; when \code{thr = NULL}, bootstrap procedure is employed for the threshold selection; when \code{thr != NULL} and \code{cp.type = 1}, \code{length(thr)} should be one, if \code{cp.type = 2}, \code{length(thr)} should match \code{length(scales)}
#' @param trim length of the intervals trimmed off around the change-point candidates; \code{trim = NULL} activates the default choice (\code{trim = round(log(dim(x)[2]))})
#' @param height maximum height of the binary tree; \code{height = NULL} activates the default choice (\code{height = floor(log(dim(x)[2], 2)/2)})
#' @param temporal used when \code{cp.type = 1}; if \code{temporal = FALSE}, rows of \code{x} are scaled by \link{mad} estimates, if \code{temporal = TRUE}, their long-run variance estimates are used
#' @param scales used when \code{cp.type = 2}; negative integers representing Haar wavelet scales to be used for computing \code{nrow(x)*(nrow(x)+1)/2} dimensional wavelet transformation of \code{x}; a small negative integer represents a fine scale
#' @param diag used when \code{cp.type = 2}; if \code{diag = TRUE}, only changes in the diagonal elements of the autocovariance matrices are searched for
#' @param B used when \code{is.null(thr)}; number of bootstrap samples for threshold selection
#' @param q used when \code{is.null(thr)}; indicates the quantile of bootstrap test statistics to be used for threshold selection
#' @param do.parallel used when \code{is.null(thr)}; number of copies of R running in parallel, if \code{do.parallel = 0}, \%do\% operator is used, see also \link{foreach}
#' @return S3 \code{bin.tree} object, which contains the following fields:
#'    \item{tree}{a \link{list} object containing information about the nodes at which change-points are detected}
#'    \item{mat}{matrix concatenation of the nodes of \code{tree}}
#'    \item{ecp}{estimated change-points}
#'    \item{thr}{threshold used to construct the tree}
#' @references
#' H. Cho (2016) Change-point detection in panel data via double CUSUM statistic. \emph{Electronic Journal of Statistics}, vol. 10, pp. 2000--2038.
#' @examples
#' x <- matrix(rnorm(10*100), nrow=10)
#' dcbs.alg(x, cp.type=1, phi=.5, temporal=FALSE, do.parallel=0)$ecp
#' \donttest{
#' x <- matrix(rnorm(100*300), nrow=100)
#' x[1:10, 151:300] <- x[1:10, 151:300] + 1
#' dcbs.alg(x, cp.type=1, phi=-1, temporal=FALSE, do.parallel=0)$ecp
#' }
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats ar rgeom mad acf cor
#' @useDynLib hdbinseg, .registration = TRUE
#' @export
dcbs.alg <- function(x, cp.type=c(1, 2)[1], phi=.5, thr=NULL, trim=NULL, height=NULL,
                     temporal=TRUE, scales=NULL, diag=FALSE,
                     B=1000, q=.01, do.parallel=4){

  n <- dim(x)[1]; T <- dim(x)[2]
  stopifnot(n >= 1 && T >=1)

  stopifnot((phi >= 0 && phi <= 1) || phi==-1)
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
  } else stopifnot(length(thr)==1)

  if(is.null(trim)) trim <- round(log(T)^1.5)
  if(is.null(height)) height <- floor(log(T, 2)/2)

  if(cp.type==1){
    tree <- list(matrix(0, 5, 0))
    mat <- matrix(NA, ncol=0, nrow=6)
    ccp <- clean.cp(x, type='dcbs', phi=phi, trim=trim, height=height)
    cx <- ccp$x
    if(!temporal) tau <- apply(cx, 1, mad) else tau <- apply(cx, 1, mad)/apply(cx, 1, function(z){(1-acf(z, lag.max=1, type='correlation', plot=FALSE)$acf[2,,1])})
    input <- x/tau

    br <- dcbs.branch(input, phi, 1, dim(input)[2], trim)
    if(is.null(thr)){
      ecp <- sort(ccp$mat[4, ccp$mat[1, ] <= 2])
      int <- cbind(c(1, ecp+1), c(ecp, T))[which.max(diff(c(0, ecp, T))), ]
      thr <- dcbs.thr(cx, interval = int, phi = phi, cp.type=1, temporal=temporal, B=B, q=q, do.parallel=do.parallel)
    }

    if(br[5] > thr){
      br[1] <- 1
      tree[[1]] <- matrix(c(tree[[1]], matrix(br, 5, 1)), 5, 1)
      mat <- cbind(mat, c(1, br))

      j <- 1
      while(length(tree)==j & j < height){
        npc <- dim(tree[[j]])[2]
        ncc <- 0; i <- 1
        while(i <= npc){
          if(tree[[j]][3, i]-tree[[j]][2, i]+1 > 2*trim+1){
            s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
            br <- dcbs.branch(input, phi, s, e, trim)
            if(br[5] > thr){
              br[1] <- 2*tree[[j]][1, i]-1
              if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
              ncc <- ncc+1
              tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
              tree[[j+1]][, ncc] <- br
              mat <- cbind(mat, c(j+1, br))
            }
          }
          if(tree[[j]][4, i]-tree[[j]][3, i] > 2*trim+1){
            s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
            br <- dcbs.branch(input, phi, s, e, trim)
            if(br[5] > thr){
              br[1] <- 2*tree[[j]][1, i]
              if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
              ncc <- ncc+1
              tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
              tree[[j+1]][, ncc] <- br
              mat <- cbind(mat, c(j+1, br))
            }
          }
          i <- i+1
        }
        j <- j+1
      }
    }
  }

  if(cp.type==2){
    sgn <- sign(cor(t(x)))
    input <- gen.input(x, scales, TRUE, diag, sgn)
    tree <- list(matrix(0, 5, 0))
    mat <- matrix(NA, ncol=0, nrow=6)

    br <- dcbs.branch(input, phi, 1, dim(input)[2], trim)
    if(is.null(thr)){
      mt <- dcbs.make.tree(input, phi, trim, 2)
      ecp <- sort(mt$mat[4, ])
      int <- cbind(c(1, ecp+1), c(ecp, T))[which.max(diff(c(0, ecp, T))), ]
      thr <- dcbs.thr(x, interval = int, phi = phi, cp.type=2, scales=scales, diag=diag, sgn=sgn, B=B, q=q, do.parallel=do.parallel)
    }

    if(br[5] > thr){
      br[1] <- 1
      tree[[1]] <- matrix(c(tree[[1]], matrix(br, 5, 1)), 5, 1)
      mat <- cbind(mat, c(1, br))

      j <- 1
      while(length(tree)==j & j < height){
        npc <- dim(tree[[j]])[2]
        ncc <- 0; i <- 1
        while(i <= npc){
          if(tree[[j]][3, i]-tree[[j]][2, i]+1 > 2*trim+1){
            s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
            br <- dcbs.branch(input, phi, s, e, trim)
            if(br[5] > thr){
              br[1] <- 2*tree[[j]][1, i]-1
              if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
              ncc <- ncc+1
              tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
              tree[[j+1]][, ncc] <- br
              mat <- cbind(mat, c(j+1, br))
            }
          }
          if(tree[[j]][4, i]-tree[[j]][3, i] > 2*trim+1){
            s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
            br <- dcbs.branch(input, phi, s, e, trim)
            if(br[5] > thr){
              br[1] <- 2*tree[[j]][1, i]
              if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
              ncc <- ncc+1
              tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
              tree[[j+1]][, ncc] <- br
              mat <- cbind(mat, c(j+1, br))
            }
          }
          i <- i+1
        }
        j <- j+1
      }
    }
  }

  rownames(mat) <- c('level_index', 'child_index', 's', 'b', 'e', 'test_stat')

  return(structure(list(tree=tree, mat=mat, ecp=sort(mat[4, ], decreasing=FALSE), thr=thr), class='bin.tree'))

}

#' Growing a branch for DCBS algorithm
#'
#' Grow a branch of the binary tree for the Double CUSUM Binary Segmentation without thresholding
#' @param input input data matrix, with each row representing the component time series or their transformation
#' @param s,e start and end of the interval of consideration at a given iteration
#' @param phi,trim see \code{\link{dcbs.alg}}
#' @return a vector containing the information about the branch, including the location of the estimated change-point and the corresponding test statistic
#' @import Rcpp
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
dcbs.branch <- function(input, phi, s, e, trim){

  N <- dim(input)[1]
  if(phi >= 0){
    stat <- func_dc(input[, s:e], phi)$res
  } else {
    stat <- apply(func_dc(input[, s:e], 0)$mat*log(N) + func_dc(input[, s:e], 0.5)$mat, 2, max)
  }
  test.stat <- max(stat[-c((1:trim), (e-s-trim+1):(e-s))])
  b <- s+min(which(stat==test.stat))-1

  c(NA, s, b, e, test.stat)

}

#' Growing a binary tree for DCBS algorithm
#'
#' Grow a binary tree of a given height via Double CUSUM Binary Segmentation without thresholding
#' @param input input data matrix, with each row representing the component time series or their transformation
#' @param phi, trim, height see \code{\link{dcbs.alg}}
#' @return S3 \code{bin.tree} object, which contains the following fields:
#'    \item{tree}{a \link{list} object containing information about the nodes at which change-points are detected}
#'    \item{mat}{matrix concatenation of the nodes of \code{tree}}
#'    \item{ecp}{estimated change-points}
#'    \item{thr}{threshold used to construct the tree}
#' @import Rcpp
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
dcbs.make.tree <- function(input, phi=.5, trim=NULL, height=NULL){

  len <- dim(input)[2]
  if(is.null(height)) height <- floor(log(len, 2)/2)
  if(is.null(trim)) trim <- round(log(len)^1.6)

  tree <- list(matrix(0, 5, 1))
  mat <- matrix(NA, ncol=0, nrow=6)

  br <- dcbs.branch(input, phi, 1, len, trim)
  br[1] <- 1
  tree[[1]][, 1] <- br
  mat <- cbind(mat, c(1, br))

  j <- 1
  while(length(tree)==j & j < height){
    npc <- dim(tree[[j]])[2]
    ncc <- 0; i <- 1
    while(i <= npc){
      if(tree[[j]][3, i]-tree[[j]][2, i]+1 > 2*trim+1){
        s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
        br <- dcbs.branch(input, phi, s, e, trim)
        br[1] <- 2*tree[[j]][1, i]-1
        if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
        ncc <- ncc+1
        tree[[j+1]] <- matrix(c(tree[[j+1]], br), 5, ncc)
        mat <- cbind(mat, c(j+1, br))
      }
      if(tree[[j]][4, i]-tree[[j]][3, i] > 2*trim+1){
        s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
        br <- dcbs.branch(input, phi, s, e, trim)
        br[1] <- 2*tree[[j]][1, i]-1
        if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 1)))
        ncc <- ncc+1
        tree[[j+1]] <- matrix(c(tree[[j+1]], br), 5, ncc)
        mat <- cbind(mat, c(j+1, br))
      }
      i <- i+1
    }
    j <- j+1
  }
  rownames(mat) <- c('level_index', 'child_index', 's', 'b', 'e', 'test_stat')

  return(structure(list(tree=tree, mat=mat, ecp=sort(mat[4, ], decreasing=FALSE), thr=0), class='bin.tree'))

}

#' Bootstrapping for threshold selection in DCBS algorithm
#'
#' Generate thresholds for DCBS algorithm via bootstrapping
#' @param z input data matrix, with each row representing the component time series
#' @param interval a vector of two containing the start and the end points of the interval from which the bootstrap test statistics are to be calculated
#' @param phi,cp.type,temporal,scales,diag,B,q,do.parallel see \code{\link{dcbs.alg}}
#' @param do.clean.cp if \code{do.clean.cp = TRUE} pre-change-point cleaning is performed
#' @param sgn if \code{diag = FALSE}, wavelet transformations of the cross-covariances are computed with the matching signs
#' @return a numeric value for the threshold
#' @import Rcpp foreach doParallel parallel iterators
#' @importFrom stats ar rgeom mad cor
#' @useDynLib hdbinseg, .registration = TRUE
#' @export
dcbs.thr <- function(z, interval = c(1, dim(z)[2]), phi=.5, cp.type=1, do.clean.cp=FALSE, temporal=TRUE, scales=NULL, diag=FALSE, sgn=NULL, B=1000, q=.01, do.parallel=4){

  if(do.parallel > 0){
    cl <- parallel::makeCluster(do.parallel); doParallel::registerDoParallel(cl)
  }
  `%mydo%` <- ifelse(do.parallel > 0, `%dopar%`, `%do%`)

  n <- dim(z)[1]; len <- dim(z)[2]
  s <- interval[1]; e <- interval[2]
  if(e - s + 1 < round(log(len)^1.5)) stop("Error: the interval too short!")
  boot.method <- 'fm'
  if(cp.type==1) mby <- round(log(n))
  if(cp.type==2){
    if(is.null(sgn)) sgn <- sign(cor(t(z[, s:e])))
    mby <- round(log(length(scales)*n*(n+1)/2))
  }
  tby <- 2

  if(cp.type==1 && do.clean.cp) z <- clean.cp(z, type='dcbs', phi=phi, trim=round(log(len)^1.5), height=round(log(len, 2)/2))$x

  if(boot.method=='fm'){
    gfm <- get.factor.model(z[, s:e])
    r <- gfm$r.hat
    if(r > 0){
      lam <- gfm$eigvec[, 1:r, drop=FALSE]*sqrt(n)
      f <- t(gfm$eigvec[, 1:r, drop=FALSE])%*%gfm$norm.x/sqrt(n)
      idio <- gfm$norm.x - lam%*%f

      burnin <- 100
      o <- floor(min(sqrt(e-s+1), 10*log(e-s+1, 10)))
      ar.coef <- matrix(0, nrow=r, ncol=o)
      residuals <- matrix(0, nrow=r, ncol=e-s+1-o)
      for(i in 1:r){
        ar.fit <- stats::ar(f[i, ], order.max=o, method='yw')
        if(length(ar.fit$ar) > 0) ar.coef[i, 1:length(ar.fit$ar)] <- ar.fit$ar
        residuals[i, ] <- ar.fit$resid[-(1:o)]
      }
      if(sum(apply(abs(ar.coef), 2, max) > 0) > 0){
        ar.coef <- ar.coef[, 1:max(which(apply(abs(ar.coef), 2, max) > 0)), drop=FALSE]; o <- ncol(ar.coef)
      } else o <- 0

      p <- min(.5, 1/mean(apply(idio, 1, function(w){g <- get.gg(w); ((g[2]/g[1])^2)^(1/3)*len^(1/5)})))
    } else boot.method <- 'sb'
  }
  if(boot.method=='sb') p <- min(.5, 1/mean(apply(z[, s:e], 1, function(w){g <- get.gg(w); ((g[2]/g[1])^2)^(1/3)*len^(1/5)})))

  ns <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo','hdbinseg')) %mydo% {
    if(boot.method=='fm'){
      bf <- residuals[, sample(dim(residuals)[2], len+o+burnin, replace=TRUE), drop=FALSE]
      if(o > 0) for(t in (o+1):(len+o+burnin)) bf[, t] <- bf[, t] + apply(ar.coef*bf[, t-(1:o)], 1, sum)
      bf <- bf[, -(1:(o+burnin))]
      bc <- lam%*%bf

      ind <- c()
      while(length(ind) < len){
        L <- stats::rgeom(1, p); I <- sample(e-s+1, 1)
        ind <- c(ind, rep(1:(e-s+1), 1+ceiling(L/(e-s+1)))[I:(I+L-1)])
      }
      ind <- ind[1:len]
      bi <- idio[, ind]

      bx <- bc + bi
      bx <- bx*gfm$sd
    }
    if(boot.method=='sb'){
      ind <- c()
      while(length(ind) < len){
        L <- stats::rgeom(1, p); I <- sample(e-s+1, 1)
        ind <- c(ind, rep(1:(e-s+1), 1+ceiling(L/(e-s+1)))[I:(I+L-1)])
      }
      ind <- ind[1:len]
      bx <- z[, ind]
    }
    if(cp.type==1){
      if(!temporal) tau <- apply(bx, 1, mad) else tau <- apply(bx, 1, mad)/apply(bx, 1, function(z){(1-acf(z, lag.max=1, type='correlation', plot=FALSE)$acf[2,,1])})
      input <- bx/tau
    }
    if(cp.type==2) input <- gen.input(bx, scales, TRUE, diag, sgn)

    if(phi >= 0){
      stat <- max(func_dc_by(input, phi, mby, tby)$res)
    } else {
      stat <- max(func_dc_by(input, 0, round(log(dim(input)[1])), 2)$mat*log(dim(input)[1]) + func_dc_by(input, 0.5, round(log(dim(input)[1])), 2)$mat)
    }
    stat
  }
  if(do.parallel > 0) parallel::stopCluster(cl)

  quantile(ns, 1-q)
}
