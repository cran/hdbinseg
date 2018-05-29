#' Generating Haar wavelet transformation of time series
#' 
#' Generate Haar wavelet transformation of time series containing information about possible change-points in the second-order structure of time series
#' @param x input data matrix, with each row representing the component time series
#' @param scales negative integers for wavelet scales, with a small negative integer representing a fine scale
#' @param sq if \code{sq = TRUE}, squared root of wavelet periodograms are used for change-point analysis
#' @param diag if \code{diag = TRUE}, only changes in the diagonal elements of the autocovariance matrices are searched for
#' @param sgn if \code{diag = FALSE}, wavelet transformations of the cross-covariances are computed with the matching signs
#' @return matrix of the square root of scaled Haar wavelet periodograms of \code{x}
#' @import Rcpp
#' @importFrom stats cor
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
gen.input <- function(x, scales, sq, diag, sgn=NULL){
  
  if(is.null(sgn)) sgn <- sign(cor(t(x)))
  input <- NULL
  for(sc in scales){
    cc <- func_coef(x, sc)
    if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
    input <- rbind(input, t(func_input(cc, sgn, sq, diag)))
  }
  
  input
  
}

#' Removing change-points in the mean
#' 
#' (Over-)estimate and remove change-points in the mean for scale estimation and bootstrapping
#' @param z input data matrix, with each row representing the component time series
#' @param type if \code{type = 'dcbs'}, a binary tree of given \code{height} is grown using DCBS algorithm without thresholding, if \code{type = 'sbs'}, the binary trees is grown using the SBS algorithm with thresholds chosen small
#' @param phi,trim,height see \code{\link{dcbs.alg}}
#' @return a list containing
#' \item{x}{\code{z} with potential change-points in the mean removed}
#' \item{z}{\code{mat} object of an S3 \code{bin.tree} object}
#' @import Rcpp
#' @importFrom stats median
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
clean.cp <- function(z, type=c('dcbs', 'sbs'), phi=.5, trim=NULL, height=NULL){
  
  len <- dim(z)[2]; n <- dim(z)[1]
  
  if(is.null(trim)) trim <- round(log(len))
  if(is.null(height)) height <- floor(log(len, 2)/2)
  
  if(type=='dcbs') mat <- dcbs.make.tree(z, phi, trim, height)$mat
  if(type=='sbs') mat <- sbs.make.tree(z, rep(1, n), apply(func_cusum(z)$acs, 1, median), trim, 2)$mat
  
  brks <- sort(unique(c(0, mat[4, ] , len)), decreasing=FALSE)
  for(l in 1:(length(brks)-1)){
    int <- (brks[l]+1):brks[l+1]		
    z[, int] <- z[, int]-apply(z[, int, drop=FALSE], 1, mean)
  }
  
  return(list(x=z, mat=mat))
}

#' Flat-top kernel
#' 
#' @param h bandwidth
#' @param len desired length of the vector containing the flat-top kernel
#' @return a vector containing half of the flat-top kernel
#' @keywords internal
tri.kern <- function(h, len=0){
  if(len==0) len <- h+1
  filter <- rep(0, len)
  i <- 1
  while (i <= len) {
    u <- (i - 1)/h
    if (u < 1/2) 
      filter[i] <- 1
    if (u >= 1/2 & u <= 1) 
      filter[i] <- 2 * (1 - u)
    if (u > 1) 
      break
    i <- i + 1
  }
  filter
}

#' Data-driven selection of the average block size
#' 
#' Computes quantities required for data-driven selection of the average block size for stationary bootstrap
#' @param z univariate time series
#' @param M bandwidth, by default, \code{M = NULL} to be automatically selected as in Politis and White (2004)
#' @param C \code{C = 2} is the default value chosen for automatic bandwidth selection in Politis and White (2004)
#' @param max.K \code{max.K = 5} is the default value chosen for automatic bandwidth selection in Politis and White (2004)
#' @return Estimates for the two quantities required for average block size selection.
#' @references 
#' D. N. Politis and H. White (2004) Automatic block-length selection for the dependent bootstrap. \emph{Econometric Reviews}, vol. 23, pp. 53--70.
#' @importFrom stats acf
#' @keywords internal
get.gg <- function(z, M=NULL, C=2, max.K=5){
  len <- length(z)
  max.K <- max(max.K, sqrt(log(len)))
  acv <- stats::acf(z, type="covariance", lag.max=len-1, plot=FALSE)$acf[,,1]
  if(is.null(M)){
    l <- 1; ind <- 0
    while(l < sqrt(len)){
      if(abs(acv[l+1])/acv[1] < C*sqrt(log(len)/len)){
        ind <- ind+1
      } else{
        if(ind>0) ind <- 0
      }
      if(ind==max.K) break
      l <- l+1
    }
    lam <- max(1/2, l-max.K); M <- 2*lam
  }
  k <- tri.kern(M)
  c(acv[1]+2*sum(k[-1]*acv[2:(M+1)]), 2*sum(k[-1]*(1:M)*acv[2:(M+1)]))
}

#' Factor model estimation via Principal Component Analysis
#' 
#' Estimates the components of the factor structure for an input time series, such as loadings and factors, as well as estimating the number of factors.
#' @param z input data matrix, with each row representing the component time series
#' @param r.max maximum allowed number of factors
#' @param ic estimator for the factor number; \code{if(ic=='ah')} eigenvalue ratio estimator of Ahn and Horenstein (2013) is used, \code{if(ic=='bn')} information criterion estimator of Bai and Ng (2002) is used
#' @param ic.op type of the estimator specified by \code{ic}
#' @param normalisation \code{if(normalisation==TRUE)} rows of \code{z} are normalised prior to factor analysis
#' @return a list containing
#'   \item{eigvec}{eigenvectors of \code{z}}
#'   \item{eigval}{eigenvalues of \code{z}}
#'   \item{norm.x}{row-wise normalised \code{z} \code{if(normalisation==TRUE)}}
#'   \item{r.hat}{estimated number of factors}
#'   \item{r.max}{maximum number of factors used}
#'   \item{ic.eval}{vector containing information criterion evaluated at \code{r = 0, 1, ..., r.max}}
#'   \item{mean}{row-wise means of \code{z}}
#'   \item{sd}{row-wise standard deviations of \code{z}}
#'   \item{ic}{input \code{ic}}
#'   \item{ic.op}{input \code{ic.op}}
#' @references
#' S. C. Ahn and A. R. Horenstein (2013) Eigenvalue ratio test for the number of factors. \emph{Econometrica}, vol. 81, pp. 1203--1227.
#' J. Bai and S. Ng (2002) Determining the number of factors in approximate factor models. \emph{Econometrica}, vol. 70, pp. 191--221.
#' @useDynLib hdbinseg, .registration = TRUE
#' @keywords internal
#' @export
get.factor.model <- function(z, r.max=NULL, ic=c('ah', 'bn')[1], ic.op=2, normalisation=TRUE){
  
  T <- ncol(z); n <- nrow(z)
  cnt <- min(n, T)
  if(is.null(r.max)) r.max <- round(sqrt(cnt))
  
  if(normalisation){
    mean.z <- apply(z, 1, mean); sd.z <- apply(z, 1, mad)
    sz <- z - matrix(rep(mean.z, T), byrow=FALSE, nrow=n)
    sz <- sz/sd.z
  } else{
    mean.z <- rep(0, n); sd.z <- rep(1, n); sz <- z
  }
  zz <- sz%*%t(sz)/T
  eig <- eigen(zz, symmetric=TRUE)
  eigvec <- eig$vectors[, eig$values > 0, drop=FALSE]
  eigval <- eig$values[eig$values > 0]
  
  ic.eval <- rep(0, 1+r.max)
  
  if(ic=='bn'){
    ic.eval[1] <- log(sum(eigval))
    pen <- (ic.op==1)*(n+T)/(n*T)*log(n*T/(n+T)) + (ic.op==2)*(n+T)/(n*T)*log(cnt) + (ic.op==3)*log(cnt)/cnt
    ic.eval[1+1:r.max] <- log(sum(eigval)-cumsum(eigval[1:r.max]))+(1:r.max)*pen
    r.hat <- min(which(ic.eval==min(ic.eval)))-1
  } else if(ic=='ah'){
    eigval0 <- min(eigval[1]+1, sum(eigval)/log(cnt))
    if(ic.op==1){ 
      ic.eval <- c(eigval0, eigval[1:r.max])/eigval[1:(r.max+1)]
    } 
    if(ic.op==2){
      v <- c(sum(eigval), sum(eigval) - cumsum(eigval[1:(r.max+1)]))
      mu <- c(eigval0, eigval)[1:(r.max+2)]/v
      ic.eval <- log(1+mu[1:(r.max+1)])/log(1+mu[-1])
    }
    r.hat <- min(which(ic.eval==max(ic.eval)))-1
  }
  
  return(list(eigvec=eigvec, eigval=eigval, norm.x=sz, r.hat=r.hat, r.max=r.max, ic.eval=ic.eval, 
              mean=mean.z, sd=sd.z, ic=ic, ic.op=ic.op))
}
