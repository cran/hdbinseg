% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbs.R
\name{sbs.alg}
\alias{sbs.alg}
\title{Sparsified Binary Segmentation}
\usage{
sbs.alg(
  x,
  cp.type = c(1, 2)[1],
  thr = NULL,
  trim = NULL,
  height = NULL,
  tau = NULL,
  temporal = TRUE,
  scales = NULL,
  diag = FALSE,
  B = 1000,
  q = 0.01,
  do.parallel = 4
)
}
\arguments{
\item{x}{input data matrix, with each row representing the component time series}

\item{cp.type}{\code{cp.type = 1} specifies change-points in the mean, \code{cp.type = 2} specifies change-points in the second-order structure}

\item{thr}{pre-defined threshold values; when \code{thr = NULL}, bootstrap procedure is employed for the threshold selection; when \code{thr != NULL} and \code{cp.type = 1}, \code{length(thr)} should match \code{nrow(x)}, if \code{cp.type = 2}, \code{length(thr)} should match \code{nrow(x)*(nrow(x)+1)/2*length(scales)}}

\item{trim}{length of the intervals trimmed off around the change-point candidates; \code{trim = NULL} activates the default choice (\code{trim = round(log(dim(x)[2]))})}

\item{height}{maximum height of the binary tree; \code{height = NULL} activates the default choice (\code{height = floor(log(dim(x)[2], 2)/2)})}

\item{tau}{a vector containing the scaling constant for each row of \code{x}; if \code{tau = NULL}, a data-driven choice is made which takes into account the presence of possibly multiple mean shifts and temporal dependence when \code{temporal = TRUE}}

\item{temporal}{used when \code{cp.type = 1}; if \code{temporal = FALSE}, rows of \code{x} are scaled by \link{mad} estimates, if \code{temporal = TRUE}, their long-run variance estimates are used}

\item{scales}{used when \code{cp.type = 2}; negative integers representing Haar wavelet scales to be used for computing \code{nrow(x)*(nrow(x)+1)/2} dimensional wavelet transformation of \code{x}; a small negative integer represents a fine scale}

\item{diag}{used when \code{cp.type = 2}; if \code{diag = TRUE}, only changes in the diagonal elements of the autocovariance matrices are searched for}

\item{B}{used when \code{is.null(thr)}; number of bootstrap samples for threshold selection}

\item{q}{used when \code{is.null(thr)}; quantile of bootstrap test statistics to be used for threshold selection}

\item{do.parallel}{used when \code{is.null(thr)}; number of copies of R running in parallel, if \code{do.parallel = 0}, \%do\% operator is used, see also \link{foreach}}
}
\value{
S3 \code{bin.tree} object, which contains the following fields:
   \item{tree}{a \link{list} object containing information about the nodes at which change-points are detected}
   \item{mat}{matrix concatenation of the nodes of \code{tree}}
   \item{ecp}{estimated change-points}
   \item{thr}{threshold used to construct the tree}
}
\description{
Perform the Sparsified Binary Segmentation algorithm detecting change-points in the mean or second-order structure of the data.
}
\examples{
x <- matrix(rnorm(20*300), nrow = 20)
sbs.alg(x, cp.type = 2, scales = -1, diag = TRUE, do.parallel = 0)$ecp
\donttest{
x <- matrix(rnorm(100*300), nrow = 100)
x[1:10, 151:300] <- x[1:10, 151:300]*sqrt(2)
sbs.alg(x, cp.type = 2, scales = -1, diag = TRUE, do.parallel = 0)$ecp
}
}
\references{
H. Cho and P. Fryzlewicz (2014) Multiple-change-point detection for high dimensional time series via sparsified binary segmentation. \emph{JRSSB}, vol. 77, pp. 475--507.
}
