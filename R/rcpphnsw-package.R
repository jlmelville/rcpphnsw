#' Rcpp bindings for the hnswlib C++ library for approximate nearest neighbors.
#'
#' hnswlib is a library implementing the Hierarchical Navigable Small World
#' method for approximate nearest neighbor search.
#'
#' Details about hnswlib are available at the reference listed below.
#'
#' @docType package
#' @name RcppHnsw-package
#' @aliases HnswL2 Rcpp_HnswL2-class HnswCosine Rcpp_HnswCosine-class HnswIp Rcpp_HnswIp-class
#' @references
#' <https://github.com/nmslib/hnswlib>
#' @author James Melville for the R interface; Yury Malkov for hnswlib itself.
#'
#' Maintainer: James Melville <jlmelville@gmail.com>
#' @references
#' Malkov, Y. A., & Yashunin, D. A. (2016).
#' Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs.
#' *arXiv preprint* *arXiv:1603.09320*.
#' @useDynLib RcppHNSW, .registration = TRUE
#' @import Rcpp
#' @import methods
#' @export HnswL2
#' @export HnswCosine
#' @export HnswIp
NULL

## ensure module gets loaded
Rcpp::loadModule("HnswL2", TRUE)
Rcpp::loadModule("HnswCosine", TRUE)
Rcpp::loadModule("HnswIp", TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("RcppHNSW", libpath)
}
