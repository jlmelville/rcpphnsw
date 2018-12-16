#' Rcpp bindings for the HNSW C++ library for approximate nearest neighbors.
#'
#' HNSW is a library implementing the Hierarchical Navigable Small World method
#' for approximate nearest neighbor search.
#'
#' Details about HNSW are available at the reference listed below.
#'
#' @docType package
#' @name RcppHnsw-package
#' @aliases HnswL2 Rcpp_HnswL2-class HnswCosine Rcpp_HnswCosine-class HnswIp Rcpp_HnswIp-class
#' @references
#' \url{https://github.com/nmslib/hnsw}
#' @author James Melville for the R inteface; Yury Malkov for HNSW itself.
#'
#' Maintainer: James Melville <jlmelville@gmail.com>
NULL

## ensure module gets loaded
Rcpp::loadModule("HnswL2", TRUE)
Rcpp::loadModule("HnswCosine", TRUE)
Rcpp::loadModule("HnswIp", TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("RcppHNSW", libpath)
}
