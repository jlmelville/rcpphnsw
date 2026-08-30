#' Rcpp bindings for the hnswlib C++ library for approximate nearest neighbors.
#'
#' hnswlib is a library implementing the Hierarchical Navigable Small World
#' method for approximate nearest neighbor search.
#'
#' Details about hnswlib are available at the reference listed below.
#'
#' @section Data and index contract:
#'
#' RcppHNSW converts coordinates to single-precision floating-point values for
#' storage and computation. Coordinates must be finite and representable in
#' that format. Cosine items and queries must have a positive finite norm after
#' conversion.
#'
#' Module indexes assign one-based labels in insertion order. `size()` reports
#' the total number of items added, including items marked deleted. Deleted
#' items remain allocated, are excluded from search, and cannot be returned by
#' `getItems()`. Consequently, `k` must be positive and no larger than the
#' active count: total items minus deleted items. If an exception escapes after
#' insertion has begun, the index fails closed and must be discarded and
#' rebuilt or reloaded.
#'
#' A thread setting of zero or one uses serial execution; zero grain size is an
#' alias for one. Parallel construction may be nondeterministic even with a
#' fixed hnswlib seed. R's `set.seed()` does not control hnswlib.
#'
#' @section Raw index checkpoints:
#'
#' The Module `save()` method and filename constructors expose hnswlib's raw
#' index format. These files are operational checkpoints, not a stable
#' cross-version or cross-platform archive. Load with the exact original
#' dimension and normally the same distance class. L2 and Euclidean classes
#' may reinterpret each other's same-width checkpoints; no corresponding
#' semantic guarantee is made for cosine or inner-product reinterpretation.
#'
#' Raw checkpoints assume the contiguous labels created by RcppHNSW. Loading
#' restores deletion state but resets search `ef` to hnswlib's default of 10;
#' call `setEf()` or `hnsw_search()` to select another value. An optional load
#' capacity may enlarge the index, but cannot be smaller than the stored item
#' count.
#'
#' @docType package
#' @name RcppHnsw-package
#' @aliases HnswL2 Rcpp_HnswL2-class HnswCosine Rcpp_HnswCosine-class HnswIp
#' @aliases Rcpp_HnswIp-class HnswEuclidean Rcpp_HnswEuclidean-class
#' @aliases RcppHNSW-package
#' @references
#' <https://github.com/nmslib/hnswlib>
#' @author James Melville for the R interface; Yury Malkov for hnswlib itself.
#'
#' Maintainer: James Melville <jlmelville@gmail.com>
#' @references
#' Malkov, Y. A., & Yashunin, D. A. (2020). Efficient and robust approximate
#' nearest neighbor search using Hierarchical Navigable Small World graphs.
#' *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 42(4),
#' 824-836. \doi{10.1109/TPAMI.2018.2889473}.
"_PACKAGE"

## ensure module gets loaded
Rcpp::loadModule("HnswL2", TRUE)
Rcpp::loadModule("HnswCosine", TRUE)
Rcpp::loadModule("HnswIp", TRUE)
Rcpp::loadModule("HnswEuclidean", TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("RcppHNSW", libpath)
}
