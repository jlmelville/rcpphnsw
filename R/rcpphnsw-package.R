#' Approximate nearest neighbor search with hnswlib
#'
#' RcppHNSW provides an R interface to hierarchical navigable small-world
#' graphs. Use [hnsw_knn()] for neighbors within one dataset, or [hnsw_build()]
#' and [hnsw_search()] to query new data against an existing index.
#'
#' @section Module classes:
#'
#' `HnswEuclidean`, `HnswL2`, `HnswCosine`, and `HnswIp` provide direct access
#' to insertion, search, deletion, resizing, and saved indexes. They use
#' Euclidean distance, squared Euclidean distance, one minus cosine similarity,
#' and one minus inner product, respectively. An index returned by [hnsw_build()]
#' is one of these objects, so you can use its methods directly too. See the
#' [Module guide](https://jlmelville.github.io/rcpphnsw/articles/module-api-index-lifecycle.html)
#' for constructors and methods.
#'
#' Labels start at one and follow insertion order. Deleted items are excluded
#' from search and retrieval but still count towards `size()` and capacity.
#' If insertion fails after modifying the index, or a native resize fails,
#' discard the index and rebuild or reload it.
#'
#' @section Numeric limits:
#'
#' Coordinates are stored as single-precision floats and must be finite and
#' representable in that format. Cosine vectors are normalized to unit length
#' and must have a nonzero norm after conversion.
#'
#' To prevent distance overflow, each vector's sum of absolute coordinates is
#' limited to approximately `4.6e18`, after conversion and cosine normalization.
#' Loaded checkpoints must meet the same limit, including deleted items.
#'
#' @section Saved indexes:
#'
#' `ann$save()` writes hnswlib's raw checkpoint format, whose compatibility
#' depends on the hnswlib version and platform. Load with the original
#' dimension and distance class. You can also load an `HnswL2` checkpoint into
#' `HnswEuclidean`, or vice versa.
#'
#' Loading restores items, deletion state, and capacity, but resets search
#' `ef` to 10. Set `ef` again with `setEf()` or [hnsw_search()]. An optional
#' positive load capacity overrides the saved value when it is at least the
#' stored item count; smaller values retain the saved capacity.
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
