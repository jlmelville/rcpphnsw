## ensure module gets loaded
Rcpp::loadModule("HnswL2", TRUE)
Rcpp::loadModule("HnswCosine", TRUE)
Rcpp::loadModule("HnswIp", TRUE)

#' Find Nearest Neighbors and Distances
#'
#' A k-nearest neighbor algorithm using the HNSW library
#' (\url{https://github.com/nmslib/hnsw}).
#'
#' @section HNSW Parameters:
#'
#' Some details on the parameters used for index construction and search, based
#' on \url{https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md}:
#'
#' \itemize{
#'   \item \code{M} Controls the number of bi-directional links created for each
#'   element during index construction. Higher values lead to better results at
#'   the expense of memory consumption, which is around \code{M * 8-10} bytes
#'   per bytes per stored element. High intrinsic dimensionalities will require
#'   higher values of \code{M}. A range of \code{2 - 100} is typical, but
#'   \code{12 - 48} is ok for most use cases.
#'   \item \code{ef_construction} Size of the dynamic list used during
#'   construction. A larger value means a better quality index, but increases
#'   build time. Should be an integer value between 1 and the size of the
#'   dataset. A typical range is \code{100 - 2000}. Beyond a certain point,
#'   increasing \code{ef_construction} has no effect. A sufficient value of
#'   \code{ef_construction} can be determined by searching with \code{ef =
#'   ef_construction}, and ensuring that the recall is at least 0.9.
#'   \item \code{ef} Size of the dynamic list used during index search. Can
#'   differ from \code{ef_construction} and be any value between \code{k} (the
#'   number of neighbors sought) and the number of elements in the index being
#'   searched.
#' }
#'
#' @param X a numeric matrix of data to add. Each of the n rows is an item in
#'   the index.
#' @param k Number of neighbors to return.
#' @param distance Type of distance to calculate. One of:
#' \itemize{
#'   \item \code{"l2"} Squared L2, i.e. squared Euclidean.
#'   \item \code{"euclidean"} Euclidean.
#'   \item \code{"cosine"} Cosine.
#'   \item \code{"ip"} Inner product: 1 - sum(ai * bi), i.e. the cosine distance
#'   where the vectors are not normalized. This can lead to negative distances
#'   and other non-metric behavior.
#' }
#' @param include_self If \code{TRUE}, return the item itself as one of its
#'   \code{k}-neighbors.
#' @param M Controls the number of bi-directional links created for each element
#'   during index construction. Higher values lead to better results at the
#'   expense of memory consumption. Typical values are \code{2 - 100}, but
#'   for most datasets a range of \code{12 - 48} is suitable.
#' @param ef_construction Size of the dynamic list used during construction.
#'   A larger value means a better quality index, but increases build time.
#'   Should be an integer value between 1 and the size of the dataset.
#' @param ef Size of the dynamic list used during search. Higher values lead
#'   to improved recall at the expense of longer search time. Can take values
#'   between \code{k} and the size of the dataset and may be greater or smaller
#'   than \code{ef_construction}. Typical values are \code{100 - 2000}.
#' @param verbose If \code{TRUE}, log progress to the console.
#' @return a list containing:
#' \itemize{
#'   \item \code{idx} an n by k matrix containing the nearest neighbor indices.
#'   \item \code{dist} an n by k matrix containing the nearest neighbor
#'    distances.
#' }
#' @examples
#' get_knn(as.matrix(iris[, -5]), k = 10)
get_knn <- function(X, k = 10, distance = "euclidean", include_self = TRUE,
                    M = 16, ef_construction = 200, ef = ef_construction,
                    verbose = FALSE) {
  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  nr <- nrow(X)
  max_k <- ifelse(include_self, nr, nr - 1)
  if (k > max_k) {
    stop("k cannot be larger than ", max_k)
  }
  distance <- match.arg(distance, c("l2", "euclidean", "cosine", "ip"))

  ann <- hnsw_build(X = X, distance = distance, M = M, ef = ef_construction,
                    verbose = verbose)
  hnsw_search(X = X, ann = ann, k = k, include_self = include_self,
              ef = ef, verbose = verbose)
}

hnsw_build <- function(X, distance = "euclidean", M = 16, ef = 200,
                       verbose = FALSE) {
  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  distance <- match.arg(distance, c("l2", "euclidean", "cosine", "ip"))

  nr <- nrow(X)
  nc <- ncol(X)

  clazz <- switch(distance,
                  "l2" = RcppHNSW::HnswL2,
                  "euclidean" = RcppHNSW::HnswL2,
                  "cosine" = RcppHNSW::HnswCosine,
                  "ip" = RcppHNSW::HnswIp
                  )
  # Create the indexing object. You must say up front the number of items that
  # will be stored (nr).
  ann <- methods::new(clazz, nc, nr, M, ef)

  if (distance == "euclidean") {
    attr(ann, "distance") <- "euclidean"
  }

  tsmessage("Building HNSW index with metric '", distance, "'")
  progress <- Progress$new(max = nr, display = verbose)
  for (i in 1:nr) {
    # Items are added directly
    ann$addItem(X[i, ])
    progress$increment()
  }

  ann
}

hnsw_search <- function(X, ann, k, include_self = TRUE, ef = k,
                        verbose = FALSE) {
  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  nr <- nrow(X)

  max_k <- ifelse(include_self, nr, nr - 1)
  if (k > max_k) {
    stop("k cannot be larger than ", max_k)
  }

  if (!include_self) {
    k <- k + 1
  }

  idx <- matrix(nrow = nr, ncol = k)
  dist <- matrix(nrow = nr, ncol = k)

  ann$setEf(ef)
  tsmessage("Searching HNSW index")
  search_progress <- Progress$new(max = nr, display = verbose)
  for (i in 1:nr) {
    # Neighbors are queried by passing the vector back in
    # To get distances as well as indices, use include_distances = TRUE
    res <- ann$getNNsList(X[i, ], k, TRUE)
    idx[i, ] <- as.integer(res$item)
    dist[i, ] <- res$distance
    search_progress$increment()
  }

  if (!include_self) {
    idx <- idx[, -1, drop = FALSE]
    dist <- dist[, -1, drop = FALSE]
    k <- k - 1
  }

  if (!is.null(attr(ann, "distance")) &&
      attr(ann, "distance") == "euclidean") {
    dist <- sqrt(dist)
  }

  list(idx = idx, dist = dist)
}
