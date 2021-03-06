#' Find Nearest Neighbors and Distances
#'
#' A k-nearest neighbor algorithm using the hnswlib library
#' (\url{https://github.com/nmslib/hnswlib}).
#'
#' @section Hnswlib Parameters:
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
#' @param X a numeric matrix of data to search Each of the n rows is an item in
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
#' @param M Controls the number of bi-directional links created for each element
#'   during index construction. Higher values lead to better results at the
#'   expense of memory consumption. Typical values are \code{2 - 100}, but
#'   for most datasets a range of \code{12 - 48} is suitable. Can't be smaller
#'   than 2.
#' @param ef_construction Size of the dynamic list used during construction.
#'   A larger value means a better quality index, but increases build time.
#'   Should be an integer value between 1 and the size of the dataset.
#' @param ef Size of the dynamic list used during search. Higher values lead
#'   to improved recall at the expense of longer search time. Can take values
#'   between \code{k} and the size of the dataset and may be greater or smaller
#'   than \code{ef_construction}. Typical values are \code{100 - 2000}.
#' @param verbose If \code{TRUE}, log messages to the console.
#' @param progress If \code{"bar"} (the default), also log a progress bar when
#'   \code{verbose = TRUE}. There is a small but noticeable overhead (a few
#'   percent of run time) to tracking progress. Set \code{progress = NULL} to
#'   turn this off. Has no effect if \code{verbose = FALSE}.
#' @param n_threads Maximum number of threads to use. The exact number is
#'   determined by \code{grain_size}.
#' @param grain_size Minimum amount of work to do (rows in \code{X} to add or
#'   search for) per thread. If the number of rows in \code{X} isn't sufficient,
#'   then fewer than \code{n_threads} will be used. This is useful in cases
#'   where the overhead of context switching with too many threads outweighs
#'   the gains due to parallelism.
#' @return a list containing:
#' \itemize{
#'   \item \code{idx} an n by k matrix containing the nearest neighbor indices.
#'   \item \code{dist} an n by k matrix containing the nearest neighbor
#'    distances.
#' }
#' Every item in the dataset is considered to be a neighbor of itself, so the
#' first neighbor of item \code{i} should always be \code{i} itself. If that
#' isn't the case, then any of \code{M}, \code{ef_construction} and \code{ef}
#' may need increasing.
#' @examples
#' iris_nn_data <- hnsw_knn(as.matrix(iris[, -5]), k = 10)
#' @references
#' Malkov, Y. A., & Yashunin, D. A. (2016).
#' Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs.
#' \emph{arXiv preprint} \emph{arXiv:1603.09320}.
hnsw_knn <- function(X, k = 10, distance = "euclidean",
                     M = 16, ef_construction = 200, ef = 10,
                     verbose = FALSE, progress = "bar", n_threads = 0,
                     grain_size = 1) {
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && n_threads >= 0)
  stopifnot(is.numeric(grain_size) && length(grain_size) == 1 && grain_size >= 0)

  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  if (M < 2) {
    stop("M cannot be < 2")
  }
  ef_construction <- max(ef_construction, k)

  nr <- nrow(X)
  max_k <- nr
  if (k > max_k) {
    stop("k cannot be larger than ", max_k)
  }
  distance <- match.arg(distance, c("l2", "euclidean", "cosine", "ip"))

  ann <- hnsw_build(
    X = X, distance = distance, M = M, ef = ef_construction,
    verbose = verbose, progress = progress, n_threads = n_threads,
    grain_size = grain_size
  )
  hnsw_search(
    X = X, ann = ann, k = k, ef = ef, verbose = verbose,
    progress = progress, n_threads = n_threads, grain_size = grain_size
  )
}

#' Build an hnswlib nearest neighbor index
#'
#' @param X a numeric matrix of data to add. Each of the n rows is an item in
#'   the index.
#' @param distance Type of distance to calculate. One of:
#' \itemize{
#'   \item \code{"l2"} Squared L2, i.e. squared Euclidean.
#'   \item \code{"euclidean"} Euclidean.
#'   \item \code{"cosine"} Cosine.
#'   \item \code{"ip"} Inner product: 1 - sum(ai * bi), i.e. the cosine distance
#'   where the vectors are not normalized. This can lead to negative distances
#'   and other non-metric behavior.
#' }
#' @param M Controls the number of bi-directional links created for each element
#'   during index construction. Higher values lead to better results at the
#'   expense of memory consumption. Typical values are \code{2 - 100}, but
#'   for most datasets a range of \code{12 - 48} is suitable. Can't be smaller
#'   than 2.
#' @param ef Size of the dynamic list used during construction.
#'   A larger value means a better quality index, but increases build time.
#'   Should be an integer value between 1 and the size of the dataset.
#' @param verbose If \code{TRUE}, log messages to the console.
#' @param progress If \code{"bar"} (the default), also log a progress bar when
#'   \code{verbose = TRUE}. There is a small but noticeable overhead (a few
#'   percent of run time) to tracking progress. Set \code{progress = NULL} to
#'   turn this off. Has no effect if \code{verbose = FALSE}.
#' @param n_threads Maximum number of threads to use. The exact number is
#'   determined by \code{grain_size}.
#' @param grain_size Minimum amount of work to do (rows in \code{X} to add) per
#'   thread. If the number of rows in \code{X} isn't sufficient, then fewer than
#'   \code{n_threads} will be used. This is useful in cases where the overhead
#'   of context switching with too many threads outweighs the gains due to
#'   parallelism.
#' @return an instance of a \code{HnswL2}, \code{HnswCosine} or \code{HnswIp}
#'   class.
#' @examples
#' irism <- as.matrix(iris[, -5])
#' ann <- hnsw_build(irism)
#' iris_nn <- hnsw_search(irism, ann, k = 5)
hnsw_build <- function(X, distance = "euclidean", M = 16, ef = 200,
                       verbose = FALSE, progress = "bar", n_threads = 0,
                       grain_size = 1) {
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && n_threads >= 0)
  stopifnot(is.numeric(grain_size) && length(grain_size) == 1 && grain_size >= 0)

  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  if (M < 2) {
    stop("M cannot be < 2")
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

  tsmessage(
    "Building HNSW index with metric '", distance, "'",
    " ef = ", formatC(ef), " M = ", formatC(M), " using ", n_threads, " threads"
  )
  ann$setNumThreads(n_threads)
  ann$setGrainSize(grain_size)

  nstars <- 50
  if (verbose && nr > nstars && !is.null(progress) && progress == "bar") {
    progress_for(
      nr, nstars,
      function(chunk_start, chunk_end) {
        ann$addItems(X[chunk_start:chunk_end, , drop = FALSE])
      }
    )
  }
  else {
    ann$addItems(X)
  }

  tsmessage("Finished building index")
  ann
}

#' Search an hnswlib nearest neighbor index
#'
#' @param X A numeric matrix of data to search for neighbors.
#' @param ann an instance of a \code{HnswL2}, \code{HnswCosine} or \code{HnswIp}
#'   class.
#' @param k Number of neighbors to return. This can't be larger than the number
#'   of items that were added to the index \code{ann}. To check the size of the
#'   index, call \code{ann$size()}.
#' @param ef Size of the dynamic list used during search. Higher values lead
#'   to improved recall at the expense of longer search time. Can take values
#'   between \code{k} and the size of the dataset. Typical values are
#'   \code{100 - 2000}.
#' @param verbose If \code{TRUE}, log messages to the console.
#' @param progress If \code{"bar"} (the default), also log a progress bar when
#'   \code{verbose = TRUE}. There is a small but noticeable overhead (a few
#'   percent of run time) to tracking progress. Set \code{progress = NULL} to
#'   turn this off. Has no effect if \code{verbose = FALSE}.
#' @param n_threads Maximum number of threads to use. The exact number is
#'   determined by \code{grain_size}.
#' @param grain_size Minimum amount of work to do (rows in \code{X} to search)
#'   per thread. If the number of rows in \code{X} isn't sufficient, then fewer
#'   than \code{n_threads} will be used. This is useful in cases where the
#'   overhead of context switching with too many threads outweighs the gains due
#'   to parallelism.
#' @return a list containing:
#' \itemize{
#'   \item \code{idx} an n by k matrix containing the nearest neighbor indices.
#'   \item \code{dist} an n by k matrix containing the nearest neighbor
#'    distances.
#' }
#'
#' Every item in the dataset is considered to be a neighbor of itself, so the
#' first neighbor of item \code{i} should always be \code{i} itself. If that
#' isn't the case, then any of \code{M}, \code{ef_construction} and \code{ef}
#' may need increasing.
#'
#' @examples
#' irism <- as.matrix(iris[, -5])
#' ann <- hnsw_build(irism)
#' iris_nn <- hnsw_search(irism, ann, k = 5)
hnsw_search <- function(X, ann, k, ef = 10, verbose = FALSE, progress = "bar",
                        n_threads = 0, grain_size = 1) {
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && n_threads >= 0)
  stopifnot(is.numeric(grain_size) && length(grain_size) == 1 && grain_size >= 0)

  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  nr <- nrow(X)

  ef <- max(ef, k)

  idx <- matrix(nrow = nr, ncol = k)
  dist <- matrix(nrow = nr, ncol = k)

  ann$setEf(ef)
  ann$setNumThreads(n_threads)
  ann$setGrainSize(grain_size)
  tsmessage(
    "Searching HNSW index with ef = ", formatC(ef), " and ", n_threads,
    " threads"
  )

  nstars <- 50
  if (verbose && nr > nstars && !is.null(progress) && progress == "bar") {
    progress_for(
      nr, nstars,
      function(chunk_start, chunk_end) {
        res <- ann$getAllNNsList(X[chunk_start:chunk_end, , drop = FALSE], k, TRUE)
        idx[chunk_start:chunk_end, ] <<- as.integer(res$item)
        dist[chunk_start:chunk_end, ] <<- res$distance
      }
    )
  }
  else {
    res <- ann$getAllNNsList(X, k, TRUE)
    idx <- res$item
    dist <- res$distance
  }

  if (!is.null(attr(ann, "distance")) &&
    attr(ann, "distance") == "euclidean") {
    dist <- sqrt(dist)
  }

  tsmessage("Finished searching")
  list(idx = idx, dist = dist)
}
