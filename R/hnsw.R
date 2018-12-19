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
#' @param n_threads Number of threads to use. Default is half that recommended
#'   by RcppParallel.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads}.
#' @param verbose If \code{TRUE}, log progress to the console.
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
hnsw_knn <- function(X, k = 10, distance = "euclidean",
                    M = 16, ef_construction = 200, ef = ef_construction,
                    n_threads = 0,
                    grain_size = 1,
                    verbose = FALSE) {
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

  ann <- hnsw_build(X = X, distance = distance, M = M, ef = ef_construction,
                    verbose = verbose)
  hnsw_search(X = X, ann = ann, k = k, ef = ef,
              n_threads = n_threads, grain_size = grain_size,
              verbose = verbose)
}

#' Build a nearest neighor index
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
#' @param verbose If \code{TRUE}, log progress to the console.
#' @return an instance of a \code{HnswL2}, \code{HnswCosine} or \code{HnswIp}
#'   class.
#' @examples
#' irism <- as.matrix(iris[, -5])
#' ann <- hnsw_build(irism)
#' iris_nn <- hnsw_search(irism, ann, k = 5)
hnsw_build <- function(X, distance = "euclidean", M = 16, ef = 200,
                       verbose = FALSE) {
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

  tsmessage("Building HNSW index with metric '", distance, "'",
            " ef = ", formatC(ef), " M = ", formatC(M))
  progress <- Progress$new(max = nr, display = verbose)
  for (i in 1:nr) {
    # Items are added directly
    ann$addItem(X[i, ])
    progress$increment()
  }

  ann
}

#' Search an HNSW nearest neighbor index
#'
#' @param X A numeric matrix of data to search for neighbors.
#' @param ann an instance of a \code{HnswL2}, \code{HnswCosine} or \code{HnswIp}
#'   class.
#' @param k Number of neighbors to return.
#' @param ef Size of the dynamic list used during search. Higher values lead
#'   to improved recall at the expense of longer search time. Can take values
#'   between \code{k} and the size of the dataset. Typical values are
#'   \code{100 - 2000}.
#' @param n_threads Number of threads to use. Default is half that recommended
#'   by RcppParallel.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads}.
#' @param verbose If \code{TRUE}, log progress to the console.
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
hnsw_search <- function(X, ann, k, ef = k,
                        n_threads = 0,
                        grain_size = 1,
                        verbose = FALSE) {
  if (n_threads > 0 ) {
    res <- hnsw_search_parallel(X = X, ann = ann, k = k, ef = ef,
                         n_threads = n_threads, grain_size = grain_size,
                         verbose = verbose)
  }
  else {
    res <- hnsw_search_serial(X = X, ann = ann, k = k, ef = ef,
                              verbose = verbose)
  }
  if (!is.null(attr(ann, "distance")) &&
      attr(ann, "distance") == "euclidean") {
    res$dist <- sqrt(res$dist)
  }
  res
}

hnsw_search_serial <- function(X, ann, k, ef = k, verbose = FALSE) {
  if (!is.matrix(X)) {
    stop("X must be matrix")
  }
  nr <- nrow(X)

  max_k <- nr
  if (k > max_k) {
    stop("k cannot be larger than ", max_k)
  }

  ef <- max(ef, k)

  idx <- matrix(nrow = nr, ncol = k)
  dist <- matrix(nrow = nr, ncol = k)

  ann$setEf(ef)
  tsmessage("Searching HNSW index with ef = ", formatC(ef))
  search_progress <- Progress$new(max = nr, display = verbose)
  for (i in 1:nr) {
    # Neighbors are queried by passing the vector back in
    # To get distances as well as indices, use include_distances = TRUE
    res <- ann$getNNsList(X[i, ], k, TRUE)
    idx[i, ] <- as.integer(res$item)
    dist[i, ] <- res$distance
    search_progress$increment()
  }

  list(idx = idx, dist = dist)
}

hnsw_search_parallel <- function(X, ann, k, ef = k,
                                 n_threads =
                                   max(1, RcppParallel::defaultNumThreads() / 2),
                                 grain_size = 1,
                                 verbose = FALSE) {
  index_file <- tempfile()
  ann$save(index_file)

  tsmessage("Searching HNSW index with ef = ", formatC(ef),
            " using ", n_threads, " thread", ifelse(n_threads != 1, "s", ""))
  RcppParallel::setThreadOptions(numThreads = n_threads)

  ann_class <- class(ann)
  search_nn_func <- switch(ann_class,
                           Rcpp_HnswCosine = hnsw_cosine_nns,
                           Rcpp_HnswIp = hnsw_ip_nns,
                           Rcpp_HnswL2 = hnsw_l2_nns,
                           stop("BUG: unknown Annoy class '", ann_class, "'")
  )

  res <- search_nn_func(index_name = index_file,
                        mat = X, search_k = k,
                        ef = ef,
                        grain_size = grain_size,
                        verbose = verbose)
  unlink(index_file)
  res
}
