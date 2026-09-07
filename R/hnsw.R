#' Find approximate nearest neighbors
#'
#' Build an HNSW index and find `k` neighbors for every item in `X`.
#'
#' @param X Numeric matrix with one item per row, or per column when
#'   `byrow = FALSE`.
#' @param k Number of neighbors per item. A positive whole number no larger
#'   than the number of items in `X`.
#' @param distance Distance to use:
#'   * `"euclidean"`: Euclidean distance.
#'   * `"l2"`: squared Euclidean distance.
#'   * `"cosine"`: one minus cosine similarity.
#'   * `"ip"`: one minus inner product, `1 - sum(a * b)`; can be negative.
#' @param M Number of graph links per item during construction. Larger values
#'   improve connectivity and use more memory. Must be between 2 and 10000.
#' @param ef_construction Candidate-list size during construction.
#'   Larger values improve index quality and increase build time. Must be a
#'   positive whole number; raised to at least `k`.
#' @param ef Candidate-list size during search. Increase it to
#'   improve recall at the cost of search time. Must be a positive whole number;
#'   the effective value is at least `k` and is not capped at the index size.
#' @param verbose If `TRUE`, log messages to the console.
#' @param progress Unused; retained for compatibility.
#' @param n_threads Maximum number of threads for batch insertion or search.
#'   Zero (the default) and one run serially.
#' @param grain_size Minimum number of items per thread. Larger values limit
#'   threading overhead for small batches. Zero is treated as one.
#' @param byrow If `TRUE`, items are rows of `X`; otherwise, they are columns.
#'   Results follow the same orientation.
#' @param random_seed Seed for hnswlib's index construction. R's `set.seed()`
#'   has no effect. Use serial construction for repeatable builds; parallel
#'   insertion order can vary with a fixed seed.
#' @return A list with matrices `idx` (one-based neighbor indices) and `dist`
#'   (distances), ordered by increasing distance for each query. Both are
#'   `n × k` when `byrow = TRUE`, or `k × n` otherwise, where `n` is the number
#'   of items in `X`.
#' @details
#' If you are searching for neighbors within `X`, don't assume that the first
#' neighbor is always the item itself: approximate search can miss self-matches.
#' With inner-product distance, an item need not be its own nearest neighbor
#' even in an exact search.
#'
#' Coordinates are stored in single precision. See [RcppHnsw-package] for
#' numeric limits. For help choosing the parameters, see the
#' [parameter guide](https://jlmelville.github.io/rcpphnsw/articles/parameters-metrics-data-layout.html)
#' which also covers data layout.
#' @seealso [hnsw_build()] and [hnsw_search()] to build once and search new data.
#' @examples
#' irism <- as.matrix(iris[, -5])
#' neighbors <- hnsw_knn(irism, k = 4)
#' neighbors$idx[1:5, ]
#' neighbors$dist[1:5, ]
#' @references
#' Malkov, Y. A., & Yashunin, D. A. (2020). Efficient and robust approximate
#' nearest neighbor search using Hierarchical Navigable Small World graphs.
#' *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 42(4),
#' 824-836. \doi{10.1109/TPAMI.2018.2889473}.
hnsw_knn <- function(
  X,
  k = 10,
  distance = "euclidean",
  M = 16,
  ef_construction = 200,
  ef = 10,
  verbose = FALSE,
  progress = "bar",
  n_threads = 0,
  grain_size = 1,
  byrow = TRUE,
  random_seed = 100
) {
  byrow <- check_logical(byrow, "byrow")
  verbose <- check_logical(verbose, "verbose")
  shape <- check_input_matrix(X, byrow)
  k <- check_whole_number(k, "k", lower = 1)
  M <- check_whole_number(M, "M", lower = 2, upper = 10000)
  ef_construction <- check_whole_number(
    ef_construction,
    "ef_construction",
    lower = 1
  )
  ef <- check_whole_number(ef, "ef", lower = 1)
  n_threads <- check_whole_number(n_threads, "n_threads", lower = 0)
  grain_size <- check_whole_number(grain_size, "grain_size", lower = 0)
  random_seed <- check_random_seed(random_seed)

  ef_construction <- max(ef_construction, k)

  max_k <- shape$nitems
  if (k > max_k) {
    stop("k cannot be larger than ", max_k, call. = FALSE)
  }
  distance <-
    match.arg(distance, c("l2", "euclidean", "cosine", "ip"))

  ann <- hnsw_build(
    X = X,
    distance = distance,
    M = M,
    ef = ef_construction,
    verbose = verbose,
    progress = progress,
    n_threads = n_threads,
    grain_size = grain_size,
    byrow = byrow,
    random_seed = random_seed
  )
  hnsw_search(
    X = X,
    ann = ann,
    k = k,
    ef = ef,
    verbose = verbose,
    progress = progress,
    n_threads = n_threads,
    grain_size = grain_size,
    byrow = byrow
  )
}

#' Build a nearest neighbor index
#'
#' Build an HNSW index that you can keep and query with [hnsw_search()].
#'
#' @param X Numeric matrix to index, with one item per row, or per column when
#'   `byrow = FALSE`.
#' @inheritParams hnsw_knn
#' @param ef Candidate-list size during construction. Larger
#'   values improve index quality and increase build time. Must be a positive
#'   whole number and is not capped at the dataset size. This is the
#'   `ef_construction` parameter of [hnsw_knn()].
#' @param byrow If `TRUE`, items are rows of `X`; otherwise, they are columns.
#' @return An `HnswEuclidean`, `HnswL2`, `HnswCosine`, or `HnswIp` index,
#'   according to `distance`. Labels are one-based and follow the order in `X`.
#' @details
#' Coordinates are stored in single precision; see [RcppHnsw-package] for
#' numeric limits. If you want to add more items, resize, or save the index,
#' see the
#' [Module guide](https://jlmelville.github.io/rcpphnsw/articles/module-api-index-lifecycle.html)
#' for the available methods.
#' @seealso [hnsw_knn()] to build and search in one call.
#' @examples
#' irism <- as.matrix(iris[, -5])
#' ann <- hnsw_build(irism[1:100, ])
#' neighbors <- hnsw_search(irism[101:150, ], ann, k = 5)
hnsw_build <- function(
  X,
  distance = "euclidean",
  M = 16,
  ef = 200,
  verbose = FALSE,
  progress = "bar",
  n_threads = 0,
  grain_size = 1,
  byrow = TRUE,
  random_seed = 100
) {
  byrow <- check_logical(byrow, "byrow")
  verbose <- check_logical(verbose, "verbose")
  shape <- check_input_matrix(X, byrow)
  M <- check_whole_number(M, "M", lower = 2, upper = 10000)
  ef <- check_whole_number(ef, "ef", lower = 1)
  n_threads <- check_whole_number(n_threads, "n_threads", lower = 0)
  grain_size <- check_whole_number(grain_size, "grain_size", lower = 0)
  seed <- check_random_seed(random_seed) # nolint: object_usage_linter.

  distance <-
    match.arg(distance, c("l2", "euclidean", "cosine", "ip"))

  nitems <- shape$nitems
  ndim <- shape$ndim
  clazz <- switch(
    distance,
    "l2" = RcppHNSW::HnswL2,
    "euclidean" = RcppHNSW::HnswEuclidean,
    "cosine" = RcppHNSW::HnswCosine,
    "ip" = RcppHNSW::HnswIp
  )
  # Create the indexing object. You must say up front the number of items that
  # will be stored (nitems).
  ann <- methods::new(clazz, ndim, nitems, M, ef, seed)

  tsmessage(
    "Building HNSW index with metric '",
    distance,
    "'",
    " ef = ",
    formatC(ef),
    " M = ",
    formatC(M),
    " using ",
    n_threads,
    " threads"
  )
  ann$setNumThreads(n_threads)
  ann$setGrainSize(grain_size)

  if (byrow) {
    ann$addItems(X)
  } else {
    ann$addItemsCol(X)
  }

  tsmessage("Finished building index")
  ann
}

#' Search a nearest neighbor index
#'
#' Find neighbors in an existing index for each query in `X`.
#'
#' @param X Numeric query matrix with the same number of dimensions as the
#'   indexed data. Each row is a query, or each column when `byrow = FALSE`.
#' @param ann An index from [hnsw_build()] or a Module constructor.
#' @param k Number of neighbors per query. A positive whole number no larger
#'   than the number of undeleted items. `ann$size()` includes deleted items.
#' @inheritParams hnsw_knn
#' @return A list with matrices `idx` (one-based labels in `ann`) and `dist`
#'   (distances), ordered by increasing distance for each query. Both are
#'   `n × k` when `byrow = TRUE`, or `k × n` otherwise, where `n` is the number
#'   of queries in `X`.
#' @details
#' The index's distance measure is used for all queries. Search is approximate;
#' with inner-product distance, an item need not be its own nearest neighbor.
#' Coordinates are stored in single precision; see [RcppHnsw-package] for numeric
#' limits.
#'
#' This call also updates `ann`'s `ef`, `n_threads`, and `grain_size` settings.
#' If you then use the Module methods directly, they will use these settings.
#' @seealso [hnsw_build()] to create an index, [hnsw_knn()] to build and search
#'   in one call.
#' @examples
#' irism <- as.matrix(iris[, -5])
#' ann <- hnsw_build(irism[1:100, ])
#' neighbors <- hnsw_search(irism[101:150, ], ann, k = 5)
hnsw_search <-
  function(
    X,
    ann,
    k,
    ef = 10,
    verbose = FALSE,
    progress = "bar",
    n_threads = 0,
    grain_size = 1,
    byrow = TRUE
  ) {
    byrow <- check_logical(byrow, "byrow")
    verbose <- check_logical(verbose, "verbose")
    check_input_matrix(X, byrow, allow_empty = TRUE)
    k <- check_whole_number(k, "k", lower = 1)
    ef <- check_whole_number(ef, "ef", lower = 1)
    n_threads <- check_whole_number(n_threads, "n_threads", lower = 0)
    grain_size <- check_whole_number(grain_size, "grain_size", lower = 0)

    max_k <- ann$size()
    if (k > max_k) {
      stop("k cannot be larger than ", max_k, call. = FALSE)
    }

    ef <- max(ef, k)

    ann$setEf(ef)
    ann$setNumThreads(n_threads)
    ann$setGrainSize(grain_size)
    tsmessage(
      "Searching HNSW index with ef = ",
      formatC(ef),
      " and ",
      n_threads,
      " threads"
    )

    if (byrow) {
      res <- ann$getAllNNsList(X, k, TRUE)
    } else {
      res <- ann$getAllNNsListCol(X, k, TRUE)
    }

    dist <- res$distance
    tsmessage("Finished searching")
    list(idx = res$item, dist = dist)
  }
