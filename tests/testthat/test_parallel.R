library(RcppHNSW)
test_that("serial and parallel searches agree", {
  set.seed(1)

  x <- matrix(rnorm(n = 1280 * 10), ncol = 10)

  ind_0 <- hnsw_build(
    x,
    distance = "euclidean",
    M = 16,
    ef = 200,
    verbose = FALSE,
    progress = "bar",
    n_threads = 0
  )

  ind_1 <- hnsw_build(
    x,
    distance = "euclidean",
    M = 16,
    ef = 200,
    verbose = FALSE,
    progress = "bar",
    n_threads = 1
  )

  ind_2 <- hnsw_build(
    x,
    distance = "euclidean",
    M = 16,
    ef = 200,
    verbose = FALSE,
    progress = "bar",
    n_threads = 2
  )

  knn_0 <- hnsw_search(x, ind_0, k = 5, n_threads = 0)
  knn_1 <- hnsw_search(x, ind_1, k = 5, n_threads = 1)
  knn_2 <- hnsw_search(x, ind_2, k = 5, n_threads = 2)

  # Parallel construction is nondeterministic, but distances and neighbor sets
  # should remain close to the serial results.
  expect_lt(mean(abs(knn_0$dist - knn_2$dist)), 1e-3)
  expect_lt(mean(abs(knn_1$dist - knn_2$dist)), 1e-3)

  mean_overlap <- function(lhs, rhs) {
    mean(vapply(
      seq_len(nrow(lhs)),
      function(i) mean(lhs[i, ] %in% rhs[i, ]),
      numeric(1)
    ))
  }
  expect_gte(mean_overlap(knn_0$idx, knn_2$idx), 0.99)
  expect_gte(mean_overlap(knn_1$idx, knn_2$idx), 0.99)
})
