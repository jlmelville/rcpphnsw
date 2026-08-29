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

test_that("thread counts above the item count preserve result contracts", {
  x <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 1)
  )

  for (n_threads in c(0, 1, 2, 8)) {
    ann <- hnsw_build(
      x,
      distance = "l2",
      M = 16,
      ef = 10,
      n_threads = n_threads,
      grain_size = 1
    )
    result <- hnsw_search(
      x,
      ann,
      k = 2,
      ef = 10,
      n_threads = n_threads,
      grain_size = 1
    )

    expect_nn_result(result, x, x, k = 2, metric = "l2")
    expect_identical(result$idx[, 1], seq_len(nrow(x)))
  }
})

test_that("parallel construction preserves orientation and metric contracts", {
  x <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 0),
    c(1, 0, 1),
    c(0, 1, 1)
  )

  for (metric in c("l2", "cosine")) {
    for (byrow in c(TRUE, FALSE)) {
      input <- if (byrow) x else t(x)
      ann <- hnsw_build(
        input,
        distance = metric,
        M = 4,
        ef = 10,
        n_threads = 4,
        grain_size = 1,
        random_seed = 42,
        byrow = byrow
      )
      result <- hnsw_search(
        input,
        ann,
        k = 2,
        ef = 10,
        n_threads = 2,
        grain_size = 1,
        byrow = byrow
      )
      if (!byrow) {
        result <- list(idx = t(result$idx), dist = t(result$dist))
      }

      expect_nn_result(result, x, x, k = 2, metric = metric)
      expect_identical(result$idx[, 1], seq_len(nrow(x)))
    }
  }
})
