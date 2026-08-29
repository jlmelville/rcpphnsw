library(RcppHNSW)

test_that("high-level nearest-neighbor results implement every metric", {
  # fmt: skip
  x <- matrix(c(
     1, 0, 0,
     0, 2, 0,
     2, 1, 1,
    -1, 1, 3
  ), nrow = 4, byrow = TRUE)

  for (metric in c("l2", "euclidean", "cosine", "ip")) {
    result <- hnsw_knn(x, k = nrow(x), distance = metric)
    expect_nn_result(
      result,
      x,
      x,
      k = nrow(x),
      metric = metric,
      tolerance = 1e-5
    )
  }
})

test_that("inner product does not imply that self is nearest", {
  x <- rbind(
    c(1, 0),
    c(3, 0),
    c(0, 2)
  )

  result <- hnsw_knn(x, k = 1, distance = "ip")
  expect_nn_result(result, x, x, k = 1, metric = "ip")
  expect_identical(result$idx[1, 1], 2L)
  expect_lt(result$dist[1, 1], nn_distance(x[1, ], x[1, ], "ip"))
})
