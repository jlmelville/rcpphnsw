library(RcppHNSW)
test_that("hnsw_knn", {
  expect_error(hnsw_knn(uiris), "(?i)matrix")
  expect_error(hnsw_knn(ui10, M = 1), "M cannot")
  res <- hnsw_knn(ui10, k = 4)
  expect_equal(res$idx, self_nn_index4, ignore_attr = TRUE)
  expect_equal(
    res$dist,
    self_nn_dist4,
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
  expect_nn_result(res, ui10, ui10, k = 4, metric = "euclidean")

  res <- hnsw_knn(ui10, k = 4, distance = "l2")
  expect_equal(res$idx, self_nn_index4, ignore_attr = TRUE)
  expect_equal(
    res$dist,
    self_nn_dist4^2,
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
  expect_nn_result(res, ui10, ui10, k = 4, metric = "l2")

  res <- hnsw_knn(ui10, k = 1)
  expect_true(is.matrix(res$idx))
  expect_true(is.matrix(res$dist))
  expect_equal(res$idx, self_nn_index4[, 1], ignore_attr = TRUE)
  expect_equal(
    res$dist,
    self_nn_dist4[, 1],
    ignore_attr = TRUE,
    tolerance = 1e-6
  )

  # test byrow = TRUE
  res <- hnsw_knn(t(ui10), k = 4, byrow = FALSE)
  expect_equal(t(res$idx), self_nn_index4, ignore_attr = TRUE)
  expect_equal(
    t(res$dist),
    self_nn_dist4,
    ignore_attr = TRUE,
    tolerance = 1e-6
  )
})
