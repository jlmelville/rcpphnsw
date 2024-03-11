library(RcppHNSW)
context("Save/load index")

num_elements <- nrow(uirism)
dim <- ncol(uirism)

M <- 16
ef_construction <- 10
p <- new(HnswL2, dim, num_elements, M, ef_construction)

for (i in 1:num_elements) {
  p$addItem(uirism[i, ])
}

nn4idx <- matrix(0L, nrow = num_elements, ncol = 4)
nn4dist <- matrix(0.0, nrow = num_elements, ncol = 4)

for (i in 1:num_elements) {
  res <- p$getNNsList(uirism[i, ], k = 4, TRUE)
  nn4idx[i, ] <- res$item
  nn4dist[i, ] <- res$distance
}

temp_file <- tempfile()
on.exit(unlink(temp_file), add = TRUE)
p$save(temp_file)

nn4idx_aftersave <- matrix(0L, nrow = num_elements, ncol = 4)
nn4dist_aftersave <- matrix(0.0, nrow = num_elements, ncol = 4)
for (i in 1:num_elements) {
  res_aftersave <- p$getNNsList(uirism[i, ], k = 4, TRUE)
  nn4idx_aftersave[i, ] <- res_aftersave$item
  nn4dist_aftersave[i, ] <- res_aftersave$distance
}
expect_equal(nn4idx, nn4idx_aftersave)
expect_equal(nn4dist, nn4dist_aftersave)

pload <- new(HnswL2, dim, temp_file)
nn4idx_afterload <- matrix(0L, nrow = num_elements, ncol = 4)
nn4dist_afterload <- matrix(0.0, nrow = num_elements, ncol = 4)
for (i in 1:num_elements) {
  res_afterload <- pload$getNNsList(uirism[i, ], k = 4, TRUE)
  nn4idx_afterload[i, ] <- res_afterload$item
  nn4dist_afterload[i, ] <- res_afterload$distance
}
expect_equal(nn4idx, nn4idx_afterload)
expect_equal(nn4dist, nn4dist_afterload)

# 21: no way to use hnsw_search to get Euclidean distances after save/load
test_that("euclidean search is more consistent with save/load", {
  ann <- hnsw_build(ui10, distance = "euclidean")
  iris_nn <- hnsw_search(ui10, ann, k = 4)
  expect_equal(iris_nn$dist, self_nn_dist4, tol = 1e-6)
  expect_equal(iris_nn$idx, self_nn_index4)

  temp_file <- tempfile()
  on.exit(unlink(temp_file), add = TRUE)
  ann$save(temp_file)

  ann2 <- methods::new(RcppHNSW::HnswEuclidean, 4, temp_file)
  iris_nn2 <- hnsw_search(ui10, ann, k = 4)
  expect_equal(iris_nn2$dist, self_nn_dist4, tol = 1e-6)
  expect_equal(iris_nn2$idx, self_nn_index4)
})
