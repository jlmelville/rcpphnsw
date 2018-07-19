library(RcppHNSW)
context("get_knn")

res <- get_knn(ui10, k = 4, include_self = FALSE)
expect_equal(res$idx, nn_index4, check.attributes = FALSE)
expect_equal(res$dist, nn_dist4, check.attributes = FALSE, tol = 1e-6)

res <- get_knn(ui10, k = 4, distance = "l2", include_self = FALSE)
expect_equal(res$idx, nn_index4, check.attributes = FALSE)
expect_equal(res$dist, nn_dist4 ^ 2, check.attributes = FALSE, tol = 1e-6)


res <- get_knn(ui10, k = 4, include_self = TRUE)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)


