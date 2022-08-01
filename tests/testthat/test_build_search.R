library(RcppHNSW)
context("build and search functions")

expect_error(hnsw_build(uiris), "(?i)matrix")
expect_error(hnsw_build(ui10, M = 1), "M cannot")
index <- hnsw_build(ui10)

expect_error(hnsw_search(uiris, index), "(?i)matrix")
res <- hnsw_search(ui10, index, k = 4)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

# test byrow = FALSE
ui10_col <- t(ui10)

index2 <- hnsw_build(ui10_col, byrow = FALSE)
res_col <- hnsw_search(ui10_col, index2, k = 4, byrow = FALSE)
expect_equal(t(res_col$idx), self_nn_index4, check.attributes = FALSE)
expect_equal(t(res_col$dist), self_nn_dist4,
  check.attributes = FALSE,
  tol = 1e-6
)
