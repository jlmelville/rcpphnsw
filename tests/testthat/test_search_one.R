library(RcppHNSW)
context("search with small test set")

# #5: when searching with a pre-built index, k <= "training" set size,
# not the test set size.
index <- hnsw_build(ui10)

res <- hnsw_search(ui10[1, , drop = FALSE], index, k = 4)

expect_equal(res$idx, self_nn_index4[1, , drop = FALSE], tol = 1e-6)
expect_equal(res$dist, self_nn_dist4[1, , drop = FALSE], tol = 1e-6)
