library(RcppHNSW)
context("search with small test set")

# #5: when searching with a pre-built index, k <= "training" set size,
# not the test set size.
index <- hnsw_build(ui10)
# Added a size method so the user can check the index size if needed
expect_equal(index$size(), 10)

res <- hnsw_search(ui10[1, , drop = FALSE], index, k = 4)

expect_equal(res$idx, self_nn_index4[1, , drop = FALSE], tol = 1e-6)
expect_equal(res$dist, self_nn_dist4[1, , drop = FALSE], tol = 1e-6)

# An error will be thrown if not enough neighbors can be retrieved
# in this case due to k > index$size()
expect_error(hnsw_search(ui10[1, , drop = FALSE], index, k = 500))
