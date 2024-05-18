library(RcppHNSW)
context("verbosity levels")

set.seed(1337)
data <- matrix(rnorm(100 * 10), nrow = 100)

# results should be the same
# 23: also counts as a test of setting the seed to get repeatable results
set.seed(1337)
res_p <- hnsw_knn(data, k = 10, distance = "euclidean", verbose = TRUE, progress = "bar")
set.seed(1337)
res_v <- hnsw_knn(data, k = 10, distance = "euclidean", verbose = TRUE, progress = NULL)
set.seed(1337)
res_q <- hnsw_knn(data, k = 10, distance = "euclidean", verbose = FALSE)
expect_equal(sum(res_p$idx - res_v$idx), 0)
expect_equal(sum(res_p$idx - res_q$idx), 0)
