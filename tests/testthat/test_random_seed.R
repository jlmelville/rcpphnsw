library(RcppHNSW)
context("random seeds")

set.seed(42)
expected_next <- stats::runif(1)
set.seed(42)
invisible(hnsw_build(ui10))
observed_next <- stats::runif(1)
expect_identical(observed_next, expected_next)

default_res <- hnsw_knn(ui10, k = 4)
seeded_res <- hnsw_knn(ui10, k = 4, random_seed = 100)
expect_equal(seeded_res, default_res)

expect_error(hnsw_build(ui10, random_seed = -1), "random_seed")
expect_error(hnsw_knn(ui10, k = 4, random_seed = 1.5), "random_seed")
expect_error(hnsw_build(ui10, random_seed = NULL), "random_seed")
