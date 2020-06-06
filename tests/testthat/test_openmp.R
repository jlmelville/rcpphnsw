library(RcppHNSW)
context("test OpenMP")

set.seed(1)

x <- matrix(rnorm(n = 1280 * 10), ncol = 10)

ind_1 <- hnsw_build(
  x,
  distance = "euclidean",
  M = 16,
  ef = 200,
  verbose = FALSE,
  progress = "bar",
  n_threads = 1
)

ind_4 <- hnsw_build(
  x,
  distance = "euclidean",
  M = 16,
  ef = 200,
  verbose = FALSE,
  progress = "bar",
  n_threads = 4
)

knn_1 = hnsw_search(x, ind_1, k = 5)
knn_4 = hnsw_search(x, ind_4, k = 5)
# Seems index which was built using more than 1 thread is not determenistic
# but in general the difference between index built with 1 thread and
# many threads should be small
expect_lt(mean(knn_1$dist - knn_4$dist), 1e-4)

# same check for indices
expect_lt(mean(knn_1$idx != knn_4$idx), 1e-2)
