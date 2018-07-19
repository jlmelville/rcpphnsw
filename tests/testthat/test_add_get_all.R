library(RcppHNSW)
context("add/get matrix")

RcppParallel::setThreadOptions(numThreads = 1)

ann <- new(HnswL2, ncol(ui10), nrow(ui10), M = 200, ef = 16)
ann$addItems(ui10)
res <- ann$getAllNNsList(ui10, k = 4, include_distances = TRUE)

expect_equal(res$item + 1, self_nn_index4, check.attributes = FALSE)
expect_equal(sqrt(res$dist), self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

items <- ann$getAllNNs(ui10, k = 4)
expect_equal(items + 1, self_nn_index4, check.attributes = FALSE)
