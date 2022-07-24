library(RcppHNSW)
context("add/get matrix")

ann <- new(HnswL2, ncol(ui10), nrow(ui10), M = 200, ef = 16)

# error thrown if empty index is searched
expect_error(ann$getAllNNs(ui10, k = 4), "(?i)unable to find")

ann$addItems(ui10)
res <- ann$getAllNNsList(ui10, k = 4, include_distances = TRUE)
expect_equal(res$item, self_nn_index4, check.attributes = FALSE)
expect_equal(sqrt(res$dist), self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

res <- ann$getAllNNsList(ui10, k = 4, include_distances = FALSE)
expect_equal(res$item, self_nn_index4, check.attributes = FALSE)

items <- ann$getAllNNs(ui10, k = 4)
expect_equal(items, self_nn_index4, check.attributes = FALSE)
