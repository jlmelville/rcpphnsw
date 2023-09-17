library(RcppHNSW)
context("get items")

ann <- new(HnswL2, ncol(ui10), nrow(ui10), M = 200, ef = 16)
# error thrown if empty index is searched
expect_error(ann$getItems(c(1, 10)), "(?i)invalid index")

ann$addItems(ui10)

expect_equivalent(ann$getItems(c(1, 10)), ui10[c(1, 10), ], tol = 1.e-7)

# error thrown if too many items are requested
expect_error(ann$getItems(c(1, 100)), "(?i)invalid index")


# repeat with cosine
ann <- new(HnswCosine, ncol(ui10), nrow(ui10), M = 200, ef = 16)
# error thrown if empty index is searched
expect_error(ann$getItems(c(1, 10)), "(?i)invalid index")

ann$addItems(ui10)

# vectors are returned l2 row-normalized in the cosine case
ui10_rl2 <-
  sweep(ui10, 1, sqrt(rowSums(apply(ui10, 2, function(x) {
    x * x
  }))), `/`)
expect_equivalent(ann$getItems(c(1, 10)), ui10_rl2[c(1, 10), ], tol = 1.e-7)

# error thrown if too many items are requested
expect_error(ann$getItems(c(1, 100)), "(?i)invalid index")
