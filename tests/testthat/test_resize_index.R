library(RcppHNSW)
test_that("resize index", {
  num_elements <- nrow(uirism)
  dim <- ncol(uirism)

  M <- 16
  ef_construction <- 10
  p <- new(HnswL2, dim, num_elements, M, ef_construction)

  for (i in 1:num_elements) {
    p$addItem(uirism[i, ])
  }

  idx <- rep(0, num_elements)
  for (i in 1:num_elements) {
    idx[i] <- p$getNNs(uirism[i, ], k = 1)
  }

  recall <- mean(idx == 1:num_elements)
  expect_equal(recall, 1)

  num_elements <- nrow(uirism)
  dim <- ncol(uirism)
  p <- new(HnswL2, dim, floor(num_elements / 2), 16, 10)

  for (i in 1:(floor(num_elements / 2))) {
    p$addItem(uirism[i, ])
  }

  p$resizeIndex(num_elements)

  for (i in (floor(num_elements / 2) + 1):num_elements) {
    p$addItem(uirism[i, ])
  }

  idx <- rep(0, num_elements)
  for (i in 1:num_elements) {
    idx[i] <- p$getNNs(uirism[i, ], k = 1)
  }
  serde_recall <- mean(idx == 1:num_elements)
  expect_equal(serde_recall, recall)
})

test_that("zero capacity is rejected before resizing an empty index", {
  classes <- list(
    l2 = RcppHNSW::HnswL2,
    euclidean = RcppHNSW::HnswEuclidean,
    cosine = RcppHNSW::HnswCosine,
    ip = RcppHNSW::HnswIp
  )

  for (name in names(classes)) {
    path <- tempfile(fileext = ".hnsw")
    on.exit(unlink(path), add = TRUE)
    ann <- methods::new(classes[[name]], 2, 2, 16, 10)

    expect_error(ann$resizeIndex(0), "^new_size ", info = name)
    expect_no_error(ann$addItem(c(1, 0)))
    expect_identical(ann$getNNs(c(1, 0), 1), 1, info = name)
    expect_no_error(ann$save(path))

    rm(ann)
    expect_silent(gc())
  }
})

test_that("resize capacity includes deleted items", {
  path <- tempfile(fileext = ".hnsw")
  on.exit(unlink(path), add = TRUE)
  ann <- methods::new(RcppHNSW::HnswL2, 2, 3, 16, 10)
  ann$addItems(diag(2))
  ann$markDeleted(1)

  expect_error(ann$resizeIndex(1), "whole-number range 2 to")
  expect_identical(ann$size(), 2)
  expect_identical(ann$getNNs(c(0, 1), 1), 2)
  expect_no_error(ann$save(path))
  expect_no_error(ann$resizeIndex(2))
})
