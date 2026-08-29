library(RcppHNSW)

test_that("insertion preflight failures leave the index usable", {
  items <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 0)
  )
  ann <- methods::new(RcppHNSW::HnswL2, 3, nrow(items), 16, 10)
  ann$setNumThreads(2)

  expect_error(ann$setNumThreads(-1), "^n_threads ")
  expect_error(ann$addItem(c(1, 0)), "incorrect dimensions")
  expect_error(
    ann$addItems(matrix(1, nrow = 2, ncol = 2)),
    "incorrect dimensions"
  )
  invalid <- items[1:2, , drop = FALSE]
  invalid[2, 1] <- Inf
  expect_error(ann$addItems(invalid), "representable")
  expect_identical(ann$size(), 0)

  ann$addItems(items[1:2, , drop = FALSE])
  expect_error(
    ann$addItems(items[2:4, , drop = FALSE]),
    "Index is too small"
  )
  expect_identical(ann$size(), 2)

  ann$addItems(items[3:4, , drop = FALSE])
  expect_identical(ann$size(), 4)
  expect_equal(ann$getItems(1:4), items, tolerance = 1e-6)
})

test_that("cosine batch preflight completes before insertion starts", {
  items <- rbind(c(1, 0, 0), c(0, 1, 0))
  ann <- methods::new(RcppHNSW::HnswCosine, 3, 2, 16, 10)
  ann$setNumThreads(2)

  invalid <- items
  invalid[2, ] <- 0
  expect_error(ann$addItems(invalid), "positive finite norm")
  expect_identical(ann$size(), 0)

  expect_no_error(ann$addItems(items))
  expect_identical(ann$size(), 2)
})

test_that("non-mutating failures leave index operations available", {
  items <- diag(3)
  ann <- hnsw_build(items, distance = "l2", n_threads = 2)

  expect_error(
    ann$getAllNNs(matrix(1, nrow = 1, ncol = 2), 1),
    "incorrect dimensions"
  )
  expect_error(ann$getNNs(items[1, ], 4), "active item count")
  expect_error(ann$setEf(0), "^ef ")
  expect_error(ann$save(""), "^path_to_index ")

  expect_identical(ann$size(), 3)
  expect_no_error(ann$setEf(10))
  expect_no_error(ann$setNumThreads(2))
  expect_no_error(ann$setGrainSize(1))
  expect_equal(ann$getNNs(items[1, ], 1), 1)
  expect_equal(ann$getItems(1:3), items, tolerance = 1e-6)
})
