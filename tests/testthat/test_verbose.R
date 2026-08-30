library(RcppHNSW)
test_that("verbosity levels", {
  set.seed(1337)
  data <- matrix(rnorm(100 * 10), nrow = 100)

  expected_messages <- c(
    "Building HNSW index with metric 'euclidean' ef = 200 M = 16 using 0 threads",
    "Finished building index",
    "Searching HNSW index with ef = 10 and 0 threads",
    "Finished searching"
  )
  strip_timestamp <- function(messages) {
    messages <- sub("^[0-9]{2}:[0-9]{2}:[0-9]{2} ", "", messages)
    sub("\n$", "", messages)
  }

  # results should be the same
  messages_p <- capture_messages(
    res_p <- hnsw_knn(
      data,
      k = 10,
      distance = "euclidean",
      verbose = TRUE,
      progress = "bar"
    )
  )
  messages_v <- capture_messages(
    res_v <- hnsw_knn(
      data,
      k = 10,
      distance = "euclidean",
      verbose = TRUE,
      progress = NULL
    )
  )
  messages_q <- capture_messages(
    res_q <- hnsw_knn(data, k = 10, distance = "euclidean", verbose = FALSE)
  )

  expect_equal(strip_timestamp(messages_p), expected_messages)
  expect_equal(strip_timestamp(messages_v), expected_messages)
  expect_length(messages_q, 0)
  expect_equal(sum(res_p$idx - res_v$idx), 0)
  expect_equal(sum(res_p$idx - res_q$idx), 0)
})
