nn_distance <- function(query, reference, metric) {
  metric <- match.arg(metric, c("l2", "euclidean", "cosine", "ip"))

  if (metric == "l2") {
    return(sum((query - reference)^2))
  }
  if (metric == "euclidean") {
    return(sqrt(sum((query - reference)^2)))
  }
  if (metric == "cosine") {
    return(
      1 -
        sum(query * reference) /
          sqrt(sum(query^2) * sum(reference^2))
    )
  }

  1 - sum(query * reference)
}

nn_distances_at_ids <- function(query, reference, idx, metric) {
  distances <- matrix(NA_real_, nrow = nrow(idx), ncol = ncol(idx))
  for (i in seq_len(nrow(idx))) {
    for (j in seq_len(ncol(idx))) {
      distances[i, j] <- nn_distance(
        query[i, ],
        reference[idx[i, j], ],
        metric
      )
    }
  }
  distances
}

expect_equal_abs <- function(object, expected, tolerance) {
  testthat::expect_identical(dim(object), dim(expected))
  difference <- abs(as.numeric(object) - as.numeric(expected))
  testthat::expect_true(all(is.finite(difference)))
  testthat::expect_lte(max(difference, 0), tolerance)
}

expect_nn_result <- function(
  result,
  query,
  reference,
  k,
  metric,
  tolerance = 1e-6,
  sorted = TRUE
) {
  testthat::expect_named(result, c("idx", "dist"))
  testthat::expect_true(is.matrix(result$idx))
  testthat::expect_true(is.matrix(result$dist))
  expected_dim <- c(nrow(query), as.integer(k))
  testthat::expect_identical(dim(result$idx), expected_dim)
  testthat::expect_identical(dim(result$dist), expected_dim)
  testthat::expect_true(all(result$idx >= 1L & result$idx <= nrow(reference)))
  testthat::expect_true(all(apply(result$idx, 1, anyDuplicated) == 0L))

  if (sorted) {
    testthat::expect_true(all(apply(result$dist, 1, diff) >= -tolerance))
  }

  expected_distances <- nn_distances_at_ids(
    query,
    reference,
    result$idx,
    metric
  )
  expect_equal_abs(result$dist, expected_distances, tolerance)
  invisible(result)
}
