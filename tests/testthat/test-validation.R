library(RcppHNSW)

test_that("high-level whole-number controls reject invalid values", {
  x <- matrix(seq_len(12), nrow = 4)
  ann <- hnsw_build(x)

  controls <- list(
    list(
      api = "hnsw_knn",
      name = "k",
      lower = 1,
      upper = nrow(x),
      call = function(value) hnsw_knn(x, k = value)
    ),
    list(
      api = "hnsw_knn",
      name = "M",
      lower = 2,
      upper = 10000,
      call = function(value) hnsw_knn(x, k = 1, M = value)
    ),
    list(
      api = "hnsw_knn",
      name = "ef_construction",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) {
        hnsw_knn(x, k = 1, ef_construction = value)
      }
    ),
    list(
      api = "hnsw_knn",
      name = "ef",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) hnsw_knn(x, k = 1, ef = value)
    ),
    list(
      api = "hnsw_knn",
      name = "n_threads",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_knn(x, k = 1, n_threads = value)
    ),
    list(
      api = "hnsw_knn",
      name = "grain_size",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_knn(x, k = 1, grain_size = value)
    ),
    list(
      api = "hnsw_knn",
      name = "random_seed",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_knn(x, k = 1, random_seed = value)
    ),
    list(
      api = "hnsw_build",
      name = "M",
      lower = 2,
      upper = 10000,
      call = function(value) hnsw_build(x, M = value)
    ),
    list(
      api = "hnsw_build",
      name = "ef",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) hnsw_build(x, ef = value)
    ),
    list(
      api = "hnsw_build",
      name = "n_threads",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_build(x, n_threads = value)
    ),
    list(
      api = "hnsw_build",
      name = "grain_size",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_build(x, grain_size = value)
    ),
    list(
      api = "hnsw_build",
      name = "random_seed",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_build(x, random_seed = value)
    ),
    list(
      api = "hnsw_search",
      name = "k",
      lower = 1,
      upper = ann$size(),
      call = function(value) hnsw_search(x, ann, k = value)
    ),
    list(
      api = "hnsw_search",
      name = "ef",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) hnsw_search(x, ann, k = 1, ef = value)
    ),
    list(
      api = "hnsw_search",
      name = "n_threads",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_search(x, ann, k = 1, n_threads = value)
    ),
    list(
      api = "hnsw_search",
      name = "grain_size",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) hnsw_search(x, ann, k = 1, grain_size = value)
    )
  )

  invalid_values <- function(lower, upper) {
    values <- list(
      character = "1",
      logical = TRUE,
      complex = 1 + 0i,
      null = NULL,
      empty = numeric(),
      length_two = c(lower, lower),
      missing = NA_real_,
      not_a_number = NaN,
      positive_infinity = Inf,
      negative_infinity = -Inf,
      fractional = lower + 0.5,
      below_range = lower - 1,
      above_range = as.double(upper) + 1,
      above_integer_max = as.double(.Machine$integer.max) + 1
    )
    if (lower > 0) {
      values$zero <- 0
    }
    values
  }

  for (control in controls) {
    invalid <- invalid_values(control$lower, control$upper)
    for (i in seq_along(invalid)) {
      warnings <- character()
      result <- withCallingHandlers(
        tryCatch(control$call(invalid[[i]]), error = function(e) e),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      label <- paste(control$api, control$name, names(invalid)[i])
      is_error <- inherits(result, "error")
      error_message <- if (is_error) conditionMessage(result) else ""

      expect_true(is_error, info = label)
      expect_true(
        grepl(paste0("^", control$name, " "), error_message),
        info = label
      )
      expect_true(length(warnings) == 0, info = label)
    }
  }
})

test_that("high-level logical controls require one non-missing flag", {
  x <- matrix(seq_len(12), nrow = 4)
  ann <- hnsw_build(x)
  controls <- list(
    list(
      api = "hnsw_knn",
      name = "byrow",
      call = function(value) hnsw_knn(x, k = 1, byrow = value)
    ),
    list(
      api = "hnsw_knn",
      name = "verbose",
      call = function(value) hnsw_knn(x, k = 1, verbose = value)
    ),
    list(
      api = "hnsw_build",
      name = "byrow",
      call = function(value) hnsw_build(x, byrow = value)
    ),
    list(
      api = "hnsw_build",
      name = "verbose",
      call = function(value) hnsw_build(x, verbose = value)
    ),
    list(
      api = "hnsw_search",
      name = "byrow",
      call = function(value) hnsw_search(x, ann, k = 1, byrow = value)
    ),
    list(
      api = "hnsw_search",
      name = "verbose",
      call = function(value) hnsw_search(x, ann, k = 1, verbose = value)
    )
  )
  invalid <- list(
    null = NULL,
    empty = logical(),
    length_two = c(TRUE, FALSE),
    missing = NA,
    numeric = 1,
    character = "TRUE"
  )

  for (control in controls) {
    for (i in seq_along(invalid)) {
      result <- tryCatch(control$call(invalid[[i]]), error = function(e) e)
      label <- paste(control$api, control$name, names(invalid)[i])
      is_error <- inherits(result, "error")
      error_message <- if (is_error) conditionMessage(result) else ""

      expect_true(is_error, info = label)
      expect_true(
        grepl(paste0("^", control$name, " "), error_message),
        info = label
      )
    }
  }
})

test_that("whole-number boundary values remain accepted", {
  x <- matrix(seq_len(12), nrow = 4)

  expect_warning(
    hnsw_knn(
      x,
      k = 1,
      M = 2,
      ef_construction = 1,
      ef = 1,
      n_threads = 0,
      grain_size = 0,
      random_seed = 0
    ),
    NA
  )

  ann <- expect_warning(
    hnsw_build(
      x,
      M = 10000,
      ef = .Machine$integer.max,
      n_threads = .Machine$integer.max,
      grain_size = .Machine$integer.max,
      random_seed = .Machine$integer.max
    ),
    NA
  )
  empty <- matrix(numeric(), nrow = 0, ncol = ncol(x))
  expect_warning(
    hnsw_search(
      empty,
      ann,
      k = ann$size(),
      ef = .Machine$integer.max,
      n_threads = .Machine$integer.max,
      grain_size = .Machine$integer.max
    ),
    NA
  )
})

test_that("rectangular high-level calls respect item orientation", {
  # fmt: skip
  reference <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    0, 2, 0,
    0, 0, 3
  ), nrow = 4, byrow = TRUE)

  row_knn <- hnsw_knn(reference, k = 4)
  expect_nn_result(row_knn, reference, reference, 4, "euclidean")

  column_knn <- hnsw_knn(t(reference), k = 4, byrow = FALSE)
  expect_nn_result(
    list(idx = t(column_knn$idx), dist = t(column_knn$dist)),
    reference,
    reference,
    4,
    "euclidean"
  )

  expect_error(
    hnsw_knn(reference, k = 4, byrow = FALSE),
    "^k cannot be larger than 3$"
  )

  row_index <- hnsw_build(reference)
  row_query <- reference[c(1, 3), , drop = FALSE]
  row_search <- hnsw_search(row_query, row_index, k = 2)
  expect_nn_result(row_search, row_query, reference, 2, "euclidean")

  column_index <- hnsw_build(t(reference), byrow = FALSE)
  column_search <- hnsw_search(
    t(row_query),
    column_index,
    k = 2,
    byrow = FALSE
  )
  expect_nn_result(
    list(idx = t(column_search$idx), dist = t(column_search$dist)),
    row_query,
    reference,
    2,
    "euclidean"
  )
})

test_that("high-level construction rejects empty items and dimensions", {
  cases <- list(
    list(
      name = "row empty items",
      X = matrix(numeric(), nrow = 0, ncol = 3),
      byrow = TRUE,
      message = "item"
    ),
    list(
      name = "row empty dimensions",
      X = matrix(numeric(), nrow = 3, ncol = 0),
      byrow = TRUE,
      message = "dimension"
    ),
    list(
      name = "column empty items",
      X = matrix(numeric(), nrow = 3, ncol = 0),
      byrow = FALSE,
      message = "item"
    ),
    list(
      name = "column empty dimensions",
      X = matrix(numeric(), nrow = 0, ncol = 3),
      byrow = FALSE,
      message = "dimension"
    )
  )

  for (case in cases) {
    calls <- list(
      hnsw_build = function() hnsw_build(case$X, byrow = case$byrow),
      hnsw_knn = function() hnsw_knn(case$X, k = 1, byrow = case$byrow)
    )
    for (name in names(calls)) {
      result <- tryCatch(calls[[name]](), error = function(e) e)
      label <- paste(name, case$name)
      is_error <- inherits(result, "error")
      error_message <- if (is_error) conditionMessage(result) else ""

      expect_true(is_error, info = label)
      expect_true(grepl(case$message, error_message), info = label)
    }
  }
})

test_that("empty high-level queries preserve orientation", {
  reference <- matrix(seq_len(12), nrow = 4)
  row_index <- hnsw_build(reference)
  row_empty <- matrix(numeric(), nrow = 0, ncol = ncol(reference))
  row_result <- hnsw_search(row_empty, row_index, k = 2)
  expect_identical(dim(row_result$idx), c(0L, 2L))
  expect_identical(dim(row_result$dist), c(0L, 2L))

  column_reference <- t(reference)
  column_index <- hnsw_build(column_reference, byrow = FALSE)
  column_empty <- matrix(
    numeric(),
    nrow = nrow(column_reference),
    ncol = 0
  )
  column_result <- hnsw_search(
    column_empty,
    column_index,
    k = 2,
    byrow = FALSE
  )
  expect_identical(dim(column_result$idx), c(2L, 0L))
  expect_identical(dim(column_result$dist), c(2L, 0L))

  expect_error(
    hnsw_search(matrix(numeric(), nrow = 2, ncol = 0), row_index, k = 1),
    "dimension"
  )
  expect_error(
    hnsw_search(
      matrix(numeric(), nrow = 0, ncol = 2),
      column_index,
      k = 1,
      byrow = FALSE
    ),
    "dimension"
  )
})
