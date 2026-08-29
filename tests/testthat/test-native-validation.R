library(RcppHNSW)

test_that("Module whole-number adapters reject original invalid R values", {
  clazz <- RcppHNSW::HnswL2
  ann <- methods::new(clazz, 3, 3, 16, 10, 100)
  ann$addItems(diag(3))
  query <- c(1, 0, 0)

  controls <- list(
    list(
      name = "dimension",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) methods::new(clazz, value, 1, 16, 10)
    ),
    list(
      name = "capacity",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) methods::new(clazz, 3, value, 16, 10)
    ),
    list(
      name = "M",
      lower = 2,
      upper = 10000,
      call = function(value) methods::new(clazz, 3, 1, value, 10)
    ),
    list(
      name = "ef_construction",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) methods::new(clazz, 3, 1, 16, value)
    ),
    list(
      name = "random_seed",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) methods::new(clazz, 3, 1, 16, 10, value)
    ),
    list(
      name = "ef",
      lower = 1,
      upper = .Machine$integer.max,
      call = function(value) ann$setEf(value)
    ),
    list(
      name = "n_threads",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) ann$setNumThreads(value)
    ),
    list(
      name = "grain_size",
      lower = 0,
      upper = .Machine$integer.max,
      call = function(value) ann$setGrainSize(value)
    ),
    list(
      name = "k",
      lower = 1,
      upper = ann$size(),
      call = function(value) ann$getNNs(query, value)
    ),
    list(
      name = "label",
      lower = 1,
      upper = ann$size(),
      call = function(value) ann$markDeleted(value)
    ),
    list(
      name = "new_size",
      lower = ann$size(),
      upper = .Machine$integer.max,
      call = function(value) ann$resizeIndex(value)
    )
  )

  invalid_values <- function(lower, upper) {
    list(
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
      label <- paste(control$name, names(invalid)[i])
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

  expect_no_error(methods::new(clazz, 3, 1, 2, 1, 0))
  expect_no_error(methods::new(clazz, 3, 1, 10000, .Machine$integer.max))
  expect_no_error(methods::new(clazz, 3, 1, 16, 1, .Machine$integer.max))
  expect_no_error(ann$setEf(.Machine$integer.max))
  expect_no_error(ann$setNumThreads(0))
  expect_no_error(ann$setNumThreads(.Machine$integer.max))
  expect_no_error(ann$setGrainSize(0))
  expect_no_error(ann$setGrainSize(.Machine$integer.max))
  expect_no_error(ann$resizeIndex(ann$size()))
})

test_that("every Module search method validates k before narrowing", {
  ann <- methods::new(RcppHNSW::HnswL2, 3, 2, 16, 10)
  items <- diag(1, nrow = 2, ncol = 3)
  ann$addItems(items)
  query <- items[1, ]
  calls <- list(
    getNNs = function() ann$getNNs(query, 1.5),
    getNNsList = function() ann$getNNsList(query, 1.5, TRUE),
    getAllNNs = function() ann$getAllNNs(items, 1.5),
    getAllNNsList = function() ann$getAllNNsList(items, 1.5, TRUE),
    getAllNNsCol = function() ann$getAllNNsCol(t(items), 1.5),
    getAllNNsListCol = function() ann$getAllNNsListCol(t(items), 1.5, TRUE)
  )

  for (name in names(calls)) {
    expect_error(calls[[name]](), "^k ", info = name)
  }
})

test_that("Module logical adapters require one non-missing flag", {
  ann <- methods::new(RcppHNSW::HnswL2, 3, 2, 16, 10)
  items <- diag(1, nrow = 2, ncol = 3)
  ann$addItems(items)
  invalid <- list(
    null = NULL,
    empty = logical(),
    length_two = c(TRUE, FALSE),
    missing = NA,
    numeric = 1,
    character = "TRUE"
  )
  calls <- list(
    getNNsList = function(value) ann$getNNsList(items[1, ], 1, value),
    getAllNNsList = function(value) ann$getAllNNsList(items, 1, value),
    getAllNNsListCol = function(value) {
      ann$getAllNNsListCol(t(items), 1, value)
    }
  )

  for (name in names(calls)) {
    for (i in seq_along(invalid)) {
      expect_error(
        calls[[name]](invalid[[i]]),
        "^include_distances ",
        info = paste(name, names(invalid)[i])
      )
    }
  }

  expect_named(ann$getNNsList(items[1, ], 1, FALSE), "item")
  expect_named(
    ann$getNNsList(items[1, ], 1, TRUE),
    c("item", "distance")
  )
})

test_that("Module path adapters reject invalid R strings", {
  clazz <- RcppHNSW::HnswL2
  ann <- methods::new(clazz, 3, 1, 16, 10)
  invalid <- list(
    null = NULL,
    empty_vector = character(),
    length_two = c("a", "b"),
    missing = NA_character_,
    empty_string = "",
    numeric = 1
  )

  for (i in seq_along(invalid)) {
    expect_error(
      methods::new(clazz, 3, invalid[[i]]),
      "^path_to_index ",
      info = paste("load", names(invalid)[i])
    )
    expect_error(
      ann$save(invalid[[i]]),
      "^path_to_index ",
      info = paste("save", names(invalid)[i])
    )
  }

  expect_error(
    methods::new(clazz, 0, tempfile()),
    "^dimension ",
    fixed = FALSE
  )
  expect_error(
    methods::new(clazz, 3, tempfile(), 1.5),
    "^capacity ",
    fixed = FALSE
  )
})

test_that("all Module data paths enforce exact vector and matrix shapes", {
  ann <- methods::new(RcppHNSW::HnswL2, 3, 4, 16, 10)
  items <- diag(2, 3)
  ann$addItems(items)
  wrong_row <- matrix(1, nrow = 2, ncol = 2)
  wrong_col <- t(wrong_row)

  calls <- list(
    addItem = function() ann$addItem(c(1, 2)),
    addItems = function() ann$addItems(wrong_row),
    addItemsCol = function() ann$addItemsCol(wrong_col),
    getNNs = function() ann$getNNs(c(1, 2), 1),
    getNNsList = function() ann$getNNsList(c(1, 2), 1, TRUE),
    getAllNNs = function() ann$getAllNNs(wrong_row, 1),
    getAllNNsList = function() ann$getAllNNsList(wrong_row, 1, TRUE),
    getAllNNsCol = function() ann$getAllNNsCol(wrong_col, 1),
    getAllNNsListCol = function() {
      ann$getAllNNsListCol(wrong_col, 1, TRUE)
    }
  )

  for (name in names(calls)) {
    expect_error(calls[[name]](), "incorrect dimensions", info = name)
  }

  dimensional_item <- array(c(1, 0, 0), dim = c(1, 1, 3))
  expect_error(ann$addItem(dimensional_item), "incorrect dimensions")
  expect_error(ann$getNNs(dimensional_item, 1), "incorrect dimensions")
  expect_error(
    ann$getNNsList(dimensional_item, 1, TRUE),
    "incorrect dimensions"
  )
  expect_error(ann$getItems(array(1, dim = c(1, 1, 1))), "numeric vector")
})

test_that("all Module data paths reject invalid float coordinates", {
  ann <- methods::new(RcppHNSW::HnswL2, 3, 4, 16, 10)
  valid <- diag(2, 3)
  ann$addItems(valid)
  invalid_vector <- c(1, Inf, 0)
  invalid_rows <- rbind(c(1, 0, 0), invalid_vector)
  invalid_columns <- t(invalid_rows)

  add_ann <- methods::new(RcppHNSW::HnswL2, 3, 4, 16, 10)
  calls <- list(
    addItem = function() add_ann$addItem(invalid_vector),
    addItems = function() add_ann$addItems(invalid_rows),
    addItemsCol = function() add_ann$addItemsCol(invalid_columns),
    getNNs = function() ann$getNNs(invalid_vector, 1),
    getNNsList = function() ann$getNNsList(invalid_vector, 1, TRUE),
    getAllNNs = function() ann$getAllNNs(invalid_rows, 1),
    getAllNNsList = function() ann$getAllNNsList(invalid_rows, 1, TRUE),
    getAllNNsCol = function() ann$getAllNNsCol(invalid_columns, 1),
    getAllNNsListCol = function() {
      ann$getAllNNsListCol(invalid_columns, 1, TRUE)
    }
  )

  for (name in names(calls)) {
    expect_error(calls[[name]](), "representable", info = name)
  }
  expect_identical(add_ann$size(), 0)

  overflow <- c(.Machine$double.xmax, 0, 0)
  expect_error(ann$getNNs(overflow, 1), "representable")
  expect_error(
    ann$getAllNNs(matrix(overflow, nrow = 1), 1),
    "representable"
  )
})

test_that("high-level paths reject non-finite and float-overflow values", {
  valid <- rbind(c(1, 0, 0), c(0, 1, 0))
  ann <- hnsw_build(valid, distance = "l2")
  invalid <- list(
    missing = NA_real_,
    not_a_number = NaN,
    positive_infinity = Inf,
    negative_infinity = -Inf,
    float_overflow = .Machine$double.xmax
  )

  for (i in seq_along(invalid)) {
    x <- valid
    x[1, 1] <- invalid[[i]]
    label <- names(invalid)[i]
    expect_error(hnsw_build(x, distance = "l2"), "representable", info = label)
    expect_error(
      hnsw_knn(x, k = 1, distance = "l2"),
      "representable",
      info = label
    )
    expect_error(hnsw_search(x, ann, k = 1), "representable", info = label)
  }
})

test_that("cosine paths reject zero vectors after float conversion", {
  valid <- rbind(c(1, 0, 0), c(0, 1, 0))
  ann <- hnsw_build(valid, distance = "cosine")
  invalid <- list(
    zero = c(0, 0, 0),
    float_underflow = rep(.Machine$double.xmin, 3)
  )

  for (i in seq_along(invalid)) {
    item <- invalid[[i]]
    rows <- rbind(valid[1, ], item)
    columns <- t(rows)
    label <- names(invalid)[i]

    expect_error(
      hnsw_build(rows, distance = "cosine"),
      "positive finite norm",
      info = paste("hnsw_build", label)
    )
    expect_error(
      hnsw_knn(rows, k = 1, distance = "cosine"),
      "positive finite norm",
      info = paste("hnsw_knn", label)
    )
    expect_error(
      hnsw_search(matrix(item, nrow = 1), ann, k = 1),
      "positive finite norm",
      info = paste("hnsw_search", label)
    )

    add_ann <- methods::new(RcppHNSW::HnswCosine, 3, 2, 16, 10)
    expect_error(add_ann$addItem(item), "positive finite norm", info = label)
    expect_error(add_ann$addItems(rows), "positive finite norm", info = label)
    expect_error(
      add_ann$addItemsCol(columns),
      "positive finite norm",
      info = label
    )
    expect_identical(add_ann$size(), 0)

    expect_error(ann$getNNs(item, 1), "positive finite norm", info = label)
    expect_error(
      ann$getNNsList(item, 1, TRUE),
      "positive finite norm",
      info = label
    )
    expect_error(
      ann$getAllNNs(rows, 1),
      "positive finite norm",
      info = label
    )
    expect_error(
      ann$getAllNNsList(rows, 1, TRUE),
      "positive finite norm",
      info = label
    )
    expect_error(
      ann$getAllNNsCol(columns, 1),
      "positive finite norm",
      info = label
    )
    expect_error(
      ann$getAllNNsListCol(columns, 1, TRUE),
      "positive finite norm",
      info = label
    )
  }
})

test_that("search k uses active rather than historical item count", {
  items <- diag(3)
  ann <- hnsw_build(items, distance = "l2")
  ann$markDeleted(1)

  expect_error(hnsw_search(items[2, , drop = FALSE], ann, k = 3), "active")
  expect_error(ann$getNNs(items[2, ], 3), "active")
  expect_error(ann$getAllNNs(items[2, , drop = FALSE], 3), "active")
  expect_error(
    ann$getAllNNsListCol(t(items[2, , drop = FALSE]), 3, TRUE),
    "active"
  )
  expect_length(ann$getNNs(items[2, ], 2), 2)
})

test_that("getItems validates identifiers before subtraction", {
  items <- diag(3)
  ann <- methods::new(RcppHNSW::HnswL2, 3, 3, 16, 10)
  ann$addItems(items)
  invalid <- list(
    zero = 0,
    negative = -1,
    missing_integer = NA_integer_,
    missing_double = NA_real_,
    not_a_number = NaN,
    infinity = Inf,
    fractional = 1.5,
    too_large = 4,
    logical = TRUE,
    character = "1",
    complex = 1 + 0i
  )

  for (i in seq_along(invalid)) {
    expect_error(ann$getItems(invalid[[i]]), info = names(invalid)[i])
  }

  duplicated <- ann$getItems(c(2, 2))
  expect_equal(duplicated[1, ], duplicated[2, ])
  expect_identical(dim(ann$getItems(integer())), c(0L, 3L))

  ann$markDeleted(1)
  expect_error(ann$getItems(1), "Label not found")
})

test_that("positive-capacity empty indexes and empty queries remain valid", {
  empty_ann <- methods::new(RcppHNSW::HnswL2, 3, 2, 16, 10)
  expect_identical(empty_ann$size(), 0)
  expect_no_error(empty_ann$addItems(matrix(numeric(), nrow = 0, ncol = 3)))
  expect_no_error(
    empty_ann$addItemsCol(matrix(numeric(), nrow = 3, ncol = 0))
  )
  expect_identical(empty_ann$size(), 0)

  ann <- methods::new(RcppHNSW::HnswL2, 3, 2, 16, 10)
  ann$addItems(diag(1, nrow = 2, ncol = 3))
  row_empty <- matrix(numeric(), nrow = 0, ncol = 3)
  col_empty <- matrix(numeric(), nrow = 3, ncol = 0)

  expect_identical(dim(ann$getAllNNs(row_empty, 1)), c(0L, 1L))
  expect_identical(dim(ann$getAllNNsList(row_empty, 1, TRUE)$item), c(0L, 1L))
  expect_identical(
    dim(ann$getAllNNsList(row_empty, 1, TRUE)$distance),
    c(0L, 1L)
  )
  expect_identical(dim(ann$getAllNNsCol(col_empty, 1)), c(1L, 0L))
  expect_identical(
    dim(ann$getAllNNsListCol(col_empty, 1, TRUE)$item),
    c(1L, 0L)
  )
  expect_identical(
    dim(ann$getAllNNsListCol(col_empty, 1, TRUE)$distance),
    c(1L, 0L)
  )

  expect_error(methods::new(RcppHNSW::HnswL2, 0, 1, 16, 10), "dimension")
  expect_error(methods::new(RcppHNSW::HnswL2, 3, 0, 16, 10), "capacity")
  expect_error(ann$getAllNNs(matrix(numeric(), nrow = 2, ncol = 0), 1))
  expect_error(ann$getAllNNsCol(matrix(numeric(), nrow = 0, ncol = 2), 1))
})

test_that("integer coordinate storage remains accepted", {
  items <- diag(3L)
  ann <- methods::new(RcppHNSW::HnswL2, 3L, 3L, 16L, 10L)
  expect_no_error(ann$addItems(items))
  expect_equal(ann$getNNs(items[1, ], 1L), 1)
  expect_no_error(hnsw_build(items, distance = "l2"))
})
