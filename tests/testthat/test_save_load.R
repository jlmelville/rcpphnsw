library(RcppHNSW)
test_that("save and load an L2 index", {
  num_elements <- nrow(uirism)
  dim <- ncol(uirism)

  M <- 16
  ef_construction <- 10
  p <- new(HnswL2, dim, num_elements, M, ef_construction)

  for (i in 1:num_elements) {
    p$addItem(uirism[i, ])
  }

  nn4idx <- matrix(0L, nrow = num_elements, ncol = 4)
  nn4dist <- matrix(0.0, nrow = num_elements, ncol = 4)

  for (i in 1:num_elements) {
    res <- p$getNNsList(uirism[i, ], k = 4, TRUE)
    nn4idx[i, ] <- res$item
    nn4dist[i, ] <- res$distance
  }

  temp_file <- tempfile()
  on.exit(unlink(temp_file), add = TRUE)
  p$save(temp_file)

  nn4idx_aftersave <- matrix(0L, nrow = num_elements, ncol = 4)
  nn4dist_aftersave <- matrix(0.0, nrow = num_elements, ncol = 4)
  for (i in 1:num_elements) {
    res_aftersave <- p$getNNsList(uirism[i, ], k = 4, TRUE)
    nn4idx_aftersave[i, ] <- res_aftersave$item
    nn4dist_aftersave[i, ] <- res_aftersave$distance
  }
  expect_equal(nn4idx, nn4idx_aftersave)
  expect_equal(nn4dist, nn4dist_aftersave)

  pload <- new(HnswL2, dim, temp_file)
  nn4idx_afterload <- matrix(0L, nrow = num_elements, ncol = 4)
  nn4dist_afterload <- matrix(0.0, nrow = num_elements, ncol = 4)
  for (i in 1:num_elements) {
    res_afterload <- pload$getNNsList(uirism[i, ], k = 4, TRUE)
    nn4idx_afterload[i, ] <- res_afterload$item
    nn4dist_afterload[i, ] <- res_afterload$distance
  }
  expect_equal(nn4idx, nn4idx_afterload)
  expect_equal(nn4dist, nn4dist_afterload)
})

# 21: no way to use hnsw_search to get Euclidean distances after save/load
test_that("euclidean search is more consistent with save/load", {
  ann <- hnsw_build(ui10, distance = "euclidean")
  iris_nn <- hnsw_search(ui10, ann, k = 4)
  expect_equal(iris_nn$dist, self_nn_dist4, tolerance = 1e-6)
  expect_equal(iris_nn$idx, self_nn_index4)

  temp_file <- tempfile()
  on.exit(unlink(temp_file), add = TRUE)
  ann$save(temp_file)

  ann2 <- methods::new(RcppHNSW::HnswEuclidean, 4, temp_file)
  ann$markDeleted(1)
  mutated_nn <- hnsw_search(ui10, ann, k = 4)
  expect_false(any(mutated_nn$idx == 1L))

  iris_nn2 <- hnsw_search(ui10, ann2, k = 4)
  expect_equal(iris_nn2$dist, self_nn_dist4, tolerance = 1e-6)
  expect_equal(iris_nn2$idx, self_nn_index4)
})

test_that("every Module class searches a loaded raw index", {
  items <- rbind(
    c(1, 0, 0),
    c(0, 2, 0),
    c(1, 1, 1),
    c(-1, 1, 2)
  )
  classes <- list(
    l2 = RcppHNSW::HnswL2,
    euclidean = RcppHNSW::HnswEuclidean,
    cosine = RcppHNSW::HnswCosine,
    ip = RcppHNSW::HnswIp
  )

  for (name in names(classes)) {
    path <- tempfile(fileext = ".hnsw")
    on.exit(unlink(path), add = TRUE)
    ann <- methods::new(classes[[name]], ncol(items), nrow(items), 16, 50)
    ann$addItems(items)
    before <- ann$getAllNNsList(items, 2, TRUE)
    ann$save(path)

    loaded <- methods::new(
      classes[[name]],
      ncol(items),
      path,
      nrow(items) + 1
    )
    after <- loaded$getAllNNsList(items, 2, TRUE)

    expect_identical(loaded$size(), as.double(nrow(items)), info = name)
    expect_identical(after$item, before$item, info = name)
    expect_equal(after$distance, before$distance, tolerance = 1e-6, info = name)
    for (wrong_dimension in c(ncol(items) - 1, ncol(items) + 1)) {
      expect_error(
        methods::new(classes[[name]], wrong_dimension, path),
        "data size is incompatible with the requested space",
        info = paste(name, wrong_dimension)
      )
    }

    if (name == "l2") {
      added_item <- c(2, 2, 2)
      loaded$addItem(added_item)
      expect_identical(loaded$size(), as.double(nrow(items) + 1L))
      expect_equal(
        loaded$getItems(nrow(items) + 1L),
        matrix(added_item, nrow = 1)
      )
      expect_identical(loaded$getNNs(added_item, 1), nrow(items) + 1)
    }
  }
})

test_that("L2 checkpoints retain Euclidean reinterpretation", {
  items <- rbind(
    c(1, 0, 0),
    c(0, 2, 0),
    c(1, 1, 1),
    c(-1, 1, 2)
  )
  path <- tempfile(fileext = ".hnsw")
  on.exit(unlink(path), add = TRUE)
  l2 <- methods::new(RcppHNSW::HnswL2, ncol(items), nrow(items), 16, 50)
  l2$addItems(items)
  l2$save(path)

  euclidean <- methods::new(RcppHNSW::HnswEuclidean, ncol(items), path)
  l2_result <- l2$getAllNNsList(items, 2, TRUE)
  euclidean_result <- euclidean$getAllNNsList(items, 2, TRUE)

  expect_identical(euclidean_result$item, l2_result$item)
  expect_equal(
    euclidean_result$distance,
    sqrt(l2_result$distance),
    tolerance = 1e-6
  )
})

test_that("raw loading rejects missing, truncated, and corrupt files", {
  expect_error(
    methods::new(RcppHNSW::HnswL2, 3, tempfile()),
    "Cannot open file"
  )
  expect_error(
    methods::new(RcppHNSW::HnswL2, 3, ""),
    "^path_to_index "
  )

  valid_path <- tempfile(fileext = ".hnsw")
  truncated_path <- tempfile(fileext = ".hnsw")
  invalid_offsets_path <- tempfile(fileext = ".hnsw")
  invalid_width_path <- tempfile(fileext = ".hnsw")
  paths <- c(
    valid_path,
    truncated_path,
    invalid_offsets_path,
    invalid_width_path
  )
  on.exit(unlink(paths), add = TRUE)

  ann <- methods::new(RcppHNSW::HnswL2, 3, 3, 16, 50)
  ann$addItems(diag(3))
  ann$save(valid_path)
  serialized <- readBin(
    valid_path,
    what = "raw",
    n = file.info(valid_path)$size
  )

  writeBin(serialized[seq_len(20)], truncated_path)
  expect_error(
    methods::new(RcppHNSW::HnswL2, 3, truncated_path),
    "invalid serialized offsets|data size is incompatible|corrupted"
  )

  size_t_bytes <- .Machine$sizeof.pointer
  label_offset <- 4L * size_t_bytes + seq_len(size_t_bytes)
  data_offset <- 5L * size_t_bytes + seq_len(size_t_bytes)

  invalid_offsets <- serialized
  invalid_offsets[label_offset] <- as.raw(0)
  writeBin(invalid_offsets, invalid_offsets_path)
  expect_error(
    methods::new(RcppHNSW::HnswL2, 3, invalid_offsets_path),
    "invalid serialized offsets"
  )

  invalid_width <- serialized
  invalid_width[label_offset] <- invalid_width[data_offset]
  writeBin(invalid_width, invalid_width_path)
  expect_error(
    methods::new(RcppHNSW::HnswL2, 3, invalid_width_path),
    "data size is incompatible with the requested space"
  )
})

test_that("raw save failures reach R without poisoning the index", {
  items <- diag(3)
  ann <- methods::new(RcppHNSW::HnswL2, 3, nrow(items), 16, 50)
  ann$addItems(items)

  expect_error(ann$save(""), "^path_to_index ")
  missing_parent <- file.path(tempfile(), "index.hnsw")
  expect_error(ann$save(missing_parent), "Cannot open file for writing")
  expect_false(file.exists(missing_parent))

  directory_path <- tempfile()
  dir.create(directory_path)
  on.exit(unlink(directory_path, recursive = TRUE), add = TRUE)
  expect_error(ann$save(directory_path), "Cannot open file for writing")

  if (file.exists("/dev/full")) {
    expect_error(ann$save("/dev/full"), "Failed to write index to file")
  }

  expect_identical(ann$size(), as.double(nrow(items)))
  expect_identical(ann$getNNs(items[1, ], 1), 1)
})
