library(RcppHNSW)
context("Save/load index")

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
