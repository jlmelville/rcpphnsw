library(RcppHNSW)
context("Resize index")

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
p <- new(HnswL2, dim, num_elements / 2, 16, 10)

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
