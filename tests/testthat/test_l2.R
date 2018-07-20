library(RcppHNSW)
context("L2 distances")

num_elements <- nrow(uirism)
dim <- ncol(uirism)

p <- new(HnswL2, dim, num_elements, ef = 10, M = 16)

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
p <- new(HnswL2, dim, num_elements, ef = 10, M = 16)

for (i in 1:(floor(num_elements / 2))) {
  p$addItem(uirism[i, ])
}

filename <- "first_half.bin"
p$save(filename)
rm(p)
p <- new(HnswL2, dim, filename)
unlink(filename)

for (i in (floor(num_elements / 2) + 1):num_elements) {
  p$addItem(uirism[i, ])
}

idx <- rep(0, num_elements)
for (i in 1:num_elements) {
  idx[i] <- p$getNNs(uirism[i, ], k = 1)
}
serde_recall <- mean(idx == 1:num_elements)
expect_equal(serde_recall, recall)
