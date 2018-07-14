## RcppHNSW
[![Travis-CI Build Status](https://travis-ci.org/jlmelville/rcpphnsw.svg?branch=master)](https://travis-ci.org/jlmelville/rcpphnsw)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jlmelville/rcpphnsw?branch=master&svg=true)](https://ci.appveyor.com/project/jlmelville/rcpphnsw)

Rcpp bindings for [HNSW](https://github.com/nmslib/hnsw).

### Status

Don't get too excited. Barely any of HNSW has been exposed so far. But it can
calculate nearest neighbors via Euclidean distance.

### HNSW

HNSW is a header-only C++ library for finding approximate nearest neighbors
(ANN) via Hierarchical Navigable Small Worlds
[(Yashunin and Malkov, 2016)](https://arxiv.org/abs/1603.09320). 
It is part of the [nmslib](https://github.com/nmslib/nmslib]) project. 

### RcppHNSW

An R package that interfaces with HNSW, taking enormous amounts of inspiration from 
[Dirk Eddelbuettel](https://github.com/eddelbuettel)'s 
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy) package which did the same for
the [Annoy](https://github.com/spotify/annoy) ANN C++ library. 

One difference is that I use
[roxygen2](https://cran.r-project.org/package=roxygen2) to generate the man
pages. The `NAMESPACE` is still built manually, however (I don't believe you can
`export` the classes currently).

### Installing

```R
devtools::install_github("jlmelville/RcppHNSW")
```

### Example

```R
data <- as.matrix(iris[, -5])

# Create a new index using the L2 (squared Euclidean) distance
# nr and nc are the number of rows and columns of the data to be added, respectively
# ef and M determines speed vs accuracy trade off
ann <- methods::new(RcppHNSW::HnswL2, nc = ncol(data), nr = nrow(data), 
                    ef = 200, M = 16)

# Add items to index
for (i in 1:nr) {
  ann$addItem(data[i, ])
}

# Find 4 nearest neighbors of row 1
# indexes are in res$item, distances in res$distance
res <- ann$getNNsList(data[1, ], k = 4, include_distances = TRUE)

# function interface returns results for all rows in nr x k matrices
all_knn <- RcppHNSW::get_knn(data, k = 4)
```

Here's a rough equivalent of the serialization/deserialization example from
the [hnsw README](https://github.com/nmslib/hnsw#python-bindings-example), but 
without any multithreading:

```R
library("RcppHNSW")

dim <- 16
num_elements <- 10000

# Generate sample data
data <- matrix(stats::runif(num_elements * dim), nrow = num_elements)

# Create index
p <- new(HnswL2, dim, num_elements, ef = 10, M = 16)

# Split data into two batches
data1 <- data[1:(num_elements / 2), ]
data2 <- data[(num_elements / 2 + 1):num_elements, ]

message("Adding first batch of ", nrow(data1), " elements")
for (i in 1:nrow(data1)) {
  p$addItem(data1[i, ])
}

# Query the elements for themselves and measure recall:
idx <- rep(0, nrow(data1))
for (i in 1:nrow(data1)) {
  idx[i] <- p$getNNs(data1[i, ], k = 1)
}

message("Recall for the first batch: ", formatC(mean(idx == 0:(nrow(data1) - 1))))

filename <- "first_half.bin"
# Serialize index
p$save(filename)

# Reinitialize and load the index
rm(p)
message("Loading index from ", filename)
p <- new(HnswL2, dim, filename)

message("Adding the second batch of ", nrow(data2), " elements")
for (i in 1:nrow(data2)) {
  p$addItem(data2[i, ])
}

# Query the elements for themselves and measure recall:
idx <- rep(0, num_elements)
for (i in 1:num_elements) {
  # Neighbors are queried by passing the vector back in
  idx[i] <- p$getNNs(data[i, ], k = 1)
}
message("Recall for two batches: ", formatC(mean(idx == 0:(num_elements - 1))))
```

### Differences from Python Bindings

* No multi-threading support.
* You can only add and retrieve items one at a time, i.e. to add or retrieve
results for a matrix of data, you will need to manually loop over it.
* Arbitrary integer labelling is not supported. Items are labelled 
`0, 1, 2 ... N`.
* The interface roughly follows the Python one but deviates with naming and also
rolls the declaration and initialization of the index into one call.

### Note

I had to add a non-portable flag to `PKG_CPPFLAGS` (`-march=native`). This should be
the only status warning in `R CMD check`.

### License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
