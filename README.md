# RcppHNSW

[![R-CMD-check](https://github.com/jlmelville/rcpphnsw/workflows/R-CMD-check/badge.svg)](https://github.com/jlmelville/rcpphnsw/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/jlmelville/rcpphnsw/master.svg)](https://codecov.io/github/jlmelville/rcpphnsw?branch=master)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
[![Dependencies](https://tinyverse.netlify.app/badge/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
[![CRAN Monthly Downloads](https://cranlogs.r-pkg.org/badges/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/RcppHNSW)
[![Last Commit](https://img.shields.io/github/last-commit/jlmelville/rcpphnsw)](https://github.com/jlmelville/rcpphnsw)

Rcpp bindings for [hnswlib](https://github.com/nmslib/hnswlib).

## Status

*May 26 2026* RcppHNSW 0.7.0 is released to CRAN to fix an undefined behavior
issue in the underlying `hnswlib` library. Also, the `random_seed` parameter
is now exposed in the API.

## hnswlib

hnswlib is a header-only C++ library for finding approximate nearest neighbors
(ANN) via Hierarchical Navigable Small Worlds
[(Malkov and Yashunin, 2020)](https://doi.org/10.1109/TPAMI.2018.2889473).
It is part of the [nmslib](https://github.com/nmslib/nmslib) project.

## The RcppHNSW Package

An R package that interfaces with hnswlib, taking enormous amounts of inspiration
from [Dirk Eddelbuettel](https://github.com/eddelbuettel)'s
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy) package which did the
same for the [Annoy](https://github.com/spotify/annoy) ANN C++ library.

One difference is that I use
[roxygen2](https://cran.r-project.org/package=roxygen2) to generate the man
pages. The `NAMESPACE` is still built manually, however (I don't believe you can
`export` the classes currently).

## Installing

From CRAN:

```R
install.packages("RcppHNSW")
```

Development versions from github:

```R
pak::pak("jlmelville/RcppHNSW")
```

## Function example

```R
irism <- as.matrix(iris[, -5])

# function interface returns results for all rows in nr x k matrices
all_knn <- RcppHNSW::hnsw_knn(irism, k = 4, distance = "l2")
# other distance options: "euclidean", "cosine" and "ip" (inner product distance)

# for high-dimensional data you may see a speed-up if you store the data
# where each *column* is an item to be indexed and searched. Set byrow = FALSE
# for this.
# Admittedly, the iris dataset is *not* high-dimensional
iris_by_col <- t(irism)
all_knn <- RcppHNSW::hnsw_knn(iris_by_col, k = 4, distance = "l2", byrow = FALSE)

# process can be split into two steps, so you can build with one set of data
# and search with another
ann <- RcppHNSW::hnsw_build(irism[1:100, ])
iris_nn <- RcppHNSW::hnsw_search(irism[101:150, ], ann, k = 5)
```

## Class Example

As noted in the "Module arguments are positional" section below, use local
variables or comments to document the meaning of class-method arguments; do
not rely on argument names for matching.

```R
library(RcppHNSW)
data <- as.matrix(iris[, -5])

# Create a new index using the L2 (squared Euclidean) distance
# nitems and dim are the number of rows and columns of the data, respectively
# M and ef determine the speed-versus-accuracy trade-off
# You must specify the maximum number of items to add to the index when it
# is created. But you can increase this number: see the next example
M <- 16
ef <- 200
dim <- ncol(data)
nitems <- nrow(data)
ann <- new(HnswL2, dim, nitems, M, ef)

# Add items to index
for (i in 1:nitems) {
  ann$addItem(data[i, ])
}

# Find 4 nearest neighbors of row 1
# indexes are in res$item, distances in res$distance
# The positional arguments are query, k, and whether to include distances
res <- ann$getNNsList(data[1, ], 4, TRUE)

# It's more efficient to use the batch methods if you have all the data you
# need at once
ann2 <- new(HnswL2, dim, nitems, M, ef)
ann2$addItems(data)
# Retrieve the 4 nearest neighbors for every item in data
res2 <- ann2$getAllNNsList(data, 4, TRUE)
# labels of the data are in res2$item, distances in res2$distance

# If you are able to store your data column-wise, then the overhead of copying
# the data into a form usable by hnsw can be noticeably reduced
data_by_col <- t(data)
ann3 <- new(HnswL2, dim, nitems, M, ef)
ann3$addItemsCol(data_by_col)
# Retrieve the 4 nearest neighbors for every item in data_by_col
res3 <- ann3$getAllNNsListCol(data_by_col, 4, TRUE)
# The returned nearest-neighbor data matrices are also returned column-wise
all(res2$item == t(res3$item) & res2$distance == t(res3$distance))

# Save the index
ann$save("iris.hnsw")

# load it back in: you do need to know the dimension of the original data
ann4 <- new(HnswL2, dim, "iris.hnsw")
# new index should behave like the original
all(ann$getNNs(data[1, ], 4) == ann4$getNNs(data[1, ], 4))
unlink("iris.hnsw")

# other distance classes:
# Cosine: HnswCosine
# Inner Product: HnswIp
# Euclidean: HnswEuclidean
```

Here's a rough equivalent of the serialization/deserialization example from
the
[hnswlib README](https://github.com/nmslib/hnswlib#python-bindings-examples),
but using the recently-added `resizeIndex` method to increase the size of the
index after its initial specification, avoiding having to read from or write
to disk:

```R
library("RcppHNSW")
set.seed(12345)

dim <- 16
num_elements <- 100000

# Generate sample data
data <- matrix(stats::runif(num_elements * dim), nrow = num_elements)

# Split data into two batches
data1 <- data[1:(num_elements / 2), ]
data2 <- data[(num_elements / 2 + 1):num_elements, ]

# Create index
M <- 16
ef <- 10
# Set the initial index size to the size of the first batch
p <- new(HnswL2, dim, num_elements / 2, M, ef)

message("Adding first batch of ", nrow(data1), " elements")
p$addItems(data1)

# Query the elements for themselves and measure recall:
idx <- p$getAllNNs(data1, 1)
message("Recall for the first batch: ", formatC(mean(idx == 1:nrow(data1))))

# Increase the total capacity, so that it will handle the new data
p$resizeIndex(num_elements)

message("Adding the second batch of ", nrow(data2), " elements")
p$addItems(data2)

# Query the elements for themselves and measure recall:
idx <- p$getAllNNs(data, 1)
# You can get distances with:
# res <- p$getAllNNsList(data, 1, TRUE)
# res$distance contains the distance matrix, res$item stores the indexes

message("Recall for two batches: ", formatC(mean(idx == 1:num_elements)))
```

Although there's no longer any need for this, for completeness, here's how you
would use `save` and `new` to achieve the same effect without `resizeIndex`:

```R
filename <- "first_half.bin"
# Serialize index
p$save(filename)

# Reinitialize and load the index
rm(p)
message("Loading index from ", filename)
# Increase the total capacity, so that it will handle the new data
p <- new(HnswL2, dim, filename, num_elements)
unlink(filename)
```

## API

### **MODULE ARGUMENTS ARE POSITIONAL**

Rcpp Module arguments are matched by position. Supplying names does not make
reordering safe. Use explanatory variables and pass them in the documented
order, especially for constructors with several integer arguments:

```R
dim <- 10
num_elements <- 100
M <- 16
ef_construction <- 200
index <- new(HnswL2, dim, num_elements, M, ef_construction)
```

### OK onto the API

* `new(HnswL2, dim, max_elements, M, ef_construction)` creates a new
index using the squared L2 distance (i.e. square of the Euclidean distance),
with `dim` dimensions and a maximum size of `max_elements` items. `M` and
`ef_construction` determine the speed versus accuracy trade-off. Other classes
for different distances are `HnswCosine` for cosine, `HnswIp` for inner
product (`1 - sum(a * b)`), and `HnswEuclidean` for Euclidean distance.
* `new(HnswL2, dim, max_elements, M, ef_construction, random_seed)` is the same
as the previous constructor, but with a specified random seed. Omitting
`random_seed` uses hnswlib's default of `100`.
* `new(HnswL2, dim, filename)` load a previously saved index (see `save` below)
with `dim` dimensions from the specified `filename`.
* `new(HnswL2, dim, filename, max_elements)` load a previously saved index (see
`save` below) with `dim` dimensions from the specified `filename`, and a new
maximum capacity of `max_elements`. This is a way to increase the capacity of
the index without a complete rebuild.
* `setEf(ef)` set search parameter `ef`.
* `setNumThreads(num_threads)` Use (at most) this number of threads when adding
items (via `addItems`) and searching the index (via `getAllNNs` and
`getAllNNsList`). Zero and one both select serial execution. See also the
`setGrainSize` method.
* `setGrainSize(grain_size)` The minimum amount of work to do (adding or
searching items) per thread. Zero is treated as one. If you don't have enough
work for all the threads specified by `setNumThreads` to process `grain_size`
items per thread, then fewer threads will be used. This is useful for cases
where the cost of context switching between a larger number of threads would
outweigh the performance gain from parallelism. For example, if you have 100
items to process and ask for four threads, then 25 items will be processed per
thread. Setting `grain_size` to 50 instead will use two threads with 50 items
each.
* `addItem(v)` add vector `v` to the index. Internally, each vector gets an
increasing integer label, with the first vector added getting the label `1`, the
second `2` and so on. These labels are returned in `getNNs` and related methods
to identify which vector in the index are neighbors.
* `addItems(m)` add the row vectors of the matrix `m` to the index. Internally,
each row vector gets an increasing integer label, with the first row added
getting the label `1`, the second `2` and so on. These labels are returned in
`getNNs` and related methods to identify which vector in the index are
neighbors. The number of threads specified by `setNumThreads` is used for
building the index. Parallel construction may be nondeterministic even with a
fixed `random_seed`; R's `set.seed()` does not control hnswlib.
* `addItemsCol(m)` Like `addItems` but adds the *column* vectors of `m` to the
index. Storing data column-wise makes copying the data for use by `hnsw` more
efficient.
* `save(filename)` saves a raw index checkpoint to the specified `filename`. To
load an index, use the `new(HnswL2, dim, filename)` constructor (see above).
* `getItems(ids)` returns a matrix where each row is the data vector from the
index associated with integer indices in the vector of `ids`. For cosine
similarity, the l2 row-normalized vectors are returned. `ids` are one-indexed,
i.e. to get the first and tenth vectors that were added to the index, use
`getItems(c(1, 10))`, not `getItems(c(0, 9))`.
* `getNNs(v, k)` returns a vector of the labels of the `k`-nearest neighbors of
the vector `v`. Labels are integers numbered from one, representing the
insertion order into the index, e.g. the label `1` represents the first item
added to the index. `k` must be positive and cannot exceed the active (not
deleted) item count. If `k` neighbors can't be found, an error will be thrown.
This can also mean that `ef` or `M` have been set too small.
* `getNNsList(v, k, include_distances)` returns a list containing a
vector named `item` with the labels of the `k`-nearest neighbors of the vector
`v`. Labels are integers numbered from one, representing the insertion order
into the index, e.g. the label `1` represents the first item added to the index.
If `include_distances = TRUE` then also return a vector `distance` containing
the distances. If `k` neighbors can't be found, an error is thrown.
* `getAllNNs(m, k)` returns a matrix of the labels of the `k`-nearest neighbors
of each row vector in `m`. Labels are integers numbered from one, representing
the insertion order into the index, e.g. the label `1` represents the first item
added to the index. If `k` neighbors can't be found, an error is thrown. The
number of threads specified by `setNumThreads` is used for searching.
* `getAllNNsList(m, k, include_distances)` returns a list containing a
matrix named `item` with the labels of the `k`-nearest neighbors of each row
vector in `m`. Labels are integers numbered from one, representing the insertion
order into the index, e.g. the label `1` represents the first item added to the
index. If `include_distances = TRUE` then also return a matrix `distance`
containing the distances. If `k` neighbors can't be found, an error is thrown.
The number of threads specified by `setNumThreads` is used for searching.
* `getAllNNsCol(m, k)` like `getAllNNs` but each item to be searched in `m` is
stored by *column*, not row. In addition the returned matrix of `k`-nearest
neighbors is also stored column-wise: i.e. the dimension of the return value
matrix is `k x n` where `n` is the number of items (columns) in `m`. By passing
the data column-wise, some overhead associated with copying data to and from
`hnsw` can be reduced.
* `getAllNNsListCol(m, k, include_distances)` is like `getAllNNsList` but each
item to be searched in `m` is stored by *column*, not row. The matrices in the
returned list are also stored column-wise: their dimensions are `k x n`, where
`n` is the number of items (columns) in `m`. By passing the data column-wise,
some overhead associated with copying data to and from `hnsw` can be reduced.
* `size()` returns the total number of items added to the index, including
deleted items. It can therefore be larger than the maximum valid `k` after
deletion.
* `markDeleted(i)` marks the item with label `i` (the `i`th item added to the
index) as deleted. This means that the item will not be returned in any further
searches of the index and `getItems()` will no longer return it. Deletion does
not reduce memory use or reclaim capacity, and `size()` still includes it.
* `resizeIndex(max_elements)` changes the maximum capacity of the index to
`max_elements`.

All input coordinates are converted to and stored as single-precision floats.
They must be finite and representable in that format. Cosine items and queries
must have a positive finite norm after conversion. If an exception escapes
after insertion has begun, the index becomes unusable and must be discarded
and rebuilt or reloaded.

Raw files created by `save()` are operational hnswlib checkpoints, not a stable
cross-version or cross-platform archive. Load with the exact original
dimension and normally the same class. `HnswL2` and `HnswEuclidean` can load
each other's same-width checkpoints; no equivalent semantic guarantee is made
for cosine or inner-product reinterpretation. Checkpoints assume RcppHNSW's
contiguous insertion-order labels. Loading restores deletion state but resets
search `ef` to 10, so call `setEf()` before direct Module search if another
value is required.

## Differences from Python Bindings

* Arbitrary integer labeling is not supported. Where labels are used, e.g. in
the return value of `getNNsList` or as input in `markDeleted` or `getItems`, the
labels represent the order in which the items were added to the index, using
1-indexing to be consistent with R. So in the Python bindings, the first item in
the index has a default of label `0`, but here it will have label `1`.
* The interface roughly follows the Python one but deviates with naming and also
rolls the declaration and initialization of the index into one call. And as
noted above, you must pass arguments by position, not keyword.

## License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
