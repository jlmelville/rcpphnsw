# RcppHNSW

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jlmelville/rcpphnsw?branch=master&svg=true)](https://ci.appveyor.com/project/jlmelville/rcpphnsw)
[![R-CMD-check](https://github.com/jlmelville/rcpphnsw/workflows/R-CMD-check/badge.svg)](https://github.com/jlmelville/rcpphnsw/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/jlmelville/rcpphnsw/master.svg)](https://codecov.io/github/jlmelville/rcpphnsw?branch=master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
[![CRAN Monthly Downloads](https://cranlogs.r-pkg.org/badges/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/RcppHNSW)
[![Last Commit](https://img.shields.io/github/last-commit/jlmelville/rcpphnsw)](https://github.com/jlmelville/rcpphnsw)

Rcpp bindings for [hnswlib](https://github.com/nmslib/hnswlib).

## Status

*July 18 2022* RcppHNSW 0.4.1 is released. Unfortunately, there are valgrind
problems with the version of hnswlib used in RcppHNSW 0.4.0, so that has been
rolled back.

*July 16 2022* RcppHNSW 0.4.0 is released. This release matches [hnswlib version
0.6.2](https://github.com/nmslib/hnswlib/releases/tag/v0.6.2), but otherwise
adds no new features. Some minor CRAN check NOTEs are fixed and there is also a
minor license change: previously the license was GPLv3. From this version, it
now supports GPLv3 or later.

*September 6 2020* RcppHNSW 0.3.0 is now available on CRAN, with multi-threading
support.

*August 30 2020*. Although not yet on CRAN, support for building and searching
an index in parallel (via the `n_threads` function argument and `setNumThreads`
object method) has been added to the current development version (available via
`devtools::install_github`). Thanks to 
[Dmitriy Selivanov](https://github.com/dselivanov) for a lot of the work on
this.

*September 20 2019*. RcppHNSW 0.2.0 is now available on CRAN, up to date with
hnswlib at <https://github.com/nmslib/hnswlib/commit/c5c38f0>, with new methods:
`size`, `resizeIndex` and `markDeleted`. Also, a bug that prevented searching
with datasets smaller than `k` has been fixed. Thanks to
[Yuxing Liao](https://github.com/yxngl) for spotting that.

*January 21 2019*. RcppHNSW is now available on CRAN.

*October 20 2018*. By inserting some preprocessor symbols into hnswlib, these
bindings no longer require a non-portable compiler flag and hence will pass `R
CMD CHECK` without any warnings: previously you would be warned about
`-march=native`. The price paid is not using specialized functions for the
distance calculations that are architecture-specific. I have not checked how bad
the performance hit is. The old settings remain in `src/Makevars` and
`src/Makevars.win` (commented out), if you want to build the project from
source directly. Otherwise, [Release
0.0.0.9000](https://github.com/jlmelville/rcpphnsw/releases/tag/v0.0.0.9000) is
the last version with the old behavior, which can be installed with something
like:

```R
devtools::install_github("jlmelville/rcpphnsw@v0.0.0.9000")
```

## hnswlib

hnswlib is a header-only C++ library for finding approximate nearest neighbors
(ANN) via Hierarchical Navigable Small Worlds
[(Yashunin and Malkov, 2016)](https://arxiv.org/abs/1603.09320).
It is part of the [nmslib](https://github.com/nmslib/nmslib]) project.

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
devtools::install_github("jlmelville/RcppHNSW")
```

## Function example

```R
irism <- as.matrix(iris[, -5])

# function interface returns results for all rows in nr x k matrices
all_knn <- RcppHNSW::hnsw_knn(irism, k = 4, distance = "l2")
# other distance options: "euclidean", "cosine" and "ip" (inner product distance)

# for high-dimensional data you may see a speed-up if you store the data
# where each *column* is an item to be indexed and searched. Set byrow = TRUE
# for this.
# Admittedly, the iris dataset is *not* high-dimensional
iris_by_col <- t(irism)
all_knn <- RcppHNSW::hnsw_knn(iris_by_col, k = 4, distance = "l2", byrow = FALSE)

# process can be split into two steps, so you can build with one set of data
# and search with another
ann <- hnsw_build(irism[1:100, ])
iris_nn <- hnsw_search(irism[101:150, ], ann, k = 5)
```

## Class Example

```R
library(RcppHNSW)
data <- as.matrix(iris[, -5])

# Create a new index using the L2 (squared Euclidean) distance
# nr and nc are the number of rows and columns of the data to be added, respectively
# ef and M determines speed vs accuracy trade off
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
# set include_distances = TRUE to get distances as well as index
res <- ann$getNNsList(data[1, ], k = 4, include_distances = TRUE)

# It's more efficient to use the batch methods if you have all the data you
# need at once
ann2 <- new(HnswL2, dim, nitems, M, ef)
ann2$addItems(data)
# Retrieve the 4 nearest neighbors for every item in data
res2 <- ann2$getAllNNsList(data, k = 4, include_distances = TRUE)
# labels of the data are in res$item, distances in res$distance

# If you are able to store your data column-wise, then the overhead of copying
# the data into a form usable by hnsw can be noticeably reduced
data_by_col <- t(data)
ann3 <- new(HnswL2, dim, nitems, M, ef)
ann3$addItemsCol(data_by_col)
# Retrieve the 4 nearest neighbors for every item in data_by_col
res3 <- ann3$getAllNNsListCol(data_by_col, k = 4, include_distances = TRUE)
# The returned neared neighbor data matrices are also returned column-wise
all(res2$item == t(res3$item) & res2$distance == t(res3$distance))

# other distance classes:
# Cosine: HnswCosine
# Inner Product: HnswIP
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
idx <- p$getAllNNs(data1, k = 1)
message("Recall for the first batch: ", formatC(mean(idx == 1:nrow(data1))))

# Increase the total capacity, so that it will handle the new data
p$resizeIndex(num_elements)

message("Adding the second batch of ", nrow(data2), " elements")
p$addItems(data2)

# Query the elements for themselves and measure recall:
idx <- p$getAllNNs(data, k = 1)
# You can get distances with:
# res <- p$getAllNNsList(data, k = 1, include_distances = TRUE)
# res$dist contains the distance matrix, res$item stores the indexes

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

### **DO NOT USE NAMED PARAMETERS**

Because these are wrappers around C++ code, you **cannot** use named
parameters in the calling R code. Arguments are parsed by position. This is
most annoying in constructors, which take multiple integer arguments, e.g.

```R
### DO THIS ###
dim <- 10
num_elements <- 100
M <- 200
ef_construction <- 16
index <- new(HnswL2, dim, num_elements, M, ef_construction)

### DON'T DO THIS ###
index <- new(HnswL2, dim, ef_construction = 16, M = 200, num_elements = 100)
# treated as if you wrote:
index <- new(HnswL2, dim, 16, 200, 100)
```

### OK onto the API

* `new(HnswL2, dim, max_elements, M = 16, ef_contruction = 200)` creates a new
index using the squared L2 distance (i.e. square of the Euclidean distance),
with `dim` dimensions and a maximum size of `max_elements` items. `ef` and `M`
determine the speed vs accuracy trade off. Other classes for different distances
are: `HnswCosine` for the cosine distance and `HnswIp` for the "Inner Product"
distance (like the cosine distance without normalizing).
* `new(HnswL2, dim, filename)` load a previously saved index (see `save` below)
with `dim` dimensions from the specified `filename`.
* `new(HnswL2, dim, filename, max_elements)` load a previously saved index (see
`save` below) with `dim` dimensions from the specified `filename`, and a new
maximum capacity of `max_elements`. This is a way to increase the capacity of
the index without a complete rebuild.
* `setEf(ef)` set search parameter `ef`.
* `setNumThreads(num_threads)` Use (at most) this number of threads when adding
items (via `addItems`) and searching the index (via `getAllNNs` and
`getAllNNsList`). See also the `setGrainSize` parameter.
* `setGrainSize(grain_size)` The minimum amount of work to do (adding or
searching items) per thread. If you don't have enough work for all the threads
specified by `setNumThreads` to process `grain_size` items per thread, then
fewer threads will be used. This is useful for cases where the cost of context
switching between larger number of threads would outweigh the performance gain
from parallelism. For example, if you have 100 items to process and asked for
four threads, then 25 items will be processed per thread. However, setting the 
`grain_size` to 50 will result in 50 items being processed per thread, and 
therefore only two threads being used.
* `addItem(v)` add vector `v` to the index. Internally, each vector gets an
increasing integer label, with the first vector added getting the label `1`, the
second `2` and so on. These labels are returned in `getNNs` and related methods
to identify which vector in the index are neighbors.
* `addItems(m)` add the row vectors of the matrix `m` to the index. Internally,
each row vector gets an increasing integer label, with the first row added
getting the label `1`, the second `2` and so on. These labels are returned in
`getNNs` and related methods to identify which vector in the index are
neighbors. The number of threads specified by `setNumThreads` is used for
building the index and may be non-deterministic.
* `addItemsCol(m)` Like `addItems` but adds the *column* vectors of `m` to the
index. Storing data column-wise makes copying the data for use by `hnsw` more
efficient.
* `save(filename)` saves an index to the specified `filename`. To load an index,
use the `new(HnswL2, dim, filename)` constructor (see above).
* `getNNs(v, k)` return a vector of the labels of the `k`-nearest neighbors of
the vector `v`. Labels are integers numbered from one, representing the
insertion order into the index, e.g. the label `1` represents the first item
added to the index. If `k` neighbors can't be found, an error will be thrown.
This normally means that `ef` or `M` have been set too small, but also bear in
mind that you can't return more items than were put into the index.
* `getNNsList(v, k, include_distances = FALSE)` return a list containing a
vector named `item` with the labels of the `k`-nearest neighbors of the vector
`v`. Labels are integers numbered from one, representing the insertion order
into the index, e.g. the label `1` represents the first item added to the index.
If `include_distances = TRUE` then also return a vector `distance` containing
the distances. If `k` neighbors can't be found, an error is thrown.
* `getAllNNs(m, k)` return a matrix of the labels of the `k`-nearest neighbors
of each row vector in `m`. Labels are integers numbered from one, representing
the insertion order into the index, e.g. the label `1` represents the first item
added to the index. If `k` neighbors can't be found, an error is thrown. The 
number  of threads specified by `setNumThreads` is used for searching.
* `getAllNNsList(m, k, include_distances = FALSE)` return a list containing a
matrix named `item` with the labels of the `k`-nearest neighbors of each row
vector in `m`. Labels are integers numbered from one, representing the insertion
order into the index, e.g. the label `1` represents the first item added to the
index. If `include_distances = TRUE` then also return a matrix `distance`
containing the distances. If `k` neighbors can't be found, an error is thrown.
The number  of threads specified by `setNumThreads` is used for searching.
* `getAllNNsCol(m, k)` like `getAllNNs` but each item to be searched in `m` is
stored by *column*, not row. In addition the returned matrix of `k`-nearest
neighbors is also stored column-wise: i.e. the dimension of the return value
matrix is `k x n` where `n` is the number of items (columns) in `m`. By passing
the data column-wise, some overhead associated with copying data to and from
`hnsw` can be reduced.
* `getAllNNsListCol(m, k)` like `getAllNNsList` but each item to be searched in
`m` is stored by *column*, not row. In addition, the matrices in the returned 
list are also stored column-wise: i.e. the dimension of the return value matrix
is `k x n` where `n` is the number of items (columns) in `m`. By passing the 
data column-wise, some overhead associated with copying data to and from `hnsw`
can be reduced.
* `size()` returns the number of items in the index. This is an upper limit on
the number of neighbors you can expect to return from `getNNs` and the other
search methods.
* `markDeleted(i)` marks the item with label `i` (the `i`th item added to the
index) as deleted. This means that the item will not be returned in any further
searches of the index. It does not reduce the memory used by the index. Calls to
`size()` do *not* reflect the number of marked deleted items.
* `resize(max_elements)` changes the maximum capacity of the index to
`max_elements`.

## Differences from Python Bindings

* Arbitrary integer labeling is not supported. Where labels are used, e.g. in
the return value of `getNNsList` or as input in `markDeleted`, the labels
represent the order in which the items were added to the index, using 1-indexing
to be consistent with R. So in the Python bindings, the first item in the index
has a default of label `0`, but here it will have label `1`.
* The interface roughly follows the Python one but deviates with naming and also
rolls the declaration and initialization of the index into one call. And as
noted above, you must pass arguments by position, not keyword.

## License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
