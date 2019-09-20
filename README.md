## RcppHNSW
[![Travis-CI Build Status](https://travis-ci.org/jlmelville/rcpphnsw.svg?branch=master)](https://travis-ci.org/jlmelville/rcpphnsw)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/jlmelville/rcpphnsw?branch=master&svg=true)](https://ci.appveyor.com/project/jlmelville/rcpphnsw)
[![Coverage Status](https://img.shields.io/codecov/c/github/jlmelville/rcpphnsw/master.svg)](https://codecov.io/github/jlmelville/rcpphnsw?branch=master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
[![CRAN Monthly Downloads](https://cranlogs.r-pkg.org/badges/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)
![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/RcppHNSW)

Rcpp bindings for [hnswlib](https://github.com/nmslib/hnswlib).

### Status

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

### hnswlib

hnswlib is a header-only C++ library for finding approximate nearest neighbors
(ANN) via Hierarchical Navigable Small Worlds
[(Yashunin and Malkov, 2016)](https://arxiv.org/abs/1603.09320). 
It is part of the [nmslib](https://github.com/nmslib/nmslib]) project. 

### RcppHNSW

An R package that interfaces with hnswlib, taking enormous amounts of inspiration
from [Dirk Eddelbuettel](https://github.com/eddelbuettel)'s
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy) package which did the
same for the [Annoy](https://github.com/spotify/annoy) ANN C++ library.

One difference is that I use
[roxygen2](https://cran.r-project.org/package=roxygen2) to generate the man
pages. The `NAMESPACE` is still built manually, however (I don't believe you can
`export` the classes currently).

### Installing

From CRAN:

```R
install.packages("RcppHNSW")
```

Development versions from github:

```R
remotes::install_github("jlmelville/RcppHNSW")
```

### Function example

```R
# function interface returns results for all rows in nr x k matrices
all_knn <- RcppHNSW::hnsw_knn(data, k = 4, distance = "l2")
# other distance options: "euclidean", "cosine" and "ip" (inner product distance)

# or it can be split into two steps, so you can build with one set of data
# and search with another
ann <- hnsw_build(irism[1:100, ])
iris_nn <- hnsw_search(irism[101:150, ], ann, k = 5)
```

### Class Example

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
ann <- new(HnswL2, ncol(data), nrow(data), M, ef)

# Add items to index
for (i in 1:nrow(data)) {
  ann$addItem(data[i, ])
}

# Find 4 nearest neighbors of row 1
# indexes are in res$item, distances in res$distance
# set include_distances = TRUE to get distances as well as index
res <- ann$getNNsList(data[1, ], k = 4, include_distances = TRUE)

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

### API

#### **DO NOT USE NAMED PARAMETERS**

Because these are wrappers around C++ code, you **cannot** use named 
parameters in the calling R code. Arguments are parsed by position. This is
most annoying in constructors, which take multiple integer arguments, e.g.

```R
### DO THIS ###
num_elements <- 100
M <- 200
ef_construction <- 16
index <- new(HnswL2, dim, num_elements, M, ef)

### DON'T DO THIS ###
index <- new(HnswL2, dim, ef_construction = 16, M = 200, num_elements = 100)
# treated as if you wrote:
index <- new(HnswL2, dim, 16, 200, 100)
```

#### OK onto the API

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
* `addItem(v)` add vector `v` to the index. Internally, each vector gets an
increasing integer label, with the first vector added getting the label `1`, the
second `2` and so on. These labels are returned in `getNNs` and related methods
to identify which vector in the index are neighbors.
* `addItems(m)` add the row vectors of the matrix `m` to the index. Internally,
each row vector gets an increasing integer label, with the first row added
getting the label `1`, the second `2` and so on. These labels are returned in
`getNNs` and related methods to identify which vector in the index are
neighbors.
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
added to the index.. If `k` neighbors can't be found, an error is thrown.
* `getAllNNsList(m, k, include_distances = FALSE)` return a list containing a
matrix named `item` with the labels of the `k`-nearest neighbors of each row
vector in `m`. Labels are integers numbered from one, representing the insertion
order into the index, e.g. the label `1` represents the first item added to the
index. If `include_distances = TRUE` then also return a matrix `distance`
containing the distances. If `k` neighbors can't be found, an error is thrown.
* `size()` returns the number of items in the index. This is an upper limit on
the number of neighbors you can expect to return from `getNNs` and the other
search methods.
* `markDeleted(i)` marks the item with label `i` (the `i`th item added to the
index) as deleted. This means that the item will not be returned in any further
searches of the index. It does not reduce the memory used by the index. Calls to
`size()` do *not* reflect the number of marked deleted items.
* `resize(max_elements)` changes the maximum capacity of the index to 
`max_elements`.

### Differences from Python Bindings

* Multi-threading is not supported.
* Arbitrary integer labeling is not supported. Where labels are used, e.g. in
the return value of `getNNsList` or as input in `markDeleted`, the labels 
represent the order in which the items were added to the index, using 1-indexing
to be consistent with R. So in the Python bindings, the first item in the index
has a default of label `0`, but here it will have label `1`.
* The interface roughly follows the Python one but deviates with naming and also
rolls the declaration and initialization of the index into one call. And as
noted above, you must pass arguments by position, not keyword.
* I have made a change to the C++ `hnswalg.h` code to use the 
[`showUpdate` macro from RcppAnnoy](https://github.com/eddelbuettel/rcppannoy/blob/498a2c241df0fcac140d80f9ee0a6985d0f08687/inst/include/annoylib.h#L57),
rather than `std::cerr` directly.

### License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
