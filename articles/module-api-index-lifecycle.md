# Module API and index lifecycle

## Classes

| Class           | Distance                                   |
|-----------------|--------------------------------------------|
| `HnswL2`        | Squared L2, i.e. squared Euclidean.        |
| `HnswEuclidean` | Euclidean.                                 |
| `HnswCosine`    | One minus cosine similarity.               |
| `HnswIp`        | One minus inner product: `1 - sum(a * b)`. |

All four classes have the same constructors and methods.

## Constructors

Rcpp Module arguments are matched by position. Supplying names does not
make reordering safe. Use explanatory variables and pass them in the
documented order, especially for constructors with several integer
arguments:

``` r

dim <- 10
num_elements <- 100
M <- 16
ef_construction <- 200
index <- new(HnswL2, dim, num_elements, M, ef_construction)
```

Replace `HnswL2` with any of the other classes in these signatures:

``` r

new(HnswL2, dim, max_elements, M, ef_construction)
new(HnswL2, dim, max_elements, M, ef_construction, random_seed)
new(HnswL2, dim, filename)
new(HnswL2, dim, filename, max_elements)
```

The first constructor creates a new index with `dim` dimensions and a
maximum size of `max_elements` items. `M` and `ef_construction`
determine the speed versus accuracy trade-off. The second form specifies
the random seed; omitting `random_seed` uses hnswlib’s default of `100`.

The filename constructors load a previously saved index with `dim`
dimensions. The final form supplies a new maximum capacity. This is a
way to increase the capacity of the index without a complete rebuild.

## Method lookup

| Category | Methods |
|----|----|
| Add and retrieve items | `addItem()`, `addItems()`, `addItemsCol()`, `getItems()` |
| Search | `getNNs()`, `getNNsList()`, `getAllNNs()`, `getAllNNsList()`, `getAllNNsCol()`, `getAllNNsListCol()` |
| Search and parallel settings | `setEf()`, `setNumThreads()`, `setGrainSize()` |
| Index lifecycle | `size()`, `markDeleted()`, `resizeIndex()`, `ann$save()` |

### Add and retrieve items

- `addItem(v)` adds vector `v` to the index.
- `addItems(m)` adds the row vectors of the matrix `m` to the index. The
  number of threads specified by `setNumThreads()` is used for building
  the index.
- `addItemsCol(m)` is like `addItems()` but adds the *column* vectors of
  `m` to the index. Storing data column-wise makes copying the data for
  use by hnsw more efficient.
- `getItems(ids)` returns a matrix where each row is the data vector
  from the index associated with integer indices in the vector of `ids`.
  For cosine similarity, the L2 row-normalized vectors are returned.

Internally, each item gets an increasing integer label, with the first
item added getting the label `1`, the second `2` and so on. These labels
are returned in search results to identify which vectors in the index
are neighbors. `getItems()` identifiers are also one-indexed: to get the
first and tenth vectors added to the index, use `getItems(c(1, 10))`,
not `getItems(c(0, 9))`.

### Search

- `getNNs(v, k)` returns a vector of the labels of the `k`-nearest
  neighbors of vector `v`.
- `getNNsList(v, k, include_distances)` returns a list containing a
  vector named `item` with the labels of the `k`-nearest neighbors of
  vector `v`. If `include_distances = TRUE`, it also returns a vector
  named `distance`.
- `getAllNNs(m, k)` returns a matrix of the labels of the `k`-nearest
  neighbors of each row vector in `m`.
- `getAllNNsList(m, k, include_distances)` returns a list containing a
  matrix named `item` with the labels of the `k`-nearest neighbors of
  each row vector in `m`. If `include_distances = TRUE`, it also returns
  a matrix named `distance`.
- `getAllNNsCol(m, k)` is like `getAllNNs()` but each item to be
  searched in `m` is stored by *column*, not row. The returned matrix is
  also stored column-wise: its dimensions are `k x n`, where `n` is the
  number of columns in `m`.
- `getAllNNsListCol(m, k, include_distances)` is like `getAllNNsList()`
  but each item to be searched in `m` is stored by *column*. The
  matrices in the returned list are also stored column-wise, with
  dimensions `k x n`.

The number of threads specified by `setNumThreads()` is used by the
batch search methods. `k` must be positive and cannot exceed the active
(not deleted) item count. If `k` neighbors cannot be found, an error is
thrown. This can also mean that `ef` or `M` have been set too small.

### Search and parallel settings

- `setEf(ef)` sets search parameter `ef`.
- `setNumThreads(num_threads)` uses at most this number of threads when
  adding items via `addItems()` or `addItemsCol()` and searching via the
  batch search methods. Zero and one both select serial execution.
- `setGrainSize(grain_size)` sets the minimum amount of work to do per
  thread. Zero is treated as one. If there is not enough work for all
  the threads to process `grain_size` items per thread, fewer threads
  will be used.

Parallel construction may be nondeterministic even with a fixed
`random_seed`. R’s [`set.seed()`](https://rdrr.io/r/base/Random.html)
does not control hnswlib.

### Index lifecycle

- `size()` returns the total number of items added to the index,
  including deleted items. It can therefore be larger than the maximum
  valid `k` after deletion.
- `markDeleted(i)` marks the item with label `i` as deleted. The item
  will not be returned in further searches and `getItems()` will no
  longer return it. Deletion does not reduce memory use or reclaim
  capacity, and `size()` still includes it.
- `resizeIndex(max_elements)` changes the maximum capacity of the index.
- `ann$save(filename)` saves a raw index checkpoint to `filename`.

Coordinates are stored as single-precision floats. The package rejects
non-finite or out-of-range coordinates and cosine vectors with zero norm
after conversion. If an exception escapes after insertion has begun, the
index becomes unusable and must be discarded and rebuilt or reloaded.

## Raw index checkpoints

`ann$save()` uses hnswlib’s raw checkpoint format. Compatibility depends
on the hnswlib version and platform. Load with the exact original
dimension and normally the same class. Same-width `HnswL2` and
`HnswEuclidean` checkpoints support cross-loading. Use matching classes
for cosine and inner-product checkpoints.

Checkpoints assume RcppHNSW’s contiguous insertion-order labels. Loading
restores deletion state but resets search `ef` to 10, so call `setEf()`
before direct Module search if another value is required. An optional
load capacity may enlarge the index, but cannot be smaller than the
stored item count.

## Example

``` r

library(RcppHNSW)

data <- as.matrix(iris[, -5])
item_dim <- ncol(data)
M <- 16
ef_construction <- 200

ann <- new(HnswL2, item_dim, 100, M, ef_construction, 100)
ann$setNumThreads(2)
ann$setGrainSize(10)
ann$addItems(data[1:100, ])

ann$resizeIndex(nrow(data))
ann$addItems(data[101:150, ])
ann$size()
```

    ## [1] 150

``` r

result <- ann$getAllNNsList(data[1:5, ], 4, TRUE)
dim(result$item)
```

    ## [1] 5 4

``` r

dim(result$distance)
```

    ## [1] 5 4

``` r

ann$markDeleted(1)
ann$size()
```

    ## [1] 150

``` r

path <- tempfile(fileext = ".hnsw")
ann$save(path)
loaded <- new(HnswL2, item_dim, path)
loaded$setEf(50)
loaded$size()
```

    ## [1] 150

``` r

unlink(path)
```
