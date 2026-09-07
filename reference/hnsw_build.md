# Build a nearest neighbor index

Build an HNSW index that you can keep and query with
[`hnsw_search()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_search.md).

## Usage

``` r
hnsw_build(
  X,
  distance = "euclidean",
  M = 16,
  ef = 200,
  verbose = FALSE,
  progress = "bar",
  n_threads = 0,
  grain_size = 1,
  byrow = TRUE,
  random_seed = 100
)
```

## Arguments

- X:

  Numeric matrix to index, with one item per row, or per column when
  `byrow = FALSE`.

- distance:

  Distance to use:

  - `"euclidean"`: Euclidean distance.

  - `"l2"`: squared Euclidean distance.

  - `"cosine"`: one minus cosine similarity.

  - `"ip"`: one minus inner product, `1 - sum(a * b)`; can be negative.

- M:

  Number of graph links per item during construction. Larger values
  improve connectivity and use more memory. Must be between 2 and 10000.

- ef:

  Candidate-list size during construction. Larger values improve index
  quality and increase build time. Must be a positive whole number and
  is not capped at the dataset size. This is the `ef_construction`
  parameter of
  [`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md).

- verbose:

  If `TRUE`, log messages to the console.

- progress:

  Unused; retained for compatibility.

- n_threads:

  Maximum number of threads for batch insertion or search. Zero (the
  default) and one run serially.

- grain_size:

  Minimum number of items per thread. Larger values limit threading
  overhead for small batches. Zero is treated as one.

- byrow:

  If `TRUE`, items are rows of `X`; otherwise, they are columns.

- random_seed:

  Seed for hnswlib's index construction. R's
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect. Use
  serial construction for repeatable builds; parallel insertion order
  can vary with a fixed seed.

## Value

An `HnswEuclidean`, `HnswL2`, `HnswCosine`, or `HnswIp` index, according
to `distance`. Labels are one-based and follow the order in `X`.

## Details

Coordinates are stored in single precision; see
[RcppHnsw-package](https://jlmelville.github.io/rcpphnsw/reference/RcppHnsw-package.md)
for numeric limits. If you want to add more items, resize, or save the
index, see the [Module
guide](https://jlmelville.github.io/rcpphnsw/articles/module-api-index-lifecycle.html)
for the available methods.

## See also

[`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md)
to build and search in one call.

## Examples

``` r
irism <- as.matrix(iris[, -5])
ann <- hnsw_build(irism[1:100, ])
neighbors <- hnsw_search(irism[101:150, ], ann, k = 5)
```
