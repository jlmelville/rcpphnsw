# Search a nearest neighbor index

Find neighbors in an existing index for each query in `X`.

## Usage

``` r
hnsw_search(
  X,
  ann,
  k,
  ef = 10,
  verbose = FALSE,
  progress = "bar",
  n_threads = 0,
  grain_size = 1,
  byrow = TRUE
)
```

## Arguments

- X:

  Numeric query matrix with the same number of dimensions as the indexed
  data. Each row is a query, or each column when `byrow = FALSE`.

- ann:

  An index from
  [`hnsw_build()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_build.md)
  or a Module constructor.

- k:

  Number of neighbors per query. A positive whole number no larger than
  the number of undeleted items. `ann$size()` includes deleted items.

- ef:

  Candidate-list size during search. Increase it to improve recall at
  the cost of search time. Must be a positive whole number; the
  effective value is at least `k` and is not capped at the index size.

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

  If `TRUE`, items are rows of `X`; otherwise, they are columns. Results
  follow the same orientation.

## Value

A list with matrices `idx` (one-based labels in `ann`) and `dist`
(distances), ordered by increasing distance for each query. Both are
`n × k` when `byrow = TRUE`, or `k × n` otherwise, where `n` is the
number of queries in `X`.

## Details

The index's distance measure is used for all queries. Search is
approximate; with inner-product distance, an item need not be its own
nearest neighbor. Coordinates are stored in single precision; see
[RcppHnsw-package](https://jlmelville.github.io/rcpphnsw/reference/RcppHnsw-package.md)
for numeric limits.

This call also updates `ann`'s `ef`, `n_threads`, and `grain_size`
settings. If you then use the Module methods directly, they will use
these settings.

## See also

[`hnsw_build()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_build.md)
to create an index,
[`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md)
to build and search in one call.

## Examples

``` r
irism <- as.matrix(iris[, -5])
ann <- hnsw_build(irism[1:100, ])
neighbors <- hnsw_search(irism[101:150, ], ann, k = 5)
```
