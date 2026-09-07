# Find approximate nearest neighbors

Build an HNSW index and find `k` neighbors for every item in `X`.

## Usage

``` r
hnsw_knn(
  X,
  k = 10,
  distance = "euclidean",
  M = 16,
  ef_construction = 200,
  ef = 10,
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

  Numeric matrix with one item per row, or per column when
  `byrow = FALSE`.

- k:

  Number of neighbors per item. A positive whole number no larger than
  the number of items in `X`.

- distance:

  Distance to use:

  - `"euclidean"`: Euclidean distance.

  - `"l2"`: squared Euclidean distance.

  - `"cosine"`: one minus cosine similarity.

  - `"ip"`: one minus inner product, `1 - sum(a * b)`; can be negative.

- M:

  Number of graph links per item during construction. Larger values
  improve connectivity and use more memory. Must be between 2 and 10000.

- ef_construction:

  Candidate-list size during construction. Larger values improve index
  quality and increase build time. Must be a positive whole number;
  raised to at least `k`.

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

- random_seed:

  Seed for hnswlib's index construction. R's
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect. Use
  serial construction for repeatable builds; parallel insertion order
  can vary with a fixed seed.

## Value

A list with matrices `idx` (one-based neighbor indices) and `dist`
(distances), ordered by increasing distance for each query. Both are
`n × k` when `byrow = TRUE`, or `k × n` otherwise, where `n` is the
number of items in `X`.

## Details

If you are searching for neighbors within `X`, don't assume that the
first neighbor is always the item itself: approximate search can miss
self-matches. With inner-product distance, an item need not be its own
nearest neighbor even in an exact search.

Coordinates are stored in single precision. See
[RcppHnsw-package](https://jlmelville.github.io/rcpphnsw/reference/RcppHnsw-package.md)
for numeric limits. For help choosing the parameters, see the [parameter
guide](https://jlmelville.github.io/rcpphnsw/articles/parameters-metrics-data-layout.html)
which also covers data layout.

## References

Malkov, Y. A., & Yashunin, D. A. (2020). Efficient and robust
approximate nearest neighbor search using Hierarchical Navigable Small
World graphs. *IEEE Transactions on Pattern Analysis and Machine
Intelligence*, 42(4), 824-836.
[doi:10.1109/TPAMI.2018.2889473](https://doi.org/10.1109/TPAMI.2018.2889473)
.

## See also

[`hnsw_build()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_build.md)
and
[`hnsw_search()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_search.md)
to build once and search new data.

## Examples

``` r
irism <- as.matrix(iris[, -5])
neighbors <- hnsw_knn(irism, k = 4)
neighbors$idx[1:5, ]
#>      [,1] [,2] [,3] [,4]
#> [1,]    1   18    5   29
#> [2,]    2   13   46   35
#> [3,]    3   48    4    7
#> [4,]    4   48   30   31
#> [5,]    5    1   38   18
neighbors$dist[1:5, ]
#>      [,1]      [,2]      [,3]      [,4]
#> [1,]    0 0.1000000 0.1414212 0.1414212
#> [2,]    0 0.1414213 0.1414213 0.1414213
#> [3,]    0 0.1414213 0.2449490 0.2645752
#> [4,]    0 0.1414215 0.1732051 0.2236071
#> [5,]    0 0.1414212 0.1414213 0.1732050
```
