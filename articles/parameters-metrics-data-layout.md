# Parameters, metrics, and data layout

## Distance

| `distance` | Distance calculated |
|----|----|
| `"l2"` | Squared L2, i.e. squared Euclidean. |
| `"euclidean"` | Euclidean. |
| `"cosine"` | One minus cosine similarity. |
| `"ip"` | One minus inner product: `1 - sum(a * b)`. Values can be negative and need not satisfy metric properties. |

Coordinates are stored as single-precision floats. The package rejects
non-finite or out-of-range coordinates and cosine vectors with zero norm
after conversion.

## Data layout and result shape

If `byrow = TRUE` (the default), the items to be processed in `X` are
stored in each row of `X`. Otherwise, the items are stored in the
columns of `X`. Storing items in each column reduces the overhead of
copying data to a form that can be used by the `hnsw` library.

The dimensions of the matrices respect the storage (row or column-based)
of `X` as indicated by the `byrow` parameter.

| `byrow` | Items in `X` | Items in `idx` and `dist` | Result dimensions |
|---------|--------------|---------------------------|-------------------|
| `TRUE`  | Rows         | Rows                      | `n x k`           |
| `FALSE` | Columns      | Columns                   | `k x n`           |

## Construction and search parameters

Some details on the parameters used for index construction and search,
based on the [hnswlib algorithm
parameters](https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md):

| Parameter | Used by |
|----|----|
| `M` | Index construction in [`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md) and [`hnsw_build()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_build.md). |
| `ef_construction` | Index construction in [`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md). |
| `ef` | Index construction in [`hnsw_build()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_build.md) and search in [`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md) and [`hnsw_search()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_search.md). |

- `M` Controls the number of bi-directional links created for each
  element during index construction. Higher values lead to better
  results at the expense of memory consumption, which is around
  `M * 8-10` bytes per stored element. High intrinsic dimensionalities
  will require higher values of `M`. A range of `2 - 100` is typical,
  but `12 - 48` is ok for most use cases.
- `ef_construction` Size of the dynamic list used during construction. A
  larger value means a better quality index, but increases build time.
  It must be a positive whole number, but is not bounded by the size of
  the dataset. A typical range is `100 - 2000`. Beyond a certain point,
  increasing `ef_construction` has no effect. A sufficient value of
  `ef_construction` can be determined by searching with
  `ef = ef_construction`, and ensuring that the recall is at least 0.9.
- `ef` Size of the dynamic list used during index search. Can differ
  from `ef_construction`. The effective value is at least `k`, and it is
  not bounded by the number of elements in the index.

## Threads and grain size

`n_threads` Maximum number of threads to use. Zero and one both select
serial execution. For larger values, the exact number is determined by
`grain_size` and the amount of work.

`grain_size` Minimum number of items in `X` to add or search per thread.
Zero is treated as one. If the number of items in `X` isn’t sufficient,
then fewer than `n_threads` will be used. This is useful in cases where
the overhead of context switching with too many threads outweighs the
gains due to parallelism.

## Random seed, reproducibility, and approximation

`random_seed` Seed passed to hnswlib for index construction. The
default, `100`, is the underlying hnswlib default. This seed belongs to
hnswlib: calling [`set.seed()`](https://rdrr.io/r/base/Random.html) does
not affect index construction.

Parallel index construction may be nondeterministic even for a fixed
`random_seed`. Use serial construction for repeatable reconstruction.
Save the constructed index to reuse the exact graph.

HNSW search is approximate. For L2, Euclidean, and cosine distance, an
item queried against its source data has distance zero from itself, but
it can be omitted when recall is insufficient. Under inner-product
distance, an item need not be its own nearest neighbor.

## Function example

``` r

irism <- as.matrix(iris[, -5])

# function interface returns results for all rows in n x k matrices
all_knn <- RcppHNSW::hnsw_knn(irism, k = 4, distance = "l2")
dim(all_knn$idx)
```

    ## [1] 150   4

``` r

# other distance options: "euclidean", "cosine" and "ip" (inner product distance)

# for high-dimensional data you may see a speed-up if you store the data
# where each *column* is an item to be indexed and searched. Set byrow = FALSE
# for this.
# Admittedly, the iris dataset is *not* high-dimensional
iris_by_col <- t(irism)
all_knn_by_col <- RcppHNSW::hnsw_knn(
  iris_by_col,
  k = 4,
  distance = "l2",
  byrow = FALSE
)
dim(all_knn_by_col$idx)
```

    ## [1]   4 150

``` r

# process can be split into two steps, so you can build with one set of data
# and search with another
ann <- RcppHNSW::hnsw_build(irism[1:100, ])
iris_nn <- RcppHNSW::hnsw_search(irism[101:150, ], ann, k = 5)
dim(iris_nn$idx)
```

    ## [1] 50  5
