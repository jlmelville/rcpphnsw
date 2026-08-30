# Build an hnswlib nearest neighbor index

Build an hnswlib nearest neighbor index

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

  A numeric matrix of data to search for neighbors. If `byrow = TRUE`
  (the default) then each row of `X` is an item to be searched.
  Otherwise, each item should be stored in the columns of `X`.

- distance:

  Type of distance to calculate. One of:

  - `"l2"` Squared L2, i.e. squared Euclidean.

  - `"euclidean"` Euclidean.

  - `"cosine"` One minus cosine similarity.

  - `"ip"` One minus inner product: `1 - sum(a * b)`. Values can be
    negative and need not satisfy metric properties.

- M:

  Controls the number of bi-directional links created for each element
  during index construction. Higher values lead to better results at the
  expense of memory consumption. Typical values are `2 - 100`, but for
  most datasets a range of `12 - 48` is suitable. Can't be smaller than
  2.

- ef:

  Size of the dynamic list used during construction. A larger value
  means a better quality index, but increases build time. Must be a
  positive whole number and is not bounded by the size of the dataset.

- verbose:

  If `TRUE`, log messages to the console.

- progress:

  defunct and has no effect.

- n_threads:

  Maximum number of threads to use. Zero and one both select serial
  execution. For larger values, the exact number is determined by
  `grain_size` and the amount of work.

- grain_size:

  Minimum number of items in `X` to add per thread. Zero is treated as
  one. If the number of items in `X` isn't sufficient, then fewer than
  `n_threads` will be used. This is useful in cases where the overhead
  of context switching with too many threads outweighs the gains due to
  parallelism.

- byrow:

  If `TRUE` (the default), this indicates that the items in `X` to be
  indexed are stored in each row. Otherwise, the items are stored in the
  columns of `X`. Storing items in each column reduces the overhead of
  copying data to a form that can be indexed by the `hnsw` library.

- random_seed:

  Seed passed to hnswlib for index construction. The default, `100`, is
  the underlying hnswlib default. This seed belongs to hnswlib: calling
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not affect
  index construction.

## Value

an instance of an `HnswEuclidean`, `HnswL2`, `HnswCosine` or `HnswIp`
class.

## Numeric data and reproducibility

Coordinates are stored as single-precision floating-point values. The
package rejects non-finite or out-of-range coordinates and, for cosine
distance, vectors with zero norm after conversion. Parallel construction
may be nondeterministic even for a fixed `random_seed`. Zero or one
thread uses serial construction. R's random seed is not used by hnswlib.

## Examples

``` r
irism <- as.matrix(iris[, -5])
ann <- hnsw_build(irism)
iris_nn <- hnsw_search(irism, ann, k = 5)
```
