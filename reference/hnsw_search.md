# Search an hnswlib nearest neighbor index

Search an hnswlib nearest neighbor index

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

  A numeric matrix of data to search for neighbors. If `byrow = TRUE`
  (the default) then each row of `X` is an item to be searched.
  Otherwise, each item should be stored in the columns of `X`.

- ann:

  an instance of an `HnswEuclidean`, `HnswL2`, `HnswCosine` or `HnswIp`
  class.

- k:

  Positive whole number of neighbors to return. It cannot exceed the
  active (not deleted) item count. `ann$size()` reports the total number
  added and therefore may be larger than the active count after
  deletion.

- ef:

  Size of the dynamic list used during search. Higher values lead to
  improved recall at the expense of longer search time. Must be
  positive; the effective value is at least `k` and is not bounded by
  index size. Typical values are `100 - 2000`.

- verbose:

  If `TRUE`, log messages to the console.

- progress:

  defunct and has no effect.

- n_threads:

  Maximum number of threads to use. Zero and one both select serial
  execution. For larger values, the exact number is determined by
  `grain_size` and the amount of work.

- grain_size:

  Minimum amount of work to do (items in `X` to search) per thread. Zero
  is treated as one. If the number of items in `X` isn't sufficient,
  then fewer than `n_threads` will be used. This is useful in cases
  where the overhead of context switching with too many threads
  outweighs the gains due to parallelism.

- byrow:

  If `TRUE` (the default), this indicates that the items to be searched
  in `X` are stored in each row of `X`. Otherwise, the items are stored
  in the columns of `X`. Storing items in each column reduces the
  overhead of copying data to a form that can be searched by the `hnsw`
  library. Note that if `byrow = FALSE`, any matrices returned from this
  function will also store the items by column.

## Value

a list containing:

- `idx` a matrix containing the nearest neighbor indices.

- `dist` a matrix containing the nearest neighbor distances.

The dimensions of the matrices respect the storage (row or column-based)
of `X` as indicated by the `byrow` parameter. If `byrow = TRUE` (the
default) each row of `idx` and `dist` contain the neighbor information
for the item passed in the equivalent row of `X`, i.e. the dimensions
are `n x k` where `n` is the number of items in `X`. If `byrow = FALSE`,
then each column of `idx` and `dist` contain the neighbor information
for the item passed in the equivalent column of `X`, i.e. the dimensions
are `k x n`.

## Numeric data and index mutation

The package rejects non-finite or out-of-range query coordinates and
cosine queries with zero norm after conversion to single precision.
Search is approximate, and under inner-product distance an item need not
be its own nearest neighbor.

This function updates `ann` in place: the effective `ef`, `n_threads`,
and `grain_size` settings remain on the external index after the call.

## Examples

``` r
irism <- as.matrix(iris[, -5])
ann <- hnsw_build(irism)
iris_nn <- hnsw_search(irism, ann, k = 5)
```
