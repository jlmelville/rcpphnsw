# RcppHNSW

[![R-CMD-check](https://github.com/jlmelville/rcpphnsw/workflows/R-CMD-check/badge.svg)](https://github.com/jlmelville/rcpphnsw/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/jlmelville/rcpphnsw/master.svg)](https://codecov.io/github/jlmelville/rcpphnsw?branch=master)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)

Rcpp bindings for [hnswlib](https://github.com/nmslib/hnswlib), a header-only
C++ library for finding approximate nearest neighbors via Hierarchical
Navigable Small Worlds
[(Malkov and Yashunin, 2020)](https://doi.org/10.1109/TPAMI.2018.2889473).

## Installing

From CRAN:

```r
install.packages("RcppHNSW")
```

Development versions from GitHub:

```r
pak::pak("jlmelville/RcppHNSW")
```

## Choose an interface

| Task | Interface | Reference |
| --- | --- | --- |
| Find neighbors within one matrix | `hnsw_knn()` | [`hnsw_knn` help](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_knn.Rd) |
| Build an index for later searches | `hnsw_build()` | [`hnsw_build` help](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_build.Rd) |
| Search an existing index | `hnsw_search()` | [`hnsw_search` help](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_search.Rd) |
| Control index construction and lifecycle directly | `HnswL2`, `HnswEuclidean`, `HnswCosine`, or `HnswIp` | [Module API and index lifecycle](https://github.com/jlmelville/rcpphnsw/blob/master/vignettes/articles/module-api-index-lifecycle.Rmd) |

## Quick start

```r
irism <- as.matrix(iris[, -5])

# returns neighbors and distances for all rows in n x k matrices
iris_knn <- RcppHNSW::hnsw_knn(irism, k = 4, distance = "euclidean")
dim(iris_knn$idx)
```

The process can also be split into two steps, so you can build with one set of
data and search with another. See [`hnsw_build()`](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_build.Rd)
and [`hnsw_search()`](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_search.Rd).

## Distance metrics

| `distance` | Distance calculated |
| --- | --- |
| `"l2"` | Squared L2, i.e. squared Euclidean. |
| `"euclidean"` | Euclidean. |
| `"cosine"` | Cosine. |
| `"ip"` | Inner product: `1 - sum(ai * bi)`, i.e. the cosine distance where the vectors are not normalized. This can lead to negative distances and other non-metric behavior. |

## Important contracts

- HNSW search is approximate. For L2, Euclidean, and cosine distance, an item
  queried against its source data has distance zero from itself, but it can be
  omitted when recall is insufficient. Under inner-product distance, an item
  need not be its own nearest neighbor.
- All input coordinates are converted to and stored as single-precision
  floats. They must be finite and representable in that format. Cosine items
  and queries must have a positive finite norm after conversion.
- Items are stored by row by default. Column-oriented data and the corresponding
  result shapes are described in the
  [parameters, metrics, and data-layout guide](https://github.com/jlmelville/rcpphnsw/blob/master/vignettes/articles/parameters-metrics-data-layout.Rmd).
- `random_seed` belongs to hnswlib: calling `set.seed()` does not affect index
  construction. Parallel construction may be nondeterministic even with a
  fixed `random_seed`.
- Module labels are one-based and follow insertion order. `k` must be positive
  and cannot exceed the active (not deleted) item count. If an exception
  escapes after insertion has begun, the index becomes unusable and must be
  discarded and rebuilt or reloaded.
- Raw files created by the Module `ann$save()` method are operational hnswlib
  checkpoints, not a stable cross-version or cross-platform archive. See the
  [Module API and index-lifecycle guide](https://github.com/jlmelville/rcpphnsw/blob/master/vignettes/articles/module-api-index-lifecycle.Rmd).

## Documentation

Function reference:

- [`hnsw_knn()`](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_knn.Rd)
- [`hnsw_build()`](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_build.Rd)
- [`hnsw_search()`](https://github.com/jlmelville/rcpphnsw/blob/master/man/hnsw_search.Rd)

Guides:

- [Parameters, metrics, and data layout](https://github.com/jlmelville/rcpphnsw/blob/master/vignettes/articles/parameters-metrics-data-layout.Rmd)
- [Module API and index lifecycle](https://github.com/jlmelville/rcpphnsw/blob/master/vignettes/articles/module-api-index-lifecycle.Rmd)

## Project

RcppHNSW takes enormous amounts of inspiration from
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy), which provides an R
interface to the [Annoy](https://github.com/spotify/annoy) library.

See [NEWS](https://github.com/jlmelville/rcpphnsw/blob/master/NEWS.md) for
release and development changes. Source and issue tracking are on
[GitHub](https://github.com/jlmelville/rcpphnsw).

## License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
