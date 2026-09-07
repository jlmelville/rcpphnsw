# RcppHNSW

[![R-CMD-check](https://github.com/jlmelville/rcpphnsw/workflows/R-CMD-check/badge.svg)](https://github.com/jlmelville/rcpphnsw/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/jlmelville/rcpphnsw/master.svg)](https://codecov.io/github/jlmelville/rcpphnsw?branch=master)
[![CRAN Status
Badge](https://www.r-pkg.org/badges/version/RcppHNSW)](https://cran.r-project.org/package=RcppHNSW)

An R interface to [hnswlib](https://github.com/nmslib/hnswlib), a C++
library for approximate nearest neighbor search using hierarchical
navigable small-world graphs [(Malkov and Yashunin,
2020)](https://doi.org/10.1109/TPAMI.2018.2889473). It supports
Euclidean, squared Euclidean, cosine, and inner-product distances.

## Installation

From CRAN:

``` r

install.packages("RcppHNSW")
```

Or, for the development version, use [pak](https://pak.r-lib.org/):

``` r

pak::pak("jlmelville/rcpphnsw")
```

## Find nearest neighbors

To find neighbors within a dataset, use
[`hnsw_knn()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_knn.md).
Here we’ll ask for four neighbors for each item in `iris`, leaving out
the species column:

``` r

library(RcppHNSW)

irism <- as.matrix(iris[, -5])
neighbors <- hnsw_knn(irism, k = 4)
neighbors$idx[1:5, ]
neighbors$dist[1:5, ]
```

The result is a list containing two matrices: `idx` gives the neighbors’
row numbers, and `dist` gives their Euclidean distances. Each row
corresponds to an item in `irism`, with one column per neighbor.

If you want to find neighbors for new data, you can keep the index
around and query it with
[`hnsw_search()`](https://jlmelville.github.io/rcpphnsw/reference/hnsw_search.md).
Here we’ll use the first 100 items to build the index, then find
neighbors for the remaining 50:

``` r

ann <- hnsw_build(irism[1:100, ])
neighbors <- hnsw_search(irism[101:150, ], ann, k = 4)
```

If the results are a bit too approximate, try increasing the search
parameter `ef`. This gives the search more work to do, so it also takes
longer. See the [parameter
guide](https://jlmelville.github.io/rcpphnsw/articles/parameters-metrics-data-layout.html)
for tuning, distance choices, and column-oriented data.

## Documentation

- [Function
  reference](https://jlmelville.github.io/rcpphnsw/reference/index.html)
- [Parameters, metrics, and data
  layout](https://jlmelville.github.io/rcpphnsw/articles/parameters-metrics-data-layout.html)
- [Module API and index
  lifecycle](https://jlmelville.github.io/rcpphnsw/articles/module-api-index-lifecycle.html):
  add items, delete items, resize, and save an index.
- [NEWS](https://jlmelville.github.io/rcpphnsw/news/index.html)

## Project

The inspiration for this package is
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy), an R interface
to [Annoy](https://github.com/spotify/annoy). Source code and issue
tracking are on [GitHub](https://github.com/jlmelville/rcpphnsw).

Licensed under [GPL-3 or
later](https://www.gnu.org/licenses/gpl-3.0.en.html).
