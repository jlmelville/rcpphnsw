## RcppHNSW

Rcpp bindings for [HNSW](https://github.com/nmslib/hnsw).

### Status

Don't get too excited. Barely any of HNSW has been exposed so far. But it can
calculate nearest neighbors via Euclidean distance.

### HNSW

HNSW is a header-only C++ library for finding approximate nearest neighbors
(ANN) via Hierarchical Navigable Small Worlds
[(Yashunin and Malkov, 2016)](https://arxiv.org/abs/1603.09320). 
It is part of the [nmslib](https://github.com/nmslib/nmslib]) project. 

### RcppHNSW

An R package that interfaces with HNSW, taking enormous amounts of inspiration from 
[Dirk EddelBuettel](https://github.com/eddelbuettel)'s 
[RcppAnnoy](https://github.com/eddelbuettel/rcppannoy) package which did the same for
the [Annoy](https://github.com/spotify/annoy) ANN C++ library. 

One difference is that I use
[roxygen2](https://cran.r-project.org/package=roxygen2) to generate the man
pages. The `NAMESPACE` is still built manually, however (I don't believe you can
`export` the classes currently).

### Installing

```R
devtools::install_github("jlmelville/RcppHNSW")
```

### Example

```R
data <- as.matrix(iris[, -5])

# Create a new index using the L2 (squared Euclidean) distance
# nr and nc are the number of rows and columns of the data to be added, respectively
# ef and M determines speed vs accuracy trade off
ann <- methods::new(RcppHNSW::HnswL2, nc = ncol(data), nr = nrow(data), 
                    ef = 200, M = 16)

# Add item to index
for (i in 1:nr) {
  ann$addItem(data[1, ])
}

# Find 4 nearest neighbors of row 1
# indexes are in res$item, distances in res$distance
res <- ann$getNNsList(data[1, ], k = 4, include_distances = TRUE)

# function interface returns results for all rows in nr x k matrices
all_knn <- get_knn(data, k = 4)
```

### Note

I had to add a non-portable flag to `PKG_CPPFLAGS` (`-march=native`). This should be
the only status warning in `R CMD check`.


### License

[GPL-3 or later](https://www.gnu.org/licenses/gpl-3.0.en.html).
