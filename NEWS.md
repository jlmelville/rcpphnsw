# RcppHNSW (development version)

## New features

* New class: `HnswEuclidean`. This uses Euclidean distances internally and will
be returned from `hnsw_build` when `distance = "euclidean"` is specified. This
fixes an issue where if you created an index with `hnsw_build` and 
`distance = "euclidean"` (the default), then after saving, you would be unable
to reload the index and have it find Euclidean distances. You would have to
create it as an `HsnwL2` object and take the square root of the distances
yourself (<https://github.com/jlmelville/rcpphnsw/issues/21>).

# RcppHNSW 0.6.0

## New features

* Updated hnswlib to [version 0.8.0](https://github.com/nmslib/hnswlib/releases/tag/v0.8.0).

# RcppHNSW 0.5.0

## New features

* Updated hnswlib to [version 0.7.0](https://github.com/nmslib/hnswlib/releases/tag/v0.7.0).
Note that I made some very minor changes to the code to silence some compiler warnings. These
changes have been submitted up-stream to the hnswlib project.
* For high-dimensional data, there can be a noticeable CPU overhead in copying data out of the
non-contiguous memory regions when row-wise data is used. If you wish to provide data where each
*column* of the input matrix contains an item to be indexed/search then see the following additions
to the API:
  * For the class-based API: `addItemsCol`, `getAllNNsCol` and `getAllNNsListCol` are the
  column-based equivalents of `addItems`, `getAllNNs` and `getAllNNsList`, respectively. Note that
  the returned nearest neighbor data from `getAllNNsCol` and `getAllNNsListCol` are *also* stored
  by column, i.e. the matrices have dimensions `k x n` where `k` is the number of neighbors, and
  `n` the number of items in the data being searched.
  * For the function-based API, a new parameter `byrow` has been added to `hnsw_knn`, `hnsw_build`
  and `hnsw_search`. By default this is set to `TRUE` and indicates that the items in the input
  matrix are found in each row. To pass column-stored items, set `byrow = FALSE`. Any matrices
  returned by `hnsw_search` and `hnsw_knn` will now follow the convention provided by the value of
  `byrow`: i.e. if  `byrow = FALSE`, the matrices contain nearest neighbor information in each
  column.
* new method: `getItems`, which returns a matrix of the data vectors in the index with the
  specified integer identifiers. From a feature request made by [d4tum](https://github.com/d4tum)
  (<https://github.com/jlmelville/rcpphnsw/issues/18>).

## Bug fixes and minor improvements

* The `progress` parameter in the functional interface no longer does anything. When
`verbose = TRUE`, a progress bar is no longer shown.
* Due to a breaking change in roxygen2 7.0.0, there was a missing package alias in the
documentation.

# RcppHNSW 0.4.1

## Bug fixes and minor improvements

* Rolled back to
[hnswlib v0.4.0](https://github.com/nmslib/hnswlib/releases/tag/v0.4.0)
due to valgrind problems in v0.6.2

# RcppHNSW 0.4.0

## New features

* Updated hnswlib to 
[version 0.6.2](https://github.com/nmslib/hnswlib/releases/tag/v0.6.2).

## Bug fixes and minor improvements

* Minor future-proofing of licensing: RcppHNSW is now GPLv3 or later, rather
than GPLv3 only.

# RcppHNSW 0.3.0

## New features

* Multi-threading support is now available. Use the `setNumThreads` method if 
using the object-based API, and the `n_threads` parameter in the `hnsw_*` 
function API. For finer control, a `setGrainSize` and `grain_size` option is
also available in the object and function interface respectively. Thank you
to [Dmitriy Selivanov](https://github.com/dselivanov) for a lot of the work on
this.
* Updated hnswlib to 
[version 0.4.0](https://github.com/nmslib/hnswlib/releases/tag/v0.4.0).

## Bug fixes and minor improvements

* Setting `verbose = TRUE` now has incurs substantially less computational 
overhead associated with calculating the progress bar. Thank you to 
[Samuel Granjeaud](https://github.com/SamGG) for spotting the problem and coming
up with various solutions.
* New parameter: `progress`. By default this is set to `"bar"` and will show the
progress bar when `verbose = TRUE`. If you want a more terse output, set
`progress = NULL`. `progress = NULL` will eventually be the default setting:
for now, `verbose = TRUE` will get you the progress bar by default for backwards
compatibility.
* No progress bar will be shown if you have less than 50 items to process.

# RcppHNSW 0.2.0

## New features

* Updated hnswlib to <https://github.com/nmslib/hnswlib/commit/c5c38f0> 
(20 September 2019).
* A new method, `markDeleted`, that will remove an object from being retrieved
from the index.
* A new method, `resizeIndex`, that allows the index to be increased without 
having to save and reload the index.
* A new method, `size` is available for the index objects and reports the
number of items added to the index.


## Bug fixes and minor improvements

* `hnsw_search` would `stop` if the number of rows in the input matrix was 
smaller than `k`. This check has been removed. Note that the correct behavior is
to ensure that `k` is smaller than or equal to `index$size()` where `index` is
the index you are searching. Because the `size()` method is new to this version,
to preserve compatibility with old indexes, this check *hasn't* been added to
`hnsw_search`. If this matters to you, manually compare `index$size()` with `k`
before running `hnsw_search`. An error will be thrown if `k` neighbors can't be
found in the index. Thank you to [Yuxing Liao](https://github.com/yxngl) for 
spotting this and the pull request to remove the check.

# RcppHNSW 0.1.0

Initial release.
