# RcppHNSW 0.4.9000

## New features

* For high-dimensional data, there can be a noticeable CPU overhead in copying 
data out of the non-contiguous memory regions when row-wise data is used. The
following new methods have been added to the class-based API for batch
index-creation and retrieval using matrices where the items are stored
column-wise, i.e. there is one item per column: `addItemsCol` and `getAllNNsCol`
`getAllNNsListCol`. These are the column-based equivalents of `addItems`,
`getAllNNs` and `getAllNNsList`, respectively. Note that the returned nearest
neighbor data from `getAllNNsCol` and `getAllNNsListCol` are *also* stored
by column, i.e. the matrices have dimensions `k x n` where `k` is the number
of neighbors, and `n` the number of items in the data being searched.

## Bug fixes and minor improvements

* The `progress` parameter in the functional interface no longer does anything. When 
`verbose = TRUE`, a progress bar is no longer shown.

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
