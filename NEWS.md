# RcppHNSW 0.2.0

## New features

* Updated hnswlib to <https://github.com/nmslib/hnswlib/commit/15e64f69> 
(1 August 2019).
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
