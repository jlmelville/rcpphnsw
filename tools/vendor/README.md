# Maintain the vendored hnswlib headers

`inst/include/` contains the headers compiled by RcppHNSW. The hnswlib headers are reconstructed
from the archive pinned in `upstream.yml`, checked against `upstream-files.sha256`, and modified
by the patches listed in `patches/series`. The separate `pforr/` directory is left untouched.

## Verify

```sh
tools/vendor/verify-hnswlib.sh
```

This downloads the pinned archive into temporary storage and checks that reconstruction matches
all seven committed headers byte for byte. For an offline check, supply a local archive:

```sh
tools/vendor/verify-hnswlib.sh --archive hnswlib-v0.9.0.tar.gz
```

The `native-maintenance` workflow runs this verifier when the headers or maintenance tools change.
To test drift detection on a disposable copy, pass `--materialized-dir DIR`.

## Refresh after editing patches

```sh
tools/vendor/refresh-hnswlib.sh --archive hnswlib-v0.9.0.tar.gz
git diff -- inst/include
```

Refresh overwrites the materialized hnswlib headers. Run it after updating the patch queue, then
review the diff. Package installation uses the committed headers without downloading or patching.
The manifest, checksums, patches, and scripts ship in source packages as provenance.

## Patch policy

The queue accounts for every difference from the pinned upstream header set:

1. `0001` retains the package's existing checked external-label write.
2. `0002` retains the package's existing error-stream integration.
3. `0003` serializes access to the shared random-level generator. The focused TSan
   reproducer is `../diagnostics/run-hnswlib-rng-tsan.sh`.
4. `0004` snapshots the entry point while hnswlib's existing global lock is held,
   removing the second data race exposed by the same reproducer.
5. `0005` validates serialized record offsets before subtraction and rejects a
   coordinate payload width that does not exactly match the requested space.
6. `0006` makes raw index saves report stream open, write, flush, and close
   failures instead of returning after an incomplete or absent write.
7. `0007` adds the prominent dated notice required for the modified vendored
   header and points recipients to the installed provenance record and exact
   source patch series.
8. `0008` preserves the saved capacity when an empty index is loaded without an
   explicit capacity.
9. `0009` rejects zero-capacity resize before any storage is reallocated.
10. `0010` gives initial visited-list pool entries temporary ownership so a
    failed pool insertion does not leak them during construction or resize.

The queue must account for every difference from upstream. Keep patches focused; avoid unrelated
formatting changes that make future updates harder to review.

If the distinct-label insertion diagnostic reveals another independent data race, serialize
package construction and reassess the concurrency contract before extending the queue.
