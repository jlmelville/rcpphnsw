# hnswlib vendor refresh

`inst/include/` contains the materialized headers compiled by RcppHNSW. The
hnswlib portion is reconstructed from the exact v0.9.0 archive described by
`upstream.yml`, checked against `upstream-files.sha256`, and modified only by
the ordered files in `patches/series`. The separate `pforr/` directory is not
part of hnswlib and is left untouched.

Package installation never downloads or patches source. These scripts are for
maintainers updating or auditing the committed materialized headers.

Verify the current tree, downloading the pinned archive into a temporary
directory:

```sh
tools/vendor/verify-hnswlib.sh
```

An already downloaded archive can be supplied for an offline audit:

```sh
tools/vendor/verify-hnswlib.sh --archive hnswlib-v0.9.0.tar.gz
```

After intentionally changing the patch series, refresh the materialized
headers and then inspect the repository diff:

```sh
tools/vendor/refresh-hnswlib.sh --archive hnswlib-v0.9.0.tar.gz
git diff -- inst/include
```

The verifier accepts `--materialized-dir DIR` so a disposable copy can be used
to prove that unexplained drift fails without modifying the repository.

## Patch policy

The queue accounts for every difference from the pinned upstream header set:

1. `0001` retains the package's existing checked external-label write.
2. `0002` retains the package's existing error-stream integration.
3. `0003` serializes access to the shared random-level generator. The focused TSan
   reproducer is `../diagnostics/run-hnswlib-rng-tsan.sh`.
4. `0004` snapshots the entry point while hnswlib's existing global lock is held,
   removing the second data race exposed by the same reproducer.

The last two patches are deliberately independent and suitable for an upstream
report. They do not change the public API, index format, graph heuristics, or
distance calculations. If the supported distinct-new-label insertion workload
reveals another independent data race, do not extend this queue one field at a
time; serialize package construction and reassess the concurrency contract.
