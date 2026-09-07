# Native diagnostics

The `native-maintenance` GitHub Actions workflow runs these checks when the hnswlib headers or
maintenance tools change. To run them locally, use the commands below from the repository root.
Both runners build in temporary directories and remove their binaries on exit. `.Rbuildignore`
excludes this directory from source packages.

## Concurrent insertion: ThreadSanitizer

```sh
tools/diagnostics/run-hnswlib-rng-tsan.sh
```

The diagnostic seeds an index with one item, then inserts distinct labels concurrently. It checks
the repairs to random-level generation and the entry-point snapshot. Each round must increase the
index's maximum level. A race report or failed check exits nonzero.

The optional arguments are rounds, threads, points per thread, `M`, and seed:

```sh
tools/diagnostics/run-hnswlib-rng-tsan.sh 8 12 256 16 101
```

The runner defaults to `clang++` and honors `CXX` and `TSAN_OPTIONS`. Its default options disable
the separate lock-order detector to focus on data races.

If this workload reveals another independent data race, serialize package construction and reassess
the concurrency contract before adding further vendor patches.

## Resize failures: AddressSanitizer

On Linux with a linker that supports `--wrap`:

```sh
tools/diagnostics/run-hnswlib-resize-asan.sh
```

The diagnostic injects each allocation failure reached during growth and shrinkage. It checks
stored-item search, safe destruction, balanced C++ `new` allocations, and zero-capacity rejection.
AddressSanitizer reports and failed checks exit nonzero.

The runner defaults to `g++` and honors `CXX` and `ASAN_OPTIONS`. LeakSanitizer is disabled by
default because it cannot run in traced environments. AddressSanitizer still checks invalid access and
double frees; the diagnostic's allocation accounting checks for C++ allocation leaks.
