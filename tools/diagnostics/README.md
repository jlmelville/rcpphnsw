# hnswlib insertion ThreadSanitizer diagnostic

`hnswlib-rng-tsan.cpp` exercises concurrent unique-label insertion into one
vendored hnswlib index after serially seeding the first item, matching the
package construction shape. Its focused purpose is to exercise the two
localized concurrency repairs; it is not a claim that all hnswlib features or
an R process are otherwise clean under ThreadSanitizer.

Run the bounded default workload from any working directory:

```sh
tools/diagnostics/run-hnswlib-rng-tsan.sh
```

The optional positional arguments are rounds, worker threads, points per
thread, `M`, and the random seed. For example:

```sh
tools/diagnostics/run-hnswlib-rng-tsan.sh 8 12 256 16 101
```

The runner uses `clang++` by default, honors `CXX` and `TSAN_OPTIONS`, builds in
a temporary directory, and removes its binary on exit. Its default
`TSAN_OPTIONS` disables the separate lock-order detector so this witness stays
focused on data races; set `TSAN_OPTIONS` explicitly for broader diagnostics.
A clean run prints one completion summary and confirms every round changed the
index's maximum level after the serial seed. A race report exits nonzero and
should contain worker `addPoint()` stacks. The accepted pre-patch witnesses
named `getRandomLevel()` and the entry-point snapshot in `hnswalg.h`; another
independent race on this supported path triggers package-level serialization,
not another open-ended vendor repair.
