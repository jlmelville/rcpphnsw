## Release Summary

This is a minor release to fix a SAN error reported by the CRAN Team, plus some
small feature additions.

## Test environments

* ubuntu 22.04 (on rhub) devel clang-ASAN
* Fedora 42 (on rhub) devel valgrind
* local ubuntu 26.04 R 4.5.2
* ubuntu 24.04 (on github actions), R 4.5.3, R 4.6.0, devel
* Windows Server 2022 (on github actions), R 4.5.3, R 4.6.0
* local Windows 11 build, R 4.6.0
* win-builder (devel)
* mac OS X Sequoia (on github actions) R 4.6.0
* local mac OS X Tahoe R 4.6.0
* mac-builder (devel)

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

## CRAN checks

There are no ERRORs, WARNINGs or NOTES as of today, but there is a badge
indicating an issue needs fixing. I was emailed by the CRAN team about a
failure in the examples under M1 SAN
<https://www.stats.ox.ac.uk/pub/bdr/M1-SAN/RcppHNSW. The error is:

```
../inst/include/hnswalg.h:1203:16: runtime error: store to misaligned address 
0x62a000078294 for type 'labeltype *' (aka 'unsigned long *'), which requires 8 
byte alignment
```

This submission is intended to fix this error.

## Downstream dependencies

We checked 9 reverse dependencies (7 from CRAN + 2 from Bioconductor), comparing R CMD check
results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
