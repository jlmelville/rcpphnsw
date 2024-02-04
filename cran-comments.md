## Release Summary

This is a patch release for a new version of the underlying hnswlib library.

## Test environments

* ubuntu 22.04 (on github actions), R 4.2.3, R 4.3.2, devel
* local ubuntu 23.10 R 4.3.1
* Debian Linux, R-devel, GCC ASAN/UBSAN (via rhub)
* Debian Linux, R-release, GCC (via rhub)
* Debian Linux, R-release, GCC valgrind (via rhub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (via rhub)
* Fedora Linux, R-devel, clang, gfortran (via rhub)
* Windows Server 2022 (on github actions), R 4.2.3, R 4.3.2
* Windows Server 2022, R-devel, 64 bit (via rhub)
* local Windows 11 build, R 4.3.2
* win-builder (devel)
* local mac OS X Sonoma R 4.3.2
* mac OS X Monterey (on github actions) R 4.3.2

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

* checking installed package size ... NOTE
  installed size is  6.7Mb
  sub-directories of 1Mb or more:
    libs   6.4Mb

This is expected due to the use of C++ templates in hnswlib.
 
## CRAN checks

There are no ERRORs or WARNINGs.

There are three flavors with NOTEs about installed package size (r-release-macos-arm64,
r-release-macos-x86_64, r-oldrel-macos-arm64). This is expected and won't be fixed.

## Downstream dependencies

We checked 3 reverse dependencies (1 from CRAN + 2 from Bioconductor), comparing R CMD check
results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
