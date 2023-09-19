## Release Summary

This is a patch release to fix various CRAN check errors.

## Test environments

* ubuntu 22.04 (on github actions), R 4.2.3, R 4.3.1, devel
* local ubuntu 23.04 R 4.2.2
* Debian Linux, R-devel, GCC ASAN/UBSAN (via rhub)
* Debian Linux, R-release, GCC (via rhub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (via rhub)
* Fedora Linux, R-devel, clang, gfortran (via rhub)
* Windows Server 2022 (on github actions), R 4.2.3, R 4.3.1
* Windows Server 2022, R-devel, 64 bit (via rhub)
* local Windows 11 build, R 4.3.1
* win-builder (devel)
* mac OS X Monterey (on github actions) R 4.3.1

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

N  checking installed package size ...
     installed size is  6.6Mb
     sub-directories of 1Mb or more:
       libs   6.3Mb

This is expected due to the use of C++ templates in hnswlib.
 
This is spelled correctly.

## CRAN checks

There are no ERRORs or WARNINGs.

There is a NOTE:

Check: C++ specification
Result: NOTE
     Specified C++11: please drop specification unless essential

This submission fixes this.

There is a NOTE:

Check: Rd metadata
Result: NOTE
    Invalid package aliases in Rd file ‘RcppHnsw-package.Rd’:
     ‘RcppHnsw-package’

This submissions fixes this.

There are four flavors with NOTEs about installed package size (r-release-macos-arm64,
r-release-macos-x86_64, r-oldrel-macos-arm64, r-oldrel-macos-x86_64). This is expected and won't be
fixed.

## Downstream dependencies

We checked 2 reverse dependencies (0 from CRAN + 2 from Bioconductor), comparing R CMD check
results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
