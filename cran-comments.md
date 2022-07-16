## Release Summary

This is a new minor release and also fixes a 'LazyData' NOTE reported on the
CRAN checks page.

## Test environments

* ubuntu 20.04 (on github actions), R 4.1.3, R 4.2.1, devel
* local ubuntu 22.04 R 4.2.1
* Windows Server 2022 (on github actions), R 4.1.3
* local Windows 11 build, R 4.2.1
* mac OS X Big Sur (on github actions) R 4.2.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

❯ checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      libs   5.1Mb

This is expected due to the use of C++ templates in hnswlib.

There was a message about possibly mis-spelled words in DESCRIPTION:

  HNSW (2:28)
 
This is spelled correctly.

## CRAN checks

There are no ERRORs or WARNINGs.

There is a NOTE for all flavors about LazyData. This release fixes that NOTE.

There are four flavors with NOTEs about installed package size 
(r-release-macos-arm64, r-release-macos-x86_64, r-oldrel-macos-arm64, 
r-oldrel-macos-x86_64). This is expected and won't be fixed.

## Downstream dependencies

We checked 3 reverse dependencies (1 from CRAN + 2 from Bioconductor), comparing 
R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
