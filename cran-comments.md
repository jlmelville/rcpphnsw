## Release Summary

This is a new minor release.

## Test environments

* ubuntu 16.04 (on travis-ci), R 3.6.3, R 4.0.0, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 32 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.6.3, R 4.0.2
* local Windows 10 build, R 4.0.2
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 4.0.2
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was a message about possibly mis-spelled words in DESCRIPTION:

  HNSW (2:28)
 
This is spelled correctly.

With r-hub checking on Debian, there was the following message:

> * checking installed package size ... NOTE
  installed size is  7.3Mb
  sub-directories of 1Mb or more:
    libs   7.1Mb

This is expected due to the use of C++ templates in HNSW.

## CRAN checks

There are no ERRORs or WARNINGs.

There are three flavors with NOTEs about installed package size 
(r-devel-linux-x86_64-fedora-clang, r-release-macos-x86_64, 
r-oldrel-macos-x86_64). This is expected and won't be fixed.

## Downstream dependencies

There are 2 reverse dependencies (0 from CRAN + 2 from BioConductor):

* BiocNeighbors has no ERRORs or WARNINGs.
* Destiny fails to complete the R CMD check due to problems building vignettes.
This is unrelated to either the current CRAN version of RcppHNSW or this new
release.
