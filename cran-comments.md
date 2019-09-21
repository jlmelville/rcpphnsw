## Release Summary

This is a minor release to upgrade the HNSW library to a new version, add some
new methods and to fix a user-reported bug with the search function. Also, a URL
and BugReports entry has been added to the DESCRIPTION.

## Test environments

* ubuntu 14.04 (on travis-ci), R 3.4.4, R 3.6.0, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 30 (on rhub), R-devel
* debian (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.1
* local Windows 10 build, R 3.6.1
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.6.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was a message about possibly mis-spelled words in DESCRIPTION:

  HNSW (2:28)
 
This is spelled correctly.

With r-hub checking on Windows only there was a message:

"N  checking for non-standard things in the check directory
   Found the following files/directories:
     'examples_x64' 'tests_i386' 'tests_x64'
     'RcppHNSW-Ex_i386.Rout' 'RcppHNSW-Ex_x64.Rout' 'examples_i386'"

This would seem to be something to do with r-hub rather than a real problem.

## Downstream dependencies

None.
