Date: 12/17/2019
### Re-submission
This is a resubmission of tcensReg 0.1.3 to tcensReg 0.1.4 In this version I have:
 + removed inst/script/* files from the Rbuild. These were codes run in the simulation for the manuscript using the package, but are unnecessary for the package. I will retain access to them on the GitHub page but ignore them for the .tar.gz build.


## Test environments
+ local OS X install (macOS Catalina 10.15.1), R version 3.6.1
+ win-builder (devel & release)

### local OS X results
R CMD check --as-cran
There were no ERRORs or WARNINGs and 1 NOTE in local environments

1. Maintainer: ‘Justin Williams <williazo@ucla.edu>’
New submission

*Response:* This is a new package submission.

### win-builder devel & release results
There were no ERRORs or WARNINGs and 1 NOTE

*Response:* Same note as above

## Downstream dependencies
There are currently no downstream dependencies for this package
