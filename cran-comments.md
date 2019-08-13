Date: 08/12/2019
## Test environments
+ local OS X install (macOS Mojave 10.14.5), R version 3.6.1
+ NEED TO TEST travis
+ win-builder (devel - 2 NOTEs;  release - 2 NOTEs)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs in local environments

devel/release win-builder identical 2 NOTEs shown below

1.  Maintainer: 'Justin Williams <williazo@ucla.edu>'
    New Submission
    Possibly mis-spelled words in DESCRIPTION:
    Raphson (9:34)
    analytical (9:36)
2.  Non-standard file/directory found at top level:
    'README_files'
In NOTE 1. the identified words are correctly spelled.
In NOTE 2. I use the README_files directory to store images used in README.MD as
the illustrations help illustrate the concepts.

## Downstream dependencies
There are currently no downstream dependencies for this package
