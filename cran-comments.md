## Resubmission

This is a resubmission. In this version I have:

* Added a reference on the methods used in the DESCRIPTION.

* Deleted the function gamma_coef which was based on code written by another author. This function is already included  in the R package vcdExtra and I import it now from vcdExtra.

* Added  importFrom stats  "coef", "complete.cases",   "na.exclude", "optimize","pchisq", "pnorm", "qnorm", "rnorm", "runif", "uniroot","xtabs"  to the NAMESPACE file.



## Tests environment 

* local Windows 7 install, R 3.5.2 (64 bit)
* win-builder R Under development (unstable) (2019-01-31 r76038)
* ubuntu 12.04 (on travis-ci), R 3.1.2
* Mac OS X 10.13.3 (on travis-ci)


## R CMD check results
There were no errors or warnings.

There was 1  NOTE:

* checking CRAN incoming feasibility ... NOTE
* Maintainer: 'Marianne Mueller <mllm@zhaw.ch>'

* New submission

This is my first submission 

