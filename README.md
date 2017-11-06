<!-- README.md is generated from README.Rmd. Please edit that file -->
LATEerror
=========

**LATEerror** package implements the generalized method of moments estimation for the local average treatment effect when the treatment variable may be mismeasured. The estimation procedure is developed by Yanagi (2017).

Intsallation
============

Run the following code to install **LATEerror**.

``` r
install.packages("devtools")
devtools::install_github("tkhdyanagi/LATEerror")
library(LATEerror)
```

Usage
=====

``` r
LATEerror(Y, D, Z, V, weight = NULL, optimal = TRUE, equal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL)
```

See also
========

1.  See package vignettes "How to use `LATEerror`" for more details.
"# LATEerror" 
