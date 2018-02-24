<!-- README.md is generated from README.Rmd. Please edit that file -->
LATEerror
=========

This R package implements the generalized method of moments estimation for the local average treatment effect when the binary treatment variable may be mismeasured. The estimation procedure is developed by Yanagi (2017).

Intsallation
============

Run the following codes to install the package.

``` r
install.packages("devtools") # if needed
devtools::install_github("tkhdyanagi/LATEerror", build_vignettes = TRUE)
```

You can see the usage of the package by the following codes.

``` r
library("LATEerror")
vignette("LATEerror")
```
