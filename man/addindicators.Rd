% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{addindicators}
\alias{addindicators}
\title{adding indicators for exogenous variable}
\usage{
addindicators(Y, D, Z, V)
}
\arguments{
\item{Y}{vector of the outcome}

\item{D}{vector of the treatment}

\item{Z}{vector of the instrument}

\item{V}{vector of the exogenous variable}
}
\value{
dataframe with indicators for the exogenous variable
}
\description{
\code{addindicators} is an auxiliary function for preparation to implement the GMM inference on your own.
It returns a data.frame that contains vectors of the outcome, the treatment, the instrument, and the exogenous variable in addition to vectors of indicators for the exogenous variable.
The defined indicator is such as I00 = 1(Z = 0, V = 0).
See the package vignette via `vignette("LATEerror")` for more details.
}
