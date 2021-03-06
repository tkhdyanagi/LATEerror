% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moments.R
\name{mean_moments}
\alias{mean_moments}
\title{Computing the means of the moments for the GMM estimation.}
\usage{
mean_moments(par, data, equal)
}
\arguments{
\item{par}{vector of parameter values}

\item{data}{dataframe containing the outcome, treatment, instrument, exogenous variable}

\item{equal}{logical whether the misclassification probabilities do not depend on instrument}
}
\description{
Computing the means of the moments for the GMM estimation.
}
