#' computing standard errors and confidence intervals
#'
#' \code{inference} is an auxiliary function to implement the GMM inference on your own.
#' It returns estimates, standard errors, and confidence intervals
#' for the parameters in the GMM estimation.
#' See the package vignette via `vignette("LATEerror")` for more details.
#'
#' @param par vector of parameter values
#' @param data dataframe that contains Y, D, Z, and V with indicators
#' @param weight weighting matrix for pilot GMM estimation
#' @param equal logical whether misclassification probabilities do not depend on Z
#'
#' @return list that contains the following elements.
#' \item{estimate}{parameter estimates}
#' \item{se}{standard errors}
#' \item{ci}{95 percent confidence intervals}
#'
#' @importFrom MASS ginv
#' @importFrom numDeriv jacobian
#' @importFrom stats qnorm
#'
#' @export
#'
inference <- function(par, data, weight, equal) {

  # sample size
  n <- length(data$Y)

  # the number of the elements in the support of V
  K <- max(data$V)

  # the number of the parameters and the names of the parameters
  temp <- number_name_par(equal = equal, K = K)
  number_par <- temp$number_par
  name_par <- temp$name_par

  # jacobian of the moments
  G <- jacobian(mean_moments, par, data = data, equal = equal)

  # variance-covariance matrix of the moments
  var_moments_estimate <- var_moments(par = par, data = data, equal = equal)

  # computing inverse matrix
  if (qr(t(G) %*% weight %*% G)$rank < number_par) {

    covariance <- ginv(t(G) %*% weight %*% G) %*%
      t(G) %*% weight %*% var_moments_estimate %*% weight %*% G %*%
      ginv(t(G) %*% weight %*% G) / n

    warning("The asymptotic variance-covariance matrix estimation produces a singular matrix.
            The generalized inverse is used.")

  } else {

    covariance <- solve(t(G) %*% weight %*% G) %*%
      t(G) %*% weight %*% var_moments_estimate %*% weight %*% G %*%
      solve(t(G) %*% weight %*% G) / n

  }

  # standard errors
  se <- sqrt(diag(covariance))

  # 95% confidence intervals
  ci <- rbind(par + qnorm(0.025) * se, par + qnorm(0.975) * se)

  # names
  names(par) <- name_par
  names(se) <- name_par
  colnames(ci) <- name_par
  rownames(ci) <- c("95% CI lower", "95% CI upper")

  # results
  result <- list(estimate = par, se = se, ci = ci)

  return(result)

}



#' conducting over-identification test (J test)
#'
#' @param par GMM estimates
#' @param data dataframe containing the outcome, treatment, instrument, exogenous variable with indicators
#' @param equal logical whether the misclassification probabilities do not depend on instrument
#'
#' @return vector that contains a test statistic value and p-value
#'
#' @importFrom stats pchisq
#'
Jtest <- function(par, data, equal) {

  # sample size
  n <- length(data$Y)

  # the number of the elements in the support of V
  K <- max(data$V)

  # the number of the moments
  number_moment <- 3 + 4 * K

  # the number of the parameters and the names of the parameters
  temp <- number_name_par(equal = equal, K = K)
  number_par <- temp$number_par
  name_par <- temp$name_par

  # the means of the moments
  S <- mean_moments(par = par, data = data, equal = equal)

  # over-identification test statistic
  var_moments_estimate <- var_moments(par = par, data = data, equal = equal)
  over_test <- n * S %*% solve(var_moments_estimate) %*% S

  # p value
  df <- number_moment - number_par
  over_test_p <- 1 - pchisq(q = over_test, df = df)

  # results
  overidentification <- c(over_test, over_test_p)
  names(overidentification) <- c("test statistic", "p-value")

  return(overidentification)

}



#' Estimating variance-covariance matrix of the moments
#'
#' @param par GMM estimates
#' @param data dataframe containing the outcome, treatment, instrument, exogenous variable
#' @param equal logical whether the misclassification probabilities do not depend on instrument
#'
var_moments <- function(par, data, equal) {

  # the number of the elements V takes
  K <- max(data$V)

  # the number of the moments
  number_moment <- 3 + 4 * K

  # the values of the moments
  f <- moments(par = par, data = data, equal = equal)

  # variance-covariance matrix estimation
  hat_gamma <- NULL
  for (j in 1:number_moment) {
    hat_gamma <- rbind(hat_gamma, colMeans(f * f[, j]))
  }

  return(hat_gamma)

}
