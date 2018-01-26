#' GMM inference to estimate the LATE for a mismeasured binary treatment
#'
#' \code{LATEerror} implements the GMM inference to correct a measurement error
#' problem due to a mismeasured binary treatment in the LATE inference
#' developed in Yanagi (2017).
#' See the package vignette via \code{vignette("LATEerror")} for more details.
#'
#' @param Y vector of the outcome
#' @param D vector of the treatment
#' @param Z vector of the instrument
#' @param V vector of the exogenous variable
#' @param weight weighting matrix for the pilot GMM estimation
#' @param optimal logical whether implementing the optimal GMM estimation
#' @param equal logical whether the misclassification probabilities do not depend on the instrument
#' @param lower vector of lower bounds of the parameters
#' @param upper vector of upper bounds of the parameters
#' @param control list of the control parameters for \code{\link[DEoptim]{DEoptim}}
#'
#' @return list that contains the two lists: (i) \code{pilot}, (ii) \code{optimal}.\cr
#' \code{pilot} contains the following elements.
#' \item{estimate}{vector of parameter estimates}
#' \item{se}{vector of standard errors}
#' \item{ci}{matrix that contains 95 percent confidence intervals}
#' \item{optim}{list that contains the results of the optimization via \code{DEoptim}}
#' \item{weight}{weighting matrix for the pilot GMM estimation}
#' \code{optimal} contains the following elements.
#' \item{estimate}{vector of parameter estimates}
#' \item{se}{vector of standard errors}
#' \item{ci}{matrix that contains 95 percent confidence intervals}
#' \item{overidentification}{vector that contains a value of the test statistic and p-value for the over-identification test}
#' \item{weight}{weighting matrix for the optimal GMM estimation}
#' \item{optim}{list that contains the results of the optimization via \code{DEoptim}}
#'
#' @seealso \code{\link[DEoptim]{DEoptim}} for the elements in \code{control}
#'
#' @importFrom DEoptim DEoptim
#'
#' @export
#'
LATEerror <- function(Y, D, Z, V, weight = NULL, optimal = FALSE, equal = TRUE, lower = NULL, upper = NULL, control = NULL) {

  #-----------------------------------------------------------------------------
  # preparation
  #-----------------------------------------------------------------------------

  # sample size
  n <- length(Y)

  # stop with error messages
  if (n != length(D) || n != length(Z) || n != length(V)) {
    stop(message = "The lengths of Y, D, Z, and V must be the same.")
  }
  if (min(D) != 0 || max(D) != 1 || sum(D == 0) == 0 || sum(D == 1) == 0) {
    stop(message = "D must take two values of 0 and 1.")
  }
  if (min(Z) != 0 || max(Z) != 1 || sum(Z == 0) == 0 || sum(Z == 1) == 0) {
    stop(message = "Z must take two values of 0 and 1.")
  }

  # the number of the elements in the support of V
  K <- max(V)

  # the number of the parameters and the names of the parameters
  number_name <- number_name_par(equal = equal, K = K)
  number_par <- number_name$number_par
  name_par <- number_name$name_par

  # the number of the moments
  number_moment <- 3 + 4 * K

  # stop with error messages
  for (j in 1:K) {
    if (min(V) <= 0 || sum(V == j) == 0) {
      stop(message = "V must take K consecutive discrete values (1, 2, ..., K).")
    }
  }
  if (K <= 2 && equal == FALSE) {
    stop(message = "If 'equal = FALSE', V must take at least three values.")
  }
  if (!is.null(weight)) {
    if (!is.matrix(weight)) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    } else if (ncol(weight) != number_moment && nrow(weight) != number_moment) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    } else if (t(weight) != weight) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    }
  }
  if (!is.logical(optimal)) {
    stop(message = "'optimal' must be logical.")
  }

  # making data frame with the indicators
  data <- addindicators(Y = Y, D = D, Z = Z, V = V)

  # lower and upper bounds for optimization
  bound <- max(Y) - min(Y)
  if (is.null(lower)) {
    if (equal == FALSE) {
      lower <- c(-bound, 0.001, 0.001, -bound, -bound, rep(0.001, 4), rep(0.001, 2 * K))
    } else if (equal == TRUE){
      lower <- c(-bound, 0.001, 0.001, -bound, -bound, rep(0.001, 2), rep(0.001, 2 * K))
    }
  } else {
    if (!is.vector(lower)) {
      stop(message = "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")
    } else if (length(lower) != number_par) {
      stop(message = "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")
    } else if (equal == FALSE) {
      if (lower[6] + lower[7] >= 1 || lower[8] + lower[9] >= 1) {
        stop(message = "The sum of the misclassification probabilities must be smaller than one.")
      }
    } else if (equal == TRUE){
      if (lower[6] + lower[7] >= 1) {
        stop(message = "The sum of the misclassification probabilities must be smaller than one.")
      }
    }
  }
  if (is.null(upper)) {
    if (equal == FALSE) {
      upper <- c(bound, 0.999, 0.999, bound, bound, rep(0.499, 4), rep(0.999, 2 * K))
    } else if (equal == TRUE) {
      upper <- c(bound, 0.999, 0.999, bound, bound, rep(0.499, 2), rep(0.999, 2 * K))
    }
  } else {
    if (!is.vector(upper)) {
      stop(message = "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")
    } else if (length(upper) != number_par) {
      stop(message = "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")
    } else if (equal == FALSE) {
      if (upper[6] + upper[7] >= 1 || upper[8] + upper[9] >= 1) {
        stop(message = "The sum of the misclassification probabilities must be smaller than one.")
      }
    } else if (equal == TRUE) {
      if (upper[6] + upper[7] >= 1) {
        stop(message = "The sum of the misclassification probabilities must be smaller than one.")
      }
    }
  }

  # weight matrix for the pilot estimation
  if (is.null(weight)) {
    weight1 <- diag(number_moment)
  } else {
    weight1 <- weight
  }

  # control for DE algorithm
  if (!is.list(control) && !is.null(control)) {
    stop(message = "'control' must be correctly specified.")
  }
  if (is.null(control$NP)) {
    control$NP <- 10 * number_par
  }
  if (is.null(control$itermax)) {
    control$itermax <- 10000
  }
  if (is.null(control$steptol)) {
    control$steptol <- 500
  }
  if (is.null(control$reltol)) {
    control$reltol <- 1e-12
  }
  if (is.null(control$trace)) {
    control$trace <- 100
  }
  if (is.null(control$CR)) {
    control$CR <- 0.5
  }
  if (is.null(control$F)){
    control$F <- 0.8
  }
  if (is.null(control$c)){
    control$c <- 0.05
  }
  if (is.null(control$strategy)){
    control$strategy <- 2
  }
  if (is.null(control$parallelType)) {
    control$parallelType <- 1
  }

  #-----------------------------------------------------------------------------
  # pilot estimation
  #-----------------------------------------------------------------------------

  # optimization
  estDE1 <- DEoptim(fn = criterion, lower = lower, upper = upper, data = data, equal = equal, weight = weight1, control = control)

  # pilot estimates
  estimate1 <- estDE1$optim$bestmem

  # standard errors and confidence intervals
  inference1 <- inference(par = estimate1, data = data, weight = weight1, equal = equal)

  # pilot estimation result
  estimation1 <- list(estimate = inference1$estimate, se = inference1$se, ci = inference1$ci, optim = estDE1$optim, weight = weight1)

  if (optimal == FALSE) {
    result <- list(pilot = estimation1)
    return(result)
  }

  #-----------------------------------------------------------------------------
  # optimal estimation
  #-----------------------------------------------------------------------------

  # optimal weighting matrix
  var_moments_estimate1 <- var_moments(par = estimate1, data = data, equal = equal)
  weight2 <- solve(var_moments_estimate1)

  # optimization
  estDE2 <- DEoptim(fn = criterion, lower = lower, upper = upper, data = data, weight = weight2, equal = equal, control = control)

  # optimal estimates
  estimate2 <- estDE2$optim$bestmem

  # standard errors for optimal estimation
  inference2 <- inference(par = estimate2, data = data, weight = weight2, equal = equal)

  # over-identication test
  if (number_moment > number_par) {

    overidentification <- Jtest(par = estimate2, data = data, equal = equal)

  } else {

    overidentification <- NULL

  }

  # optimal estimation result
  estimation2 <- list(estimate = inference2$estimate, se = inference2$se,
                      ci = inference2$ci, overidentification = overidentification,
                      weight = weight2, optim = estDE2$optim)

  result <- list(pilot = estimation1, optimal = estimation2)

  return(result)

}
