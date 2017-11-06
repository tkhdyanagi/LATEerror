#' Implementing the GMM inference to estimate the LATE for misclassified treatment.
#'
#' @param Y A vector of outcome variable.
#' @param D A vector of treatment variable.
#' @param Z A vector of instrument.
#' @param V A vector of exogenous variable.
#' @param weight A matrix for the weighting matrix for pilot GMM estimation.
#' @param optimal Logical. Whether implementing the optimal GMM estimation?
#' @param equal Logical. Whether imposing restrictions that the misclassification probabilities do not depend on Z?
#' @param lower A vector for lower bounds of parameters.
#' @param upper A vector for upper bounds of parameters.
#' @param controlDE A list for control parameters for \code{\link[DEoptim]{DEoptim}}.
#' @param controlBFGS A list for control parameters for \code{\link[stats]{optim}} with L-BFGS-B method.
#'
#' @return The output is a list that contains the following two list.
#'   (1) A list for pilot estimation results by the pilot GMM estimation based on weight.
#'   (2) A list for optimal estimation results by the optimal GMM estimation based on the optimal weighting matrix.
#'
#' @seealso \code{\link[DEoptim]{DEoptim}} for controlDE and \code{\link[stats]{optim}} for controlBFGS.
#'
#' @importFrom DEoptim DEoptim
#' @importFrom MASS ginv
#' @importFrom numDeriv jacobian
#' @importFrom stats optim pchisq qnorm
#'
#' @export
LATEerror <- function(Y, D, Z, V, weight = NULL, optimal = TRUE, equal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL) {

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
  if (equal == FALSE) {
    number_par <- 2 * K + 9
    name_par <- c("LATE", "first-stage", "E(Z)", "tau0*", "tau1*", "m00", "m01", "m10", "m11")
    j <- 0
    for (j in 1:K) {
      name_par <- c(name_par, sprintf("p0%d*", j))
    }
    j <- 0
    for (j in 1:K) {
      name_par <- c(name_par, sprintf("p1%d*", j))
    }
  } else if (equal == TRUE) {
    number_par <- 2 * K + 7
    name_par <- c("LATE", "first-stage", "E(Z)", "tau0*", "tau1*", "m0", "m1")
    j <- 0
    for (j in 1:K){
      name_par <- c(name_par, sprintf("p0%d*", j))
    }
    j <- 0
    for (j in 1:K){
      name_par <- c(name_par, sprintf("p1%d*", j))
    }
  }

  # the number of the moments
  number_moment <- 3 + 4 * K

  # stop with error messages
  j <- 0
  for (j in 1:K) {
    if (min(V) <= 0 || sum(V == j) == 0) {
      stop(message = "V must take K consecutive discrete values (1, 2, ..., K).")
    }
  }

  if (K <= 2 && equal == FALSE) {
    stop(message = "If 'equal = FALSE', V must take at least three values.")
  }

  if (is.null(weight) == FALSE) {
    if (is.matrix(weight) == FALSE) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    } else if (ncol(weight) != number_moment && nrow(weight) != number_moment) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    } else if (t(weight) != weight) {
      stop(message = "'weight' must be NULL or a symmetric matrix.")
    }
  }

  if (is.logical(optimal) == FALSE) {
    stop(message = "'optimal' must be logical.")
  }

  # making data frame with indicators
  data <- data.frame(Y = Y, D = D, Z = Z, V = V)
  j <- 0
  for (j in 1:K) {
    varname0 <- sprintf("I0%d", j)
    data[[varname0]] <- (1 - Z) * ifelse(V == j, 1, 0)
    varname1 <- sprintf("I1%d", j)
    data[[varname1]] <- Z * ifelse(V == j, 1, 0)
  }

  # lower and upper bounds for optimization
  bound <- max(Y) - min(Y)

  if (is.null(lower)) {
    if (equal == FALSE) {
      lower <- c(-bound, 0.001, 0.001, -bound, -bound, rep(0.001, 4), rep(0.001, 2 * K))
    } else if ( equal == TRUE ){
      lower <- c(-bound, 0.001, 0.001, -bound, -bound, rep(0.001, 2), rep(0.001, 2 * K))
    }
  } else {
    if (is.vector(lower) == FALSE) {
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
    if (is.vector(upper) == FALSE) {
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

  # weight matrix for pilot estimation
  if (is.null(weight)) {
    weight1 <- diag(number_moment)
  } else {
    weight1 <- weight
  }

  # control for DE algorithm
  if (is.list(controlDE) == FALSE && is.null(controlDE) == FALSE) {
    stop(message = "'controlDE' must be correctly specified.")
  }

  if (is.null(controlDE$NP)) {
    controlDE$NP <- 10 * number_par
  }

  if (is.null(controlDE$itermax)) {
    controlDE$itermax <- 10000
  }

  if (is.null(controlDE$steptol)) {
    controlDE$steptol <- 500
  }

  if (is.null(controlDE$reltol)) {
    controlDE$reltol <- 1e-12
  }

  if (is.null(controlDE$trace)) {
    controlDE$trace <- 100
  }

  if (is.null(controlDE$CR)) {
    controlDE$CR <- 0.5
  }

  if (is.null(controlDE$F)){
    controlDE$F <- 0.8
  }

  if (is.null(controlDE$c)){
    controlDE$c <- 0.05
  }

  if (is.null(controlDE$strategy)){
    controlDE$strategy <- 2
  }

  if (is.null(controlDE$parallelType)) {
    controlDE$parallelType <- 1
  }

  # control for optim
  if (is.list(controlBFGS) == FALSE && is.null(controlBFGS) == FALSE) {
    stop(message = "'controlBFGS' must be correctly specified.")
  }

  if (is.null(controlBFGS$maxit)) {
    controlBFGS$maxit <- 1000
  }

  # pilot estimation------------------------------------------------------------
  # optimization based on DE algorithm for pilot estimation
  estDE1 <- DEoptim(fn = criterion, lower = lower, upper = upper, data = data, equal = equal, weight = weight1, control = controlDE)

  # optimization check based on L-BFGS-B method for pilot estimation
  estBFGS1 <- optim(par = estDE1$optim$bestmem, fn = criterion, lower = lower, upper = upper, method = "L-BFGS-B", control = controlBFGS, data = data, weight = weight1, equal = equal)

  # pilot estimation results
  estimate1 <- estBFGS1$par

  # computing standard errors for pilot estimation
  G1 <- jacobian(mean_moments, estimate1, data = data, equal = equal)
  var_moments_estimate1 <- var_moments(par = estimate1, data = data, equal = equal)
  if (qr(t(G1) %*% weight1 %*% G1)$rank < number_par) {

    covariance1 <- ginv(t(G1) %*% weight1 %*% G1) %*%
      t(G1) %*% weight1 %*% var_moments_estimate1 %*% weight1 %*% G1 %*%
      ginv(t(G1) %*% weight1 %*% G1) / n

    warning("A matrix for asymptotic covariance matrix estimation is singular.
            The generalized inverse is used.")

  } else {

    covariance1 <- solve(t(G1) %*% weight1 %*% G1) %*%
      t(G1) %*% weight1 %*% var_moments_estimate1 %*% weight1 %*% G1 %*%
      solve(t(G1) %*% weight1 %*% G1) / n

  }
  se1 <- sqrt(diag(covariance1))

  # computing confidence intervals for pilot estimation
  ci1 <- rbind(estimate1 + qnorm(0.025) * se1, estimate1 + qnorm(0.975) * se1)

  # reporting pilot estimation results
  names(estimate1) <- name_par
  names(se1) <- name_par
  colnames(ci1) <- name_par
  rownames(ci1) <- c("95%CI lower", "95%CI upper")

  estimation1 <- list(estimate = estimate1, se = se1, ci = ci1, DEresult = estDE1$optim, BFGSresult = estBFGS1)

  if (optimal == FALSE) {
    result <- list(pilot = estimation1)
    return(result)
  }

  # optimal estimation----------------------------------------------------------
  # weight matrix for optimal estimation
  weight2 <- solve(var_moments_estimate1)

  # optimization based on DE method for optimal estimation
  estDE2 <- DEoptim(fn = criterion, lower = lower, upper = upper, data = data, weight = weight2, equal = equal, control = controlDE)

  # optimization check by L-BFGS-B method for optimal estimation
  estBFGS2 <- optim(par = estDE2$optim$bestmem, fn = criterion, lower = lower, upper = upper, method = "L-BFGS-B", control = controlBFGS, data = data, weight = weight2, equal = equal)

  # optimal estimation results
  estimate2 <- estBFGS2$par

  # computing standard errors for optimal estimation
  G2 <- jacobian(mean_moments, estimate2, data = data, equal = equal)
  var_moments_estimate2 <- var_moments(par = estimate2, data = data, equal = equal)
  if (qr(t(G2) %*% weight2 %*% G2)$rank < number_par) {

    covariance2 <- ginv(t(G2) %*% weight2 %*% G2) %*%
      t(G2) %*% weight2 %*% var_moments_estimate2 %*% weight2 %*% G2 %*%
      ginv(t(G2) %*% weight2 %*% G2) / n

    warning("A matrix for asymptotic covariance matrix estimation is singular.
            The generalized inverse is used.")

  } else {

    covariance2 <- solve(t(G2) %*% weight2 %*% G2) %*%
      t(G2) %*% weight2 %*% var_moments_estimate2 %*% weight2 %*% G2 %*%
      solve(t(G2) %*% weight2 %*% G2) / n

  }

  se2 <- sqrt(diag(covariance2))

  # computing confidence intervals for optimal estimation
  ci2 <- rbind(estimate2 + qnorm(0.025) * se2, estimate2 + qnorm(0.975) * se2)

  # reporting optimal estimation results
  names(estimate2) <- name_par
  names(se2) <- name_par
  colnames(ci2) <- name_par
  rownames(ci2) <- c("95%CI lower", "95%CI upper")

  # conducting over-identication test
  if (number_moment > number_par) {

    S2 <- mean_moments(par = estimate2, data = data, equal = equal)
    over_test <- n * S2 %*% solve(var_moments_estimate2) %*% S2
    over_test_p <- 1 - pchisq(q = over_test, df = (number_moment - number_par))
    overidentification <- c(over_test, over_test_p)
    names(overidentification) <- c("test statistic", "p-value")

    estimation2 <- list(estimate = estimate2, se = se2, ci = ci2, overidentification = overidentification, DEresult = estDE2$optim, BFGSresult = estBFGS2)

  } else {

    estimation2 <- list(estimate = estimate2, se = se2, ci = ci2, DEresult = estDE2$optim, BFGSresult = estBFGS2)

  }

  return(list(pilot = estimation1, optimal = estimation2))

}
