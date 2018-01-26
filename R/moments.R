#' Computing the GMM criterion.
#'
#' \code{criterion} is an auxiliary function to implement the GMM inference on your own.
#' It returns a value of the GMM criterion.
#'
#' @param par vector of parameter values
#' @param data dataframe that contains Y, D, Z, and V with indicators
#' @param weight weighting matrix for GMM estimation
#' @param equal logical whether the misclassification probabilities do not depend on instrument
#'
#' @export
#'
criterion <- function(par, data, weight, equal) {

  # the means of the moments
  S <- mean_moments(par = par, data = data, equal = equal)

  # criterion
  return(t(S) %*% weight %*% S)

}



#' Computing the means of the moments for the GMM estimation.
#'
#' @param par vector of parameter values
#' @param data dataframe containing the outcome, treatment, instrument, exogenous variable
#' @param equal logical whether the misclassification probabilities do not depend on instrument
#'
mean_moments <- function(par, data, equal) {

  # the means of the moments
  result <- colMeans(moments(par = par, data = data, equal = equal))

  return(result)

}



#' Computing the values of the moments for the GMM estimation.
#'
#' @param par vector of parameter values
#' @param data dataframe containing the outcome, treatment, instrument, exogenous variable
#' @param equal logical whether the misclassification probabilities do not depend on instrument
#'
moments <- function(par, data, equal) {

  # the number of the elements V takes
  K <- max(data$V)

  if (equal == FALSE) {

    # parameters
    beta <- par[1]
    delta_p <- par[2]
    r <- par[3]
    tau0 <- par[4]
    tau1 <- par[5]
    m00 <- par[6]
    m01 <- par[7]
    m10 <- par[8]
    m11 <- par[9]
    j <- 0
    for (j in 1:K) {
      assign(sprintf("p0%d", j), par[9 + j])
      assign(sprintf("p1%d", j), par[9 + K + j])
    }

    # functions for moments
    f1 <- beta - (data$Y * data$Z / r - data$Y * (1 - data$Z) / (1 - r)) / delta_p
    f2 <- delta_p - ( (data$D * data$Z / r - m01) / (1 - m01 - m11) - (data$D * (1 - data$Z) / (1 - r) - m00) / (1 - m00 - m10) )
    f3 <- r - data$Z
    j <- 0
    for (j in 1:K) {
      I0v <- data[[sprintf("I0%d", j)]]
      I1v <- data[[sprintf("I1%d", j)]]
      p0v <- get(sprintf("p0%d", j))
      p1v <- get(sprintf("p1%d", j))
      assign(sprintf("f%d", 3 + j), (m00 + (1 - m00 - m10) * p0v - data$D) * I0v)
      assign(sprintf("f%d", 3 + K + j), (m01 + (1 - m01 - m11) * p1v - data$D) * I1v)
      assign(sprintf("f%d", 3 + 2 * K + j), (tau0 + (data$Y * data$D - (1 - m10) * p0v * tau0) / (m00 + (1 - m00 - m10) * p0v) - (data$Y * (1 - data$D) + (1 - m00) * (1 - p0v) * tau0) / (1 - (m00 + (1 - m00 - m10) * p0v))) * I0v)
      assign(sprintf("f%d", 3 + 3 * K + j), (tau1 + (data$Y * data$D - (1 - m11) * p1v * tau1) / (m01 + (1 - m01 - m11) * p1v) - (data$Y * (1 - data$D) + (1 - m01) * (1 - p1v) * tau1) / (1 - (m01 + (1 - m01 - m11) * p1v))) * I1v)
    }

    # return functions for moments
    f <- NULL
    j <- 0
    for (j in 1:(3 + 4 * K)) {
      f <- cbind(f, get(sprintf("f%d", j)))
    }

    return(f)

  } else if (equal == TRUE) {

    # parameters
    beta <- par[1]
    delta_p <- par[2]
    r <- par[3]
    tau0 <- par[4]
    tau1 <- par[5]
    m0 <- par[6]
    m1 <- par[7]
    j <- 0
    for (j in 1:K) {
      assign(sprintf("p0%d", j), par[7 + j])
      assign(sprintf("p1%d", j), par[7 + K + j])
    }
    # functions for moments
    f1 <- beta - (data$Y * data$Z / r - data$Y * (1 - data$Z) / (1 - r)) / delta_p
    f2 <- delta_p - (data$D * data$Z / r - data$D * (1 - data$Z) / (1 - r)) / (1 - m0 - m1)
    f3 <- r - data$Z
    j <- 0
    for (j in 1:K) {
      I0v <- data[[sprintf("I0%d", j)]]
      I1v <- data[[sprintf("I1%d", j)]]
      p0v <- get(sprintf("p0%d", j))
      p1v <- get(sprintf("p1%d", j))
      assign(sprintf("f%d", 3 + j), (m0 + (1 - m0 - m1) * p0v - data$D) * I0v)
      assign(sprintf("f%d", 3 + K + j), (m0 + (1 - m0 - m1) * p1v - data$D) * I1v)
      assign(sprintf("f%d", 3 + 2 * K + j), (tau0 + (data$Y * data$D - (1 - m1) * p0v * tau0) / (m0 + (1 - m0 - m1) * p0v) - (data$Y * (1 - data$D) + (1 - m0) * (1 - p0v) * tau0) / (1 - (m0 + (1 - m0 - m1) * p0v))) * I0v)
      assign(sprintf("f%d", 3 + 3 * K + j), (tau1 + (data$Y * data$D - (1 - m1) * p1v * tau1) / (m0 + (1 - m0 - m1) * p1v) - (data$Y * (1 - data$D) + (1 - m0) * (1 - p1v) * tau1) / (1 - (m0 + (1 - m0 - m1) * p1v))) * I1v)
    }

    # return functions for moments
    f <- NULL
    j < - 0
    for (j in 1:(3 + 4 * K)) {
      f <- cbind(f, get(sprintf("f%d", j)))
    }

    return(f)
  }
}
