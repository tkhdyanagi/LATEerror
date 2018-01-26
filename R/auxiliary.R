#' adding indicators for exogenous variable
#'
#' \code{addindicators} is an auxiliary function for preparation to implement the GMM inference on your own.
#' It returns a data.frame that contains vectors of the outcome, the treatment, the instrument, and the exogenous variable in addition to vectors of indicators for the exogenous variable.
#' The defined indicator is such as I00 = 1(Z = 0, V = 0).
#' See the package vignette via `vignette("LATEerror")` for more details.
#'
#' @param Y vector of the outcome
#' @param D vector of the treatment
#' @param Z vector of the instrument
#' @param V vector of the exogenous variable
#'
#' @return dataframe with indicators for the exogenous variable
#'
#' @export
#'
addindicators <- function(Y, D, Z, V) {

  # make dataframe
  data <- data.frame(Y = Y, D = D, Z = Z, V = V)

  # the number of the elements that V takes
  K <- max(data$V)

  # add indicators
  for (j in 1:K) {
    varname0 <- sprintf("I0%d", j)
    data[[varname0]] <- (1 - data$Z) * ifelse(data$V == j, 1, 0)
    varname1 <- sprintf("I1%d", j)
    data[[varname1]] <- data$Z * ifelse(data$V == j, 1, 0)
  }

  return(data)

}



#' returning the names and the numbers of the parameters
#'
#' @param equal logical whether misclassification probabilities do not depend on Z
#' @param K the number of the elements that V takes
#'
number_name_par <- function(equal, K) {

  if (equal == FALSE) {

    # the number of the parameters
    number_par <- 2 * K + 9

    # the name of the parameters
    name_par <- c("LATE", "first-stage", "E(Z)", "tau0*", "tau1*", "m00", "m01", "m10", "m11")
    for (j1 in 1:K) {
      name_par <- c(name_par, sprintf("p0%d*", j1))
    }
    for (j2 in 1:K) {
      name_par <- c(name_par, sprintf("p1%d*", j2))
    }

  } else if (equal == TRUE) {

    # the number of the parameters
    number_par <- 2 * K + 7

    # the name of the parameters
    name_par <- c("LATE", "first-stage", "E(Z)", "tau0*", "tau1*", "m0", "m1")
    for (j1 in 1:K){
      name_par <- c(name_par, sprintf("p0%d*", j1))
    }
    for (j2 in 1:K){
      name_par <- c(name_par, sprintf("p1%d*", j2))
    }

  }

  # result
  result <- list(number_par = number_par, name_par = name_par)

  return(result)

}
