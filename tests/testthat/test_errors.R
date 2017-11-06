
context("Testing errors")

# simulated data
n <- 1000
set.seed(102)
mydata <- data.frame(U1 = rnorm(n, 0, 1), U2 = rnorm(n, 0, 1))
mydata <- mydata %>%
  dplyr::mutate(Z = rbinom(n, size = 1, prob = 0.5)) %>%
  dplyr::mutate(V = rbinom(n, size = 1, prob = 0.5) + 1) %>%
  dplyr::mutate(Ds = ifelse(( -1.25 + Z + 0.5 * V - U1 >= 0), 1, 0)) %>%
  dplyr::mutate(D = Ds * rbinom(n, size = 1, prob = 0.8) + (1 - Ds) * rbinom(n, size = 1, prob = 0.2)) %>%
  dplyr::mutate(Y = 1 + Ds + U2)

# errors for data
expect_error(LATEerror(Y = NA, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "The lengths of Y, D, Z, and V must be the same.")

expect_error(LATEerror(Y = mydata$Y, D = NULL, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "The lengths of Y, D, Z, and V must be the same.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = 1, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "The lengths of Y, D, Z, and V must be the same.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = Inf, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "The lengths of Y, D, Z, and V must be the same.")

# errors for weight matrix
expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = 1, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "'weight' must be NULL or a symmetric matrix.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = diag(2), optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "'weight' must be NULL or a symmetric matrix.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = FALSE, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "'weight' must be NULL or a symmetric matrix.")

# errors for optimal
expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = 1, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "'optimal' must be logical.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = diag(9), lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = NULL), "'optimal' must be logical.")

# errors for lower and upper
expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = 1, upper = NULL, controlDE = NULL, controlBFGS = NULL), "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NA, controlDE = NULL, controlBFGS = NULL), "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = Inf, controlDE = NULL, controlBFGS = NULL), "The lengths of 'lower' and 'upper' must be identical to that of the parameters.")

# errors for controlDE
expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = TRUE, controlBFGS = NULL), "'controlDE' must be correctly specified.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = 1, controlBFGS = NULL), "'controlDE' must be correctly specified.")

# errors for controlBFGS
expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = FALSE), "'controlBFGS' must be correctly specified.")

expect_error(LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V, weight = NULL, optimal = TRUE, lower = NULL, upper = NULL, controlDE = NULL, controlBFGS = diag(10)), "'controlBFGS' must be correctly specified.")
