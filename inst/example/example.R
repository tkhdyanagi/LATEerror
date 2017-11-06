library(LATEerror)
library(dplyr)

# simulated data----------------------------------------------------------------
n <- 2000
set.seed(102)
mydata <- data.frame(U1 = rnorm(n, 0, 1), U2 = rnorm(n, 0, 1))
mydata <- mydata %>%
  dplyr::mutate(Z = rbinom(n, size = 1, prob = 0.5)) %>%
  dplyr::mutate(V = rbinom(n, size = 1, prob = 0.5) + 1) %>%
  dplyr::mutate(Ds = ifelse(( -2 + Z + V - U1 >= 0), 1, 0)) %>%
  dplyr::mutate(D = Ds * rbinom(n, size = 1, prob = 0.75) + (1 - Ds) * rbinom(n, size = 1, prob = 0.25)) %>%
  dplyr::mutate(Y = 1 + Ds + U2)

# GMM estimation----------------------------------------------------------------
set.seed(102)
gmmest <- LATEerror(Y = mydata$Y, D = mydata$D, Z = mydata$Z, V = mydata$V)

# pilot GMM estimation result---------------------------------------------------
# estimates
gmmest$pilot$estimate

# standard errors
gmmest$pilot$se

# 95% confidence intervals
gmmest$pilot$ci

# optimal GMM estimation result-------------------------------------------------
# estimates
gmmest$optim$estimate

# standard errors
gmmest$optim$se

# 95% confidence intervals
gmmest$optim$ci

# over-identification test
gmmest$optim$overidentification
