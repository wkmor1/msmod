# the jsdm in jags

jsdm_jags <- function() {
  model.file = tempfile()
  cat("model {
    for (i in 1:n) {
      Z[i, 1:J] ~ dmnorm(Mu[i, ], Tau)
      for (j in 1:J) {
        Mu[i, j] <- inprod(B_raw[j, ], X[i, ])
        Y[i, j] ~ dbern(step(Z[i, j]))
      }
    }
    for (j in 1:J) {
      sigma_[j] <- sqrt(Sigma[j, j])
      env_sigma_[j] <- sqrt(EnvSigma[j, j])
      for (k in 1:K) {
        B_raw[j, k] ~ dnorm(mu[k], tau[k])
        B[j, k] <- B_raw[j, k] / sigma_[j]
      }
      for (j_ in 1:J) {
        Rho[j, j_] <- Sigma[j, j_] / (sigma_[j] * sigma_[j_])
        EnvSigma[j, j_] <- sum(EnvSigma1[, j, j_]) + sum(EnvSigma2[, , j, j_])
        EnvRho[j, j_] <- EnvSigma[j, j_] / (env_sigma_[j] * env_sigma_[j_])
        for (k in 2:K) {
          EnvSigma1[k - 1, j, j_] <- B[j, k] * B[j_, k]
          for (k_ in 2:K) {
            EnvSigma2[k - 1, k_ - 1, j, j_] <-
              B[j, k] * B[j_, k_] * ifelse(k_ != k, covx[k, k_], 0)
          }
        }
      }
    }
    for (k in 1:K) {
      mu[k] ~ dnorm(0, 1)
      tau[k] <- pow(sigma[k], -2)
      sigma[k] ~ dnorm(0, 1)T(0,)
    }
    Tau ~ dwish(I, df)
    Sigma <- inverse(Tau)
  }", file = model.file)
  model.file
}

