msm_jags <- function(y, sites, x, species, n_species, data, type, dots)
{
  if (identical(type, 'mstm')) {
    stop('mstm models currently only implemented with method = "glmer"')
  }

  serial <- is.null(dots$n.cluster)

  if (serial) {
    jags_fn <- "R2jags::jags"
  } else {
    jags_fn <- "R2jags::jags.parallel"
    dots$export_obj_names <- base::c('J', 'K', 'Y', 'X', 'I', 'df', 'n')
    dots$envir <- base::environment()
  }

  model <-
    function() {
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
        mu[k] ~ dnorm(0, .001)
        tau[k] <- pow(sigma[k], -2)
        sigma[k] ~ dunif(0, 100)
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
    }

  J <- n_species

  Y <-
    y %>%
    dplyr::select_(data, .) %>%
    base::unlist(.) %>%
    base::matrix(ncol = n_species)

  K <-
    x %>%
    base::length(.) %>%
    magrittr::add(1)

  X <-
    data %>%
    dplyr::distinct_(sites) %>%
    dplyr::select_(.dots = x) %>%
    magrittr::inset2('(Intercept)', value = 1) %>%
    base::subset(
      select =
        K %>%
        base::c(
          {.} %>%
            magrittr::subtract(1) %>%
            base::seq_len()
        )
    ) %>%
    base::as.matrix()

  df <-
    n_species %>%
    magrittr::add(1)

  I <-
    n_species %>%
    base::diag()

  n <-
    X %>%
    base::NROW()

  inits <- function() {
    Tau <-
      stats::rWishart(1, df, I)[, ,1]
    Sigma <-
      base::solve(Tau)
    Z <-
      base::rep(0, n_species) %>%
      MASS::mvrnorm(1, ., Sigma) %>%
      base::replicate(n, .) %>%
      base::t() %>%
      base::abs() %>%
      {base::ifelse(base::as.matrix(Y), ., -.)}
    Sigma <-
      mclust::mvnXXX(Z)$parameters$variance$sigma[, , 1]
    B <- base::suppressWarnings({
      base::sapply(
        base::seq_len(base::ncol(Y)),
        function(x) {
          stats::glm(
            Y[, x] ~ X[, -1],
            family = stats::binomial(link = probit)
          ) %>%
          stats::coef() %>%
          base::unname()
        }
      ) %>% base::t()
    })
    B_raw <- B * sqrt(diag(Sigma))
    mu <- apply(B_raw, 2, mean)
    sigma <- apply(B_raw, 2, sd)
    list(Tau = solve(Sigma), Z = Z, B_raw = B_raw, mu = mu, sigma = sigma)
  }

  parameters.to.save <-
    base::c('B', 'sigma', 'Rho', 'EnvRho')

  result <-
    base::list(
      Y    = Y,
      X    = X,
      covx = stats::cov(X),
      K    = K,
      J    = J,
      n    = n,
      I    = I,
      df   = df
    ) %>%
    base::list(model) %>%
    magrittr::set_names(base::c('data', 'model.file')) %>%
    bind_if_not_in(dots, 'inits') %>%
    bind_if_not_in(dots, 'parameters.to.save') %>%
    bind_if_not_in(dots, 'n.iter', 200) %>%
    bind_if_not_in(dots, 'DIC', FALSE) %>%
    base::c(dots) %>%
    eval_with_args(jags_fn)

  return(result)
}
