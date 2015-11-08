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
        for (k in 1:K) {
          B_raw[j, k] ~ dnorm(mu[k], tau[k])
          B[j, k] <- B_raw[j, k] / sigma_[j]
        }
        for (j_ in 1:J) {
          Rho[j, j_] <- Sigma[j, j_] / pow(sigma_[j], 2)
          for (k in 1:K) {
            for (k_ in 1:K) {
              num_[k, k_, j, j_] <- B[j, k] * B[j_, k_] *
                ifelse(k_ != k, covx[k, k_], 0)
              dem1_[k, k_, j, j_] <- B[j, k] * B[j, k_] *
                ifelse(k_ != k, covx[k, k_], 0)
              dem2_[k, k_, j, j_] <- B[j_, k] * B[j_, k_] *
                ifelse(k_ != k, covx[k, k_], 0)
            }
            num[k, j, j_] <- B[j, k] * B[j_, k] * sum(num_[, , j, j_])
            dem1[k, j, j_] <- B[j, k] * B[j, k] * sum(dem1_[, , j, j_])
            dem2[k, j, j_] <- B[j_, k] * B[j_, k] * sum(dem2_[, , j, j_])
          }
          EnvRho[j, j_] <- sum(num[, j, j_]) /
            sqrt(sum(dem1[, j, j_]) * sum(dem2[, j, j_]))
        }
      }
      for (k in 1:K) {
        mu[k] ~ dnorm(0, .001)
        tau[k] <- pow(sigma[k], -2)
        sigma[k] ~ dunif(0, 100)
        for (k_ in 1:K) {
          for (i in 1:n) {
            covx_[i, k, k_] <-
              (X[i, k] - mean(X[, k])) * (X[i, k_] - mean(X[, k_]))
          }
          covx[k, k_] <- mean(covx_[, k, k_])
        }
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
    }

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
            base::seq_len(.)
        )
    )

  inits <-
    function(x)
      {
        B_raw =
          n_species %>%
          magrittr::multiply_by(K) %>%
          stats::rnorm(.) %>%
          base::matrix(n_species, K)
        base::list(
          Z =
            Y %>%
            magrittr::subtract(.5),
          B_raw = B_raw,
          mu =
            base::colMeans(B_raw),
          Tau =
            n_species %>%
            base::diag(.),
          sigma =
            K %>%
            base::rep(1, .)
        )
      }

  parameters.to.save <-
    c('B', 'sigma', 'Rho', 'EnvRho')

  if (!serial) dots$export_obj_names <- c('n_species', 'K', 'Y')

  base::list(
    Y = Y,
    X = X,
    K = K,
    J =
      n_species,
    n =
      X %>%
      base::NROW(.),
    I =
      n_species %>%
      base::diag(.),
    df =
      n_species %>%
      magrittr::add(1)
  ) %>%
  base::list(model) %>%
  magrittr::set_names(c('data', 'model.file')) %>%
  bind_if_not_in(dots, 'inits') %>%
  bind_if_not_in(dots, 'parameters.to.save') %>%
  bind_if_not_in(dots, 'n.iter', 200) %>%
  bind_if_not_in(dots, 'DIC', FALSE) %>%
  base::c(dots) %>%
  eval_with_args(jags_fn)
}
