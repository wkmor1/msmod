msm_jags <- function(y, sites, x, species, n_species, data, type, dots)
{
  if (type == 'mstm') {
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
      for (site in 1:n_sites) {
        Z[site, 1:n_species] ~ dmnorm(Mu[site, ], Tau)
        for (species in 1:n_species) {
          Mu[site, species] <- inprod(Beta_raw[species, ], X[site, ])
          Y[site, species] ~ dbern(step(Z[site, species]))
        }
      }
      for (species in 1:n_species) {
        sigma_[species] <- sqrt(Sigma[species, species])
        for (k in 1:K) {
          Beta_raw[species, k] ~ dnorm(mu[k], tau[k])
          Beta[species, k] <- Beta_raw[species, k] / sigma_[species]
        }
        for (species_ in 1:n_species) {
          Rho[species, species_] <-
            Sigma[species, species_] / pow(sigma_[species], 2)
          EnvRho[species, species_] <-
            sum(Beta[species, ] * Beta[species_, ])
        }
      }
      for (k in 1:K) {
        mu[k] ~ dnorm(0, .001)
        tau[k] <- pow(sigma[k], -2)
        sigma[k]  ~ dunif(0, 100)
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
    dplyr::select_(x) %>%
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
        Beta_raw =
          n_species %>%
          magrittr::multiply_by(K) %>%
          stats::rnorm(.) %>%
          base::matrix(n_species, K)
        base::list(
          Z =
            Y %>%
            magrittr::subtract(.5),
          Beta_raw = Beta_raw,
          mu = 
            base::colMeans(Beta_raw),
          Tau =
            n_species %>%
            base::diag(.),
          sigma =
            K %>%
            base::rep(1, .)
        )
      }

  parameters.to.save <-
    c('Beta', 'sigma', 'Rho')

  if (!serial) dots$export_obj_names <- c('n_species', 'K', 'Y') 

  base::list(
    Y = Y,
    X = X,
    K = K,
    n_species = 
      n_species,
    n_sites =
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
