msm_jags <- function(y, sites, x, species, n_species, traits, data, site_re,
                     dots)
{
  model <-
    function() {
      for (site in 1:n_sites) {
        Z[site, 1:n_species] ~ dmnorm(Mu[site, ], Tau)
        for (species in 1:n_species) {
          Mu[site, species] <- inprod(Beta_raw[species, ], X_site[site, ])
          Y[site, species] ~ dbern(step(Z[site, species]))
        }
      }
      for (species in 1:n_species) {
        sigma_[species] <- sqrt(Sigma[species, species])
        B_raw[species, 1:K] ~ dmnorm(B_hat_raw[species, ], Tau_B)
        for (k in 1:K) {
          B_hat_raw[species, k] <- sum(B_hat_raw_[species, k, ])
          Beta_raw[species, k] <- B_raw[species, k] * xi[k]
          Beta[species, k] <- Beta_raw[species, k] / sigma_[species]
          for (j in 1:J) {
            B_hat_raw_[species, k, j] <- G_raw[k, j] * X_sp[species, j]
            Gamma_hat[species, k, j] <- B_hat_raw_[species, k, j] * xi[k] /
              sigma_[species] / X_sp[species, j]
          }
        }
        for (species_ in 1:n_species) {
          Rho[species, species_] <- Sigma[species, species_] / sigma_[species] *
            sigma_[species_]
        }
      }
      for (k in 1:K) {
        xi[k] ~ dunif(0, 100)
        sigma[k] <- sqrt(Sigma_B[k, k]) * xi[k]
        for (k_ in 1:K) {
          Rho_B[k, k_] <- Sigma_B[k, k_] /
            (sqrt(Sigma_B[k, k]) * sqrt(Sigma_B[k_, k_]))
        }
        for (j in 1:J) {
          G_raw_[k, j] ~ dnorm(0, .0001)
          G_raw[k, j] <- ifelse(k == 1 && j > 1, 0, G_raw_[k, j])
          Gamma[k, j] <- mean(Gamma_hat[, k, j])
        }
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
      Tau_B ~ dwish(I_B, df_B)
      Sigma_B <- inverse(Tau_B)
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
  
  X_site <-
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
  
  J <-
    traits %>%
    base::length(.) %>%
    magrittr::add(1)
  
  X_sp <-
    data %>%
    dplyr::distinct_(species) %>%
    dplyr::select_(traits) %>%
    magrittr::inset2('(Intercept)', value = 1) %>%
    base::subset(
      select =
        J %>%
        base::c(
          {.} %>%
            magrittr::subtract(1) %>%
            base::seq_len(.)
        )
    )
  
  inits <-
    dots %>%
    base::getElement('n.chains') %>%
    base::is.null(.) %>%
    base::ifelse(
      yes = 3,
      no =
        dots %>%
        base::getElement('n.chains')
    ) %>%
    base::seq_len(.) %>%
    base::lapply(
      function(x)
      {
        xi <-
          K %>%
          stats::runif(.)
        
        B_raw <-
          n_species %>%
          magrittr::multiply_by(K) %>%
          stats::runif(-1, 1) %>%
          matrix(K, n_species) %>%
          magrittr::divide_by(xi) %>%
          base::t(.)
        
        base::list(
          Z =
            Y %>%
            magrittr::subtract(.5),
          Tau =
            n_species %>%
            base::diag(.),
          B_raw = B_raw,
          xi = xi,
          Tau_B =
            K %>%
            base::diag(.),
          G_raw_ =
            J %>%
            magrittr::multiply_by(K) %>%
            stats::runif(-1, 1) %>%
            matrix(K, J)
        )
      }
    )
  
  parameters.to.save <-
    c('Beta', 'Gamma', 'Rho_B', 'sigma', 'Rho', 'Gamma_hat')
  
  base::list(
    Y = Y,
    X_site = X_site,
    X_sp = X_sp,
    K = K,
    J = J,
    n_species = n_species,
    n_sites =
      X_site %>%
      base::NROW(.),
    I =
      n_species %>%
      base::diag(.),
    df =
      n_species %>%
      magrittr::add(1),
    I_B =
      X_sp %>%
      base::NCOL(.) %>%
      base::diag(.),
    df_B =
      X_sp %>%
      base::NCOL(.) %>%
      magrittr::add(1)
  ) %>%
    base::list(model) %>%
    magrittr::set_names(c('data', 'model.file')) %>%
    bind_if_not_in(dots, 'inits') %>%
    bind_if_not_in(dots, 'parameters.to.save') %>%
    bind_if_not_in(dots, 'n.iter', 200) %>%
    bind_if_not_in(dots, 'DIC', FALSE) %>%
    base::c(dots) %>%
    eval_with_args(R2jags::jags.parallel)
}