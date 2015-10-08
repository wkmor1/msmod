#' Fit a multispecies model
#'
#' @param y Character. Column name of response variable, species occurrences.
#' @param sites Character. Column name of site variable.
#' @param x Character. Column names of site level predictor variables.
#' @param species Character. Column name of species variable.
#' @param traits Character. Column names of species trait variables.
#' @param data A data.frame containing the variables for the model.
#' @param site_re Logical. Should a site level random effect be included.
#' @param method The used to fit the model.
#' @param ... Further arguments to pass on to model fitting functions.
#' @examples
#' msm_glmer <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
#' msm_glmer_probit <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla',
#'   eucs, family=binomial(link='probit'))
#' msm_jags <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs,
#'   method='jags')
#' @export
msm <- function(y, sites, x, species, traits, data, site_re=FALSE,
  method=c("glmer", "jags", "stan"), ...)
{
  x %<>% dplyr::select_vars_(base::names(data), .)

  traits %<>% dplyr::select_vars_(base::names(data), .)

  data %<>%
  dplyr::mutate_each_(
    dplyr::funs(
      base::scale(.) %>%
      magrittr::extract(, 1) %>%
      magrittr::divide_by(2)
      ),
    base::c(x, traits)
  )

  n_species <-
    data %>%
    dplyr::select(species) %>%
    dplyr::distinct(.) %>%
    base::nrow(.)

  dots <-
    list(...)

  base::match.arg(method) %>%
  base::switch(
    glmer =
      msm_glmer(y, sites, x, species, n_species, traits, data, site_re, dots),
    jags =
      msm_jags(y, sites, x, species, n_species, traits, data, site_re, dots),
    stan =
      msm_stan(y, sites, x, species, n_species, data , dots)
  )
}

msm_glmer <- function(y, sites, x, species, n_species, traits, data, site_re,
  dots)
{
  ' %s ~ %s + %s + (1 + %s | %s)' %>%
  base::paste0(
    site_re %>%
    dplyr::first(.) %>%
    base::ifelse(' + (1 | %s)', '')
  ) %>%
  base::sprintf(
    y,
    x %>%
    base::paste(collapse = ' + '),
    x %>%
    base::expand.grid(traits) %>%
    base::do.call(
      function(...) base::paste(..., sep = ':'), .
    ) %>%
    base::unlist(.) %>%
    base::paste(collapse = ' + '),
    x %>%
    base::paste(collapse = ' + '),
    species,
    sites
  ) %>%
  stats::formula(.) %>%
  base::list(data) %>%
  magrittr::set_names(c('formula', 'data')) %>%
  bind_if_not_in(dots, 'family', stats::binomial) %>%
  base::c(dots) %>%
  eval_with_args(lme4::glmer)
}

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
  eval_with_args(R2jags::jags)
}

msm_stan <- function(y, sites, x, species, n_species, data, dots)
{ 
  model <- "
    functions {
     int sum(int[,] a) {
        int s;
        s <- 0;
        for (i in 1:size(a))
          for (j in 1:size(a[i]))
            s <- s + a[i, j];
        return s;
      }
    }
    
    data {
      int<lower=1> K;
      int<lower=1> D;
      int<lower=0> N;
      int<lower=0, upper=1> y[N, D];
      vector[K] x[N];
    }
    
    transformed data {
      int<lower=0> N_pos;
      int<lower=1, upper=N> n_pos[sum(y)];
      int<lower=1, upper=D> d_pos[size(n_pos)];
      int<lower=0> N_neg;
      int<lower=1,upper=N> n_neg[(N * D) - size(n_pos)];
      int<lower=1,upper=D> d_neg[size(n_neg)];
    
      N_pos <- size(n_pos);
      N_neg <- size(n_neg);
    
      {
        int i;
        int j;
        i <- 1;
        j <- 1;
        for (n in 1:N) {
          for (d in 1:D) {
            if (y[n, d] == 1) {
              n_pos[i] <- n;
              d_pos[i] <- d;
              i <- i + 1;
            } else {
              n_neg[j] <- n;
              d_neg[j] <- d;
              j <- j + 1;
            }
          }
        }
      }
    }
    
    parameters {
      matrix[D, K] beta;
      cholesky_factor_corr[D] L_Omega;
      vector<lower=0>[N_pos] z_pos;
      vector<upper=0>[N_neg] z_neg;
    }
    
    transformed parameters {
      vector[D] z[N];
      for (n in 1:N_pos)
        z[n_pos[n], d_pos[n]] <- z_pos[n];
      for (n in 1:N_neg)
        z[n_neg[n], d_neg[n]] <- z_neg[n];
    }
    
    model {
      vector[D] beta_x[N];
      to_vector(beta) ~ cauchy(0, 2.5);
      L_Omega ~ lkj_corr_cholesky(1);
      for (n in 1:N)
        beta_x[n] <- beta * x[n];
      z ~ multi_normal_cholesky(beta_x, L_Omega);
    }
    
    generated quantities {
      corr_matrix[D] Omega;
      Omega <- multiply_lower_tri_self_transpose(L_Omega);  
    }
  "

  y %<>%
    dplyr::select_(data, .) %>%
    base::unlist(.) %>%
    base::as.integer(.) %>%
    base::matrix(ncol = n_species)
  
  K <-
    x %>%
    base::length(.) %>%
    magrittr::add(1)

  x <-
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
    ) %>%
    base::as.matrix(.)

  N <- base::nrow(y)
  D <- n_species
  
  old <- base::options(mc.cores = parallel::detectCores())
  base::on.exit(base::options(old), add = TRUE)
  
  base::list(K = K, D = D, N = N, y = y, x = x) %>%
  base::list(model) %>%
  magrittr::set_names(c('data', 'model_code')) %>%  
  bind_if_not_in(dots, 'chains', 4) %>%
  bind_if_not_in(dots, 'iter', 400) %>%
  base::c(dots) %>%
  eval_with_args(rstan::stan)
}

c("n_sites", "inprod", "Z", "pow", "inverse",
  "Tau", "B_raw", "G_raw_", "xi", "Tau_B") %>%
utils::globalVariables(.)
