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
  method=c('glmer', 'jags'), ...) {

  method <- base::match.arg(method)

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
  
  method %>%
  switch(
    glmer=
      msm_glmer(y, sites, x, species, n_species, traits, data, site_re, ...),
    jags=
      msm_jags(y, sites, x, species, n_species, traits, data, site_re, ...)

  )
}

msm_glmer <- function(y, sites, x, species, n_species, traits, data, site_re,
  ...) {

  args <-
    ' %s ~ %s + %s + (1 + %s | %s)' %>%
    base::paste0(
      site_re %>% 
      dplyr::first(.) %>% 
      base::ifelse(' + (1 | %s)', '')
    ) %>%
    base::sprintf(
      y,
      x %>% 
      base::paste(collapse=' + '),
      x %>% 
      base::expand.grid(traits) %>%
      base::do.call(
        function(...) base::paste(..., sep=':'), .
      ) %>%
      base::unlist(.) %>%
      base::paste(collapse=' + '),
      x %>%
      base::paste(collapse=' + '),
      species,
      sites
    ) %>%
    stats::formula(.) %>%
    base::list() %>%
    magrittr::inset2('data', data) %>%
    magrittr::inset2('family', stats::binomial)

  args %>%
  return_if_not(
    list(...) %>%
    magrittr::extract2('family') %>%
    base::is.null(.),
    args %>%
    magrittr::inset2('family', NULL)
  ) %>%
  base::do.call(lme4::glmer, .)
}

msm_jags <- function(y, sites, x, species, n_species, traits, data, site_re,
  ...) {

  Y <-
    y %>%
    dplyr::select_(data, .) %>%
    base::unlist(.) %>%
    base::matrix(ncol=n_species)

  X <-
    data %>%
    dplyr::distinct_(sites) %>%
    dplyr::select_(x) %>%
    magrittr::inset2('(Intercept)', value=1) %>%
    base::rev(.)
  
  K <-
    X %>%
    base::ncol(.)

  n_sites <-
    X %>%
    base::nrow(.)
  
  I <-
    n_species %>%
    base::diag(.)

  df <-
    n_species %>%
    magrittr::add(1)

  inits <-
    base::list(Tau=I) %>%
    magrittr::inset2(
      'Z',
       Y %>%
       magrittr::subtract(.5)
    ) %>%
    base::list(.) %>%
    base::rep(3)
  
  (function() {
    for (site in 1:n_sites) {
      Z[site, 1:n_species] ~ dmnorm(Mu[site, ], Tau[, ])
      for (species in 1:n_species) {
        Mu[site, species] <- 
          inprod(Beta_raw[species, ], X[site, ])
        Y[site, species] ~ dbern(P[site, species])
        P[site, species] <- step(Z[site, species])
      }
    }
    for (species in 1:n_species) {
      for (k in 1:K) {
        Beta_raw[species, k] ~ dnorm(mu_raw[k], tau[k])
      }
    }
    for (k in 1:K) {
      mu_raw[k] ~ dnorm(0, .0001)
      tau[k] <- pow(sigma_raw[k], -2)
      sigma_raw[k] ~ dunif(0, 100)
    }
    Tau[1:n_species, 1:n_species] ~ dwish(I[, ], df)
  }) %>%
  R2jags::jags(
    data=c('Y', 'X', 'K', 'n_sites', 'n_species', 'I', 'df'),
    inits=inits,
    parameters.to.save=c('Beta_raw', 'Tau'), 
    model.file=.,
    ...
  )

}

utils::globalVariables(c(".", "x_"))
