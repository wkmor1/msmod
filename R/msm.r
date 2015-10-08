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

c("n_sites", "inprod", "Z", "pow", "inverse",
  "Tau", "B_raw", "G_raw_", "xi", "Tau_B") %>%
utils::globalVariables(.)
