#' Fit a Multispecies Model
#'
#' @param y Character. Column name of response variable, species occurrences.
#' @param sites Character. Column name of site variable.
#' @param x Character. Column names of site level predictor variables.
#' @param species Character. Column name of species variable.
#' @param traits Character. Column names of species trait variables. Ignore
#' @param data A data.frame containing the variables for the model.
#' @param site_re Logical. Should a site level random effect be included.
#' @param type The type of model to fit either 'mstm', multispecies trait model,
#' or 'jsdm', joint species distribution model.
#' @param method The method used to fit the model.
#' @param ... Further arguments to pass on to model fitting functions.
#' @examples
#' msm_glmer <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
#' msm_glmer_probit <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla',
#'   eucs, family = binomial(link = 'probit'))
#' msm_jags <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
#'   method = 'jags')
#' msm_stan <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
#'   method = 'stan')
#' @export
msm <- function(y, sites, x, species, traits, data, site_re = FALSE,
  type = c("mstm", "jsdm"),
  method = c("glmer", "jags", "stan"),
  ...)
{
  x %<>% dplyr::select_vars_(base::names(data), .)

  unscaled_data <- x
  
  if (!missing(traits)) {
    traits %<>% dplyr::select_vars_(base::names(data), .)
    unscaled_data %<>% base::c(traits)
  }
    
  data %<>%
    dplyr::mutate_each_(
      dplyr::funs(
        base::scale(.) %>%
        magrittr::extract(, 1) %>%
        magrittr::divide_by(2)
        ),
      unscaled_data
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
      msm_glmer(y, sites, x, species, n_species, traits, data, site_re, type, dots),
    jags =
      msm_jags(y, sites, x, species, n_species, data, type, dots),
    stan =
      msm_stan(y, sites, x, species, n_species, data , type, dots)
  )
}

c("n_sites", "inprod", "Z", "pow", "inverse",
  "Tau", "B_raw", "G_raw_", "xi", "Tau_B") %>%
utils::globalVariables(.)
