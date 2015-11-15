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
#' @export
msm <- function(y, sites, x, species, traits, data, site_re = FALSE,
  type = c("mstm", "jsdm"),
  method = c("glmer", "jags", "stan"),
  ...) {
  x <- x %>% dplyr::select_vars_(base::names(data), .)

  if (!missing(traits)) {
    traits <- traits %>% dplyr::select_vars_(base::names(data), .)
    x <- x %>% base::c(traits)
  }

  data <-
    data %>%
    dplyr::mutate_each_(
      dplyr::funs(
        base::scale(.) %>%
        magrittr::extract(, 1) %>%
        magrittr::divide_by(2)
        ),
      x
    )

  n_species <-
    data %>%
    dplyr::select(species) %>%
    dplyr::distinct(.) %>%
    base::nrow(.)

  dots <-
    base::list(...)

  base::match.arg(method) %>%
  base::switch(
    glmer = msm_glmer(y, sites, x, species, n_species, traits, data, site_re,
                      type, dots),
    jags  = msm_jags(y, sites, x, species, n_species, data, type, dots),
    stan  = msm_stan(y, sites, x, species, n_species, data , type, dots)
  )
}
