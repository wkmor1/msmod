#' Fit a multispecies model
#'
#' @param y Character. Column name of response variable, species occurrences.
#' @param sites Character. Column name of site variable.
#' @param x Character. Column names of site level predictor variables.
#' @param species Character. Column name of species variable.
#' @param traits Character. Column names of species trait variables.
#' @param data A data.frame containing the variables for the model.
#' @param site_re Logical. Should a site level random effect be included.
#' @export
msm <- function (y, sites, x, species, traits, data, site_re=FALSE) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`

  x %<>% dplyr::select_vars_(base::names(data), .)
  traits %<>% dplyr::select_vars_(base::names(data), .)

  data %<>%

  dplyr::mutate_each_(dplyr::funs(base::scale(.) / 2), base::c(x, traits))

  n_species <-
    data %>%
    dplyr::select(species) %>%
    dplyr::distinct(.) %>%
    base::nrow(.)

  ' %s ~ %s + %s + (1 + %s | %s)' %>%

  base::paste0(
    site_re %>% dplyr::first(.) %>% base::ifelse(' + (1 | %s)', '')
  ) %>%

  base::sprintf(
    y,
    x %>% base::paste(collapse=' + '),
    x %>% base::expand.grid(traits) %>%
      base::do.call(function(...) base::paste(..., sep=':'), .) %>%
      base::unlist(.) %>%
      base::paste(collapse=' + '),
    x %>% base::paste(collapse=' + '),
    species,
    sites
  ) %>%

  stats::formula(.) %>%

  lme4::glmer(
    data=data,
    family  = stats::binomial,
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl   = base::list(
        maxfun = (n_species^2 / 2)^2 * 10 + 1
      )
    )
  )
}

utils::globalVariables(c(".", "x_"))
