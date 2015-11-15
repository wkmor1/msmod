msm_glmer <- function(y, sites, x, species, n_species, traits, data, site_re, type, dots)
{
  if (identical(type, "jsdm")) {
    stop('mstm models currently only implemented with method = "jags" or method = "stan"')
  }

  " %s ~ %s + %s + (1 + %s | %s)" %>%
    base::paste0(
      site_re %>%
        dplyr::first(.) %>%
        base::ifelse(" + (1 | %s)", "")
    ) %>%
    base::sprintf(
      y,
      x %>%
        base::paste(collapse = " + "),
      x %>%
        base::expand.grid(traits) %>%
        base::do.call(
          function(...) base::paste(..., sep = ":"), .
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
    bind_if_not_in(dots, 'control',
      lme4::glmerControl(optimizer = "bobyqa")) %>%
    base::c(dots) %>%
    eval_with_args("lme4::glmer")
}
