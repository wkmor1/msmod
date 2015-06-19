#' Plot the coefficients of a multispecies model
#'
#' @param x An object of class 'glmerMod'
#' @examples
#'
#' library(msm)
#' data(eucs)
#' msm_fit <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
#' coef_plot(msm_fit) 
#'
#' @export

coef_plot <- function(x) {

  x %>%
  is_glmerMod(.)

  dplyr::data_frame(
    coef=
      x %>%
      methods::slot('pp') %>%
      base::get('X', envir=.) %>%
      base::colnames(.),
    `Coefficient value`=
      x %>%
      lme4::fixef(.),
    ci=
      x %>%
      lme4::vcov.merMod(.) %>%
      methods::slot('factors') %>%
      base::getElement('correlation') %>%
      methods::slot('sd') %>%
      magrittr::multiply_by(1.96)
  ) %>%

  ggplot2::ggplot(
    ggplot2::aes(x=`Coefficient value`, y=coef)
  ) + 
  ggplot2::geom_point() +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin=`Coefficient value` - ci,
      xmax=`Coefficient value` + ci,
      height=0)  
  ) +
  ggplot2::scale_y_discrete(name='')
}

utils::globalVariables(c(".", "Coefficient value", "ci"))
