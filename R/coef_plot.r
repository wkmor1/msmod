#' Plot the coefficients of a multispecies model
#'
#' @param x An object of class "glmerMod", "rjags" or "rjags_parallel"
#' @examples
#' msm_glmer <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
#' coef_plot(msm_glmer)
#' msm_jags <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
#'   method = 'jags')
#' coef_plot(msm_jags)
#' @export

setGeneric(
  "coef_plot",
  function(x) {
    base::standardGeneric("coef_plot")
  }
)

coef_plot_ <- function(coefs) {
  ggplot2::ggplot(
    coefs,
    ggplot2::aes(x = `Coefficient value`, y = coef)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = `Coefficient value` - ci,
        xmax = `Coefficient value` + ci,
        height = 0)
    ) +
    ggplot2::scale_y_discrete(name = "")
}

#' @describeIn coef_plot coefficient plot for glmer model
setMethod(
  "coef_plot",
  base::c(x = "glmerMod"),
  function(x) {
    dplyr::data_frame(
      coef =
        x %>%
        methods::slot("pp") %>%
        base::get("X", envir = .) %>%
        base::colnames(.),
      "Coefficient value" =
        x %>%
        lme4::fixef(.),
      ci =
        x %>%
        lme4::vcov.merMod(.) %>%
        methods::slot("factors") %>%
        base::getElement("correlation") %>%
        methods::slot("sd") %>%
        magrittr::multiply_by(1.96)
    ) %>%
    coef_plot_
  }
)

coef_plot_rjags <- function(x) {
  samples <-
    magrittr::extract2(x, "BUGSoutput") %>%
    magrittr::extract2("sims.list") %>%
    magrittr::extract2("B") %>% 
    base::apply(base::c(1, 3), base::mean)
  dplyr::data_frame(
    coef =
      base::c(
        "(Intercept)",
        base::attr(x, "x")
      ),
    "Coefficient value" =
      base::colMeans(samples),
    ci =
      base::apply(samples, 2, sd) %>%
      magrittr::multiply_by(1.96)
  ) %>%
  coef_plot_
}

#' @describeIn coef_plot coefficient plot for rjags and rjags.parallel objects
setMethod(
  "coef_plot",
  base::c(x = "rjags"),
  coef_plot_rjags
)

#' @describeIn coef_plot coefficient plot for rjags and rjags.parallel objects
setMethod(
  "coef_plot",
  base::c(x = "rjags.parallel"),
  coef_plot_rjags
)

utils::globalVariables(base::c("Coefficient value", "ci"))
