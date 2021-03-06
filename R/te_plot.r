#' Plot the relationship between a trait and the environment infered from a multispecies model
#'
#' @param x An object of class 'glmerMod'
#' @param x_var Character. Name of site-level predictor variable.
#' @param trait Character. Name of species trait variable.
#' @param nsims Integer. Number of simulations for approximate Bayesian estimates of credible intervals and posterior densities.
#' @param rlim Numeric vector. The x axis limits of the confidence ribbon.
#' @examples
#' msm_fit <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
#' te_plot(msm_fit, 'logit_rock', 'ln_sla')
#' @export

setGeneric(
  "te_plot",
  function(x, x_var, trait, nsims = 200, rlim = c(-2, 2)) {
    base::standardGeneric("te_plot")
  }
)

#' @describeIn te_plot Trait environment plot for glmer model
setMethod(
  "te_plot",
  base::c(x = "glmerMod"),
  function(x, x_var, trait, nsims, rlim = c(-2, 2)) {

    frame <-
      x %>%
      methods::slot("frame")

    species <-
      frame %>%
      base::colnames(.) %>%
      base::rev(.) %>%
      dplyr::nth(
        x %>%
        methods::slot("flist") %>%
        base::length(.)
      )

    fixefs <-
      x %>%
      lme4::fixef(.)

    int <- 
      grep(
        sprintf("%1$s:%2$s|%2$s:%1$s", x_var, trait),
        names(fixefs),
        value = TRUE
      )

    k <-
      fixefs %>%
      base::names(.) %>%
      magrittr::is_in(
        x_var %>%
          c(int)
      )

    fixefs <-
      fixefs %>%
      magrittr::extract(k)

    se_fixefs <-
      x %>%
      lme4::vcov.merMod(.) %>%
      methods::slot("factors") %>%
      base::getElement("correlation") %>%
      methods::slot("sd") %>%
      magrittr::extract(k) %>%
      magrittr::set_names(
        fixefs %>%
        base::names(.)
      )

    fixef_x <-
      fixefs %>%
      magrittr::extract(x_var)

    fixef_int <-
      fixefs %>%
      magrittr::extract(int)

    se_fixef_x <-
      se_fixefs %>%
      magrittr::extract(x_var)

    se_fixef_int <-
      se_fixefs %>%
      magrittr::extract(int)

    trait_frame <-
      frame %>%
      dplyr::distinct_(species, trait) %>%
      base::lapply(base::as.vector) %>%
      dplyr::as_data_frame(.) %>%
      dplyr::arrange(species)

    trait_values <-
      trait_frame %>%
      dplyr::select_(trait) %>%
      base::unlist(.)

    species_names <-
      trait_frame %>%
      dplyr::select_(species) %>%
      base::unlist(.)

    xvals <-
      base::seq(rlim[1], rlim[2], length.out = 201)

    yvals <-
      nsims %>%
      stats::rnorm(fixef_x, se_fixef_x) %>%
      magrittr::add(
        nsims %>%
        stats::rnorm(fixef_int, se_fixef_int) %>%
        base::outer(xvals)
      ) %>%
      base::as.data.frame(.) %>%
      base::vapply(
        stats::quantile,
        2 %>%
        base::numeric(.),
        probs = .025 %>%
        base::c(.975)
      ) %>%
      base::t(.)

    sim <-
      x %>%
      arm::sim(nsims)

    nsims <-
      nsims %>%
      magrittr::multiply_by(.95) %>%
      base::as.integer(.)

    sim_fixef <-
      sim %>%
      methods::slot("fixef")

    sim <-
      sim %>%
      methods::slot("ranef") %>%
      magrittr::extract2(species) %>%
      magrittr::extract(, , x_var) %>%
      magrittr::add(
        sim_fixef %>%
        magrittr::extract(, int) %>%
        base::outer(trait_values)
      ) %>%
      magrittr::add(
        sim_fixef %>%
        magrittr::extract(, x_var)
      ) %>%
      magrittr::extract(
        nsims %>%
        base::seq_len(.) %>%
        magrittr::add(
          nsims %>%
          magrittr::multiply_by(.025) %>%
          magrittr::divide_by(.95) %>%
          as.integer(.)
        ),
      ) %>%
      base::as.vector(.) %>%
      dplyr::data_frame(
        response = .,
        trait_values =
          trait_values %>%
          base::rep(each = nsims),
        species =
          species_names %>%
          base::rep(each = nsims)
      ) %>%
      dplyr::group_by_(species)

    ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data =
        base::cbind(xvals, yvals) %>%
        base::as.data.frame(.),
      ggplot2::aes(
        x = xvals,
        ymin = `2.5%`,
        ymax = `97.5%`
      ),
      alpha = .3
    ) +
    ggplot2::coord_cartesian(
       xlim = c(-2, 2)
    ) +
    ggplot2::xlab(trait) +
    ggplot2::ylab(
      "Response to " %>% base::paste0(x_var)
    ) +
    ggplot2::geom_abline(
      intercept =
        fixefs %>%
        magrittr::extract(x_var),
      slope =
        fixefs %>%
        magrittr::extract(int)
    ) +
    ggplot2::geom_text(
      data =
        sim %>%
        dplyr::summarise(
          response =
            response %>%
            base::mean(.)
        ) %>%
        dplyr::bind_cols(
         dplyr::data_frame(trait_values)
        ),
      mapping =
        ggplot2::aes(
          x = trait_values,
          y = response,
          label = species
        ),
      alpha = .9
    ) +
    ggplot2::geom_violin(
      data = sim,
      ggplot2::aes(
        x = trait_values,
        y = response,
        group = species
      ),
      scale = "width",
      width = .11,
      position = "identity",
      alpha = .2,
      adjust = 2
    )
  }
)

utils::globalVariables(c("2.5%", "97.5%", "response"))
