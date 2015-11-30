msm_jags <- function(y, sites, x, species, n_species, data, type, dots, model) {
  if (identical(type, "mstm")) {
    stop("mstm models currently only implemented with method = \"glmer\"")
  }

  serial <- is.null(dots$n.cluster)

  if (serial) {
    jags_fn <- "R2jags::jags"
  } else {
    jags_fn <- "R2jags::jags.parallel"
    dots$export_obj_names <- base::c("J", "K", "Y", "X", "I", "df", "n")
    dots$envir <- base::environment()
  }

  J <- n_species

  Y <-
    y %>%
    dplyr::select_(data, .) %>%
    base::unlist(.) %>%
    base::matrix(ncol = n_species)

  K <-
    x %>%
    base::length(.) %>%
    magrittr::add(1)

  X <-
    data %>%
    dplyr::distinct_(sites) %>%
    dplyr::select_(.dots = x) %>%
    magrittr::inset2("(Intercept)", value = 1) %>%
    base::subset(
      select =
        K %>%
        base::c(
          {.} %>%
            magrittr::subtract(1) %>%
            base::seq_len()
        )
    ) %>%
    base::as.matrix()

  df <-
    n_species %>%
    magrittr::add(1)

  I <-
    n_species %>%
    base::diag()

  n <-
    X %>%
    base::NROW()

  inits <- function() {
    Tau <-
      stats::rWishart(1, df, I)[, , 1]
    Sigma <-
      base::solve(Tau)
    Z <-
      base::rep(0, n_species) %>%
      MASS::mvrnorm(1, ., Sigma) %>%
      base::replicate(n, .) %>%
      base::t() %>%
      base::abs() %>% {
        base::ifelse(base::as.matrix(Y), ., -.)
      }
    Sigma <- stats::cov(Z)
    B <- base::suppressWarnings({
      base::sapply(
        base::seq_len(base::ncol(Y)),
        function(x) {
          stats::glm(
            Y[, x] ~ X[, -1],
            family = stats::binomial(link = "probit")
          ) %>%
          stats::coef() %>%
          base::unname()
        }
      ) %>% base::t()
    })
    B_raw <- B * base::sqrt(base::diag(Sigma))
    mu <- base::apply(B_raw, 2, mean)
    sigma <- base::pmin(99, base::apply(B_raw, 2, sd))
    Tau <- base::solve(Sigma)
    list(Tau = Tau, Z = Z, B_raw = B_raw, mu = mu, sigma = sigma)
  }

  parameters.to.save <-
    base::c("B", "sigma", "Rho", "EnvRho")

  result <-
    base::list(
      Y    = Y,
      X    = X,
      covx = stats::cov(X),
      K    = K,
      J    = J,
      n    = n,
      I    = I,
      df   = df
    ) %>%
    base::list(model) %>%
    magrittr::set_names(base::c("data", "model.file")) %>%
    bind_if_not_in(dots, "inits") %>%
    bind_if_not_in(dots, "parameters.to.save") %>%
    bind_if_not_in(dots, "n.iter", 200) %>%
    bind_if_not_in(dots, "DIC", FALSE) %>%
    base::c(dots) %>%
    eval_with_args(jags_fn)

  base::attr(result, "x") <- base::unname(x)
  
  return(result)
  
}
