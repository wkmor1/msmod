msm_stan <- function(y, sites, x, species, n_species, data, type, dots) {
  if (identical(type, "mstm")) {
    stop("mstm models currently only implemented with method = \"glmer\"")
  }

  model <- "
    functions {
      int sum(int[,] a) {
        int s;
        s <- 0;
         for (i in 1:size(a))
          for (j in 1:size(a[i]))
            s <- s + a[i, j];
        return s;
      }
    }

    data {
      int<lower=1> K;
      int<lower=1> D;
      int<lower=0> N;
      int<lower=0, upper=1> y[N, D];
      vector[K] x[N];
    }

    transformed data {
      int<lower=0> N_pos;
      int<lower=1, upper=N> n_pos[sum(y)];
      int<lower=1, upper=D> d_pos[size(n_pos)];
      int<lower=0> N_neg;
      int<lower=1,upper=N> n_neg[(N * D) - size(n_pos)];
      int<lower=1,upper=D> d_neg[size(n_neg)];

      N_pos <- size(n_pos);
      N_neg <- size(n_neg);

      {
        int i;
        int j;
        i <- 1;
        j <- 1;
        for (n in 1:N) {
          for (d in 1:D) {
            if (y[n, d] == 1) {
              n_pos[i] <- n;
              d_pos[i] <- d;
              i <- i + 1;
            } else {
              n_neg[j] <- n;
              d_neg[j] <- d;
              j <- j + 1;
            }
          }
        }
      }
    }

    parameters {
      matrix[D, K] beta;
      cholesky_factor_corr[D] L_Omega;
      vector<lower=0>[N_pos] z_pos;
      vector<upper=0>[N_neg] z_neg;
    }

    transformed parameters {
      vector[D] z[N];
      for (n in 1:N_pos)
        z[n_pos[n], d_pos[n]] <- z_pos[n];
      for (n in 1:N_neg)
        z[n_neg[n], d_neg[n]] <- z_neg[n];
    }

    model {
      vector[D] beta_x[N];
      to_vector(beta) ~ cauchy(0, 2.5);
      L_Omega ~ lkj_corr_cholesky(1);
      for (n in 1:N)
        beta_x[n] <- beta * x[n];
        z ~ multi_normal_cholesky(beta_x, L_Omega);
    }

    generated quantities {
      corr_matrix[D] Omega;
      Omega <- multiply_lower_tri_self_transpose(L_Omega);
    }
  "
  y <-
    y %>%
    dplyr::select_(data, .) %>%
    base::unlist(.) %>%
    base::as.integer(.) %>%
    base::matrix(ncol = n_species)

  K <-
    x %>%
    base::length(.) %>%
    magrittr::add(1)

  x <-
    data %>%
    dplyr::distinct_(sites) %>%
    dplyr::select_(x) %>%
    magrittr::inset2("(Intercept)", value = 1) %>%
    base::subset(
      select =
        K %>%
        base::c(
          {.} %>%
            magrittr::subtract(1) %>%
            base::seq_len(.)
        )
    ) %>%
    base::as.matrix(.)

  N <- base::nrow(y)
  D <- n_species

  base::list(K = K, D = D, N = N, y = y, x = x) %>%
  base::list(model) %>%
  magrittr::set_names(c("data", "model_code")) %>%
  bind_if_not_in(dots, "chains", 4) %>%
  bind_if_not_in(dots, "iter", 400) %>%
  base::c(dots) %>%
  eval_with_args("rstan::stan")
}
