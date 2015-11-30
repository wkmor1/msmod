#' Calculate waic
#'
#' @param x An object of class "rjags", or "rjags.parallel"
#' @examples
#' msm_jags <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
#'   method = 'jags')
#' waic(msm_jags)
#' @export

setGeneric(
  "waic",
  function(x) {
    base::standardGeneric("waic")
  }
)

waic_rjags <-
  function(x) {
    
    data    <- x$model$data()
    samples <- x$BUGSoutput

    n      <- data$n
    nsp    <- data$J
    nsamp  <- samples$n.sims
    
    B <- samples$sims.list$B
    X <- data$X
    Y <- data$Y

    theta <- base::array(base::numeric(), dim = base::c(n, nsp, nsamp))
         
    for (i in 1:n) {
      for (j in 1:nsp) {
         theta[i, j, ] <-
           base::sapply(1:nsamp, function(k) B[k, j, ] %*% X[i, ])
      }
    }

    pY <-
      base::apply(
        theta,
        3,
        function(x) {
          x <- stats::pnorm(x)
          base::ifelse(Y, x, 1 - x)
        }
      )
    
    lpd <- base::sum(base::log(base::apply(pY, 1, base::mean)))
    
    p_waic <- base::sum(base::apply(base::log(pY), 1, stats::var))
    
    -2 * (lpd - p_waic)
    
  }

#' @describeIn waic calculate waic for rjags and rjags.parallel objects
setMethod(
  "waic",
  base::c(x = "rjags"),
  waic_rjags
)

#' @describeIn waic calculate waic for rjags and rjags.parallel objects
setMethod(
  "waic",
  base::c(x = "rjags.parallel"),
  waic_rjags
)
