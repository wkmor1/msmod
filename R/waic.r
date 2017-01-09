#' Calculate waic
#'
#' @param x An object of class "jagsUI"
#' @examples
#' msm_jags <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
#'   method = 'jags')
#' waic(msm_jags)

#' @export
waic <- function(x) UseMethod("waic")

setGeneric("waic")

waic.jagsUI <-
  function(x) {
    
    data    <- x$model$data()
    n      <- data$n
    nsp    <- data$J
    nsamp  <- x$mcmc.info$n.samples
    
    B <- x$sims.list$B
    X <- data$X
    Y <- data$Y
    
    theta <- vapply(
      1:nsamp, 
      function(x) inprod_mat(X, B[x, ,]),
      matrix(NA_real_, n, nsp)
    )
    
    pY <- apply(theta, 3, pbern_probit, Y = Y)
    
    lpd <- sum(log(apply(pY, 1, mean)))
    
    p_waic <- sum(apply(log(pY), 1, stats::var))
    
    -2 * (lpd - p_waic)
    
  }

pbern_probit <-
  function(x, Y) {
    x <- stats::pnorm(x)
    ifelse(Y, x, 1 - x)
  }

#' @describeIn waic waic for jags model
setMethod(
  "waic",
  base::c(x = "jagsUI"),
  waic.jagsUI
)
