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
    
    theta <- vapply(
      1:nsamp, 
      function(x) inprod_mat(X, B[x, ,]),
      matrix(NA_real_, n, nsp)
    )
    
    pY <- apply(theta, 3, pbern_probit, Y = Y)
    
    lpd <- sum(log(apply(pY, 1, mean)))
    
    p_waic <- sum(apply(log(pY), 1, var))
    
    -2 * (lpd - p_waic)
    
  }

pbern_probit <-
  function(x, Y) {
    x <- pnorm(x)
    ifelse(Y, x, 1 - x)
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
