#' Eucalyptus species data from the Grampians, Victoria.
#'
#' A dataset containing the occurences, predictor variables and traits of
#' species at sites in the Grampians.
#'
#' @format A data frame with 9160 rows and 13 variables:
#' \describe{
#'   \item{plot}{Plot id}
#'   \item{logit_rock}{logit proportion of cover of plot}
#'   \item{ln_mrvbf}{log of the multiresolution valley bottom flatness index}
#'   \item{ln_prec_yr}{log of annual precipitation}
#'   \item{ln_cv_temp}{log of the coefficient of variation of temperature}
#'   \item{ln_rad_d21}{log of solar radiation at the summer solstice}
#'   \item{sand}{0=soil is not sandy, 1=soil is sandy}
#'   \item{loam}{0=soil is not loam, 1=plot soil is loam}
#'   \item{present}{0=tree absent, 1=tree present}
#'   \item{ln_seed_wt}{log of seed weight in mg }
#'   \item{ln_sla}{log of specific leaf area}
#'   \item{ln_ht}{log of maximum height}
#'   \item{species}{tree species}
#' }
"eucs"