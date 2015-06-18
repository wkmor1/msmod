#' @import Matrix
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom magrittr %T>%

is_glmerMod <- function(x) {
  x %>%
  base::class(.) %>% 
  magrittr::equals("glmerMod") %>%
  magrittr::not(.) %>%
  {
    if (.) {
      'x must object of class "glmerMod"' %>%
     base::stop(call.=FALSE)
    }
  }
}
