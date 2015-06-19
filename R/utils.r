#' @import Matrix
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom magrittr %T>%

# Check if a glmerMod class
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

# magrittr like function to return something else if condition is met
return_if <- function(x, test, y) {
  if (test) y else x
  }
