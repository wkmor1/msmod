#' @import ggplot2
#' @import lme4
#' @import R2jags
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%

# to use magrittr shortcut
utils::globalVariables(".")

# Check if a glmerMod class
is_glmerMod <- function(x) {
  x %>%
  base::class(.) %>% 
  magrittr::equals("glmerMod") %>%
  magrittr::not(.) %>%
  {
    if (.) {
      "x must be object of class 'glmerMod'" %>%
     base::stop(call. = FALSE)
    }
  }
}

# magrittr like functions to return something else if condition is not met
return_if_not <- function(x, test, y) {
  if (test) y else x
}

# bind y to x if y not in z
bind_if_not_in <- function(x, y, z, out=base::get(z, parent.frame())) {
  x %>%
  return_if_not(
    y %>%
    magrittr::extract2(z) %>%
    base::is.null(.),
    x %>%
    magrittr::inset2(z, out)
  )
}

# evaluate a function with a list of named arguments (do.call doesnt keep 
# the names.)
eval_with_args <- function(args, fun) {
  fun <-
    base::substitute(fun) %>%
    base::deparse(.)
  
  base::names(args) %>%
  base::paste0(., "=") %>%
  base::paste0(names(args), collapse = ",") %>%
  base::sprintf(
    "with(args, %s(%s))",
    fun,
    .) %>%
  base::parse(text = .) %>%
  base::eval(.)
}
