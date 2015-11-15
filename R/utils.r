#' @importClassesFrom lme4 glmerMod
#' @importFrom R2jags jags
#' @importFrom rstan stan
#' @importFrom dplyr %>%

# to use magrittr shortcut
utils::globalVariables(".")

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

# evaluate a function with a list of named arguments (do.call doesn't keep
# the names.)
eval_with_args <- function(args, fun) {
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
