#' Print convergence diagnostics
#'
#' @param mod A BBT model
#'
#' @return The output from cmdstand_diagnose() from the cmdstanr package
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 data sets and 5 algorithms
#'
#' m1 <- bbtcomp(ss)
#' convergence_check(m1)
#' }
convergence_check <- function(mod){
  mod$model$cmdstan_diagnose()
}
