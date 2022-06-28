#' Compute the WAIC for the model
#'
#' @param mod A BBT model (use_log_lik must be set)
#'
#' @return the waic, as computed by the loo package
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 data sets and 5 algorithms
#'
#' m1 <- bbtcomp(ss, use_log_lik = TRUE)
#' get_waic(m1)
#' }
get_waic <- function(mod) {
  if (mod$log_lik) {
    log_lik = posterior::as_draws_matrix(mod$model$draws('log_lik'))
    return(loo::waic(log_lik))
  }
  stop("use_log_lik was not set for this model")
}


#' Compute the loo for the model
#'
#' @param mod A BBT model (use_log_lik must be set)
#'
#' @return The loo, as computed by the loo package
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 datasets and 5 algorithms
#'
#' m1 <- bbtcomp(ss, use_log_lik = TRUE)
#' get_loo(m1)
#' }
get_loo <- function(mod) {
  if (mod$log_lik) return(mod$model$loo())
  stop("use_log_lik was not set for this model")
}


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
