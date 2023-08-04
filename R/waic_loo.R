#' Compute the WAIC for the model
#'
#' @param mod A BBT model
#'
#' @return the waic, as computed by the loo package
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 data sets and 5 algorithms
#'
#' m1 <- bbtcomp(ss)
#' get_waic(m1)
#' }
get_waic <- function(mod) {
  log_lik = posterior::as_draws_matrix(mod$model$draws('log_lik'))
  return(loo::waic(log_lik))
}


#' Compute the loo for the model
#'
#' @param mod A BBT model
#'
#' @return The loo, as computed by the loo package
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 datasets and 5 algorithms
#'
#' m1 <- bbtcomp(ss)
#' get_loo(m1)
#' }
get_loo <- function(mod) return(mod$model$loo())




