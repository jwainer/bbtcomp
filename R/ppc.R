#' Plot the Posterior Predictive Check
#'
#' @param modout A BBT model
#' @param maxshow The maximum number of plots to display
#'
#' @return A ggplot2 plot
#' @export
#'
#' @examples
#' \donttest{
#' m1 = bbtcomp(ll, lrope = FALSE, deal_with_ties = "davidson")
#' ppc(m1)
#' ppc(m1, maxshow = 30)
#' }
ppc <- function(modout,
                maxshow = 20) {

  testit::assert("modout must be a BBT model",
                 is_bbt_model(modout))

  sums <- modout$wintable$table
  y <- sums$win1
  yrep <- aux_get_win1_rep(modout)
  if (modout$davidson) {
    y <- c(y, sums$ties)
    yrep <- cbind(yrep, aux_get_ties_rep(modout))
  }
  nn <- colnames(yrep)
  n <- length(nn)
  if (!is.null(maxshow) && maxshow < n) {
    keep <- sort(sample(n, maxshow))
    y <- y[keep]
    yrep <- yrep[, keep]
    nn <- nn[keep]
  }
  n <- length(nn)
  plt <- bayesplot::bayesplot_grid(
    plots = lapply(1 : n,
                   function(i) bayesplot::ppc_stat(y, yrep, stat = function(a) a[i], binwidth = 2)),
    titles = nn,
    legends = FALSE
  )
  return(plt)
}
