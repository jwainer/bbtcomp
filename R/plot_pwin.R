#' Plot the probabilities of an algorithm being better than another.
#'
#' @param modout A BBT model
#' @param selected List of algorithms names to be printed. NULL for all
#' @param control An algorithm name, to be seen as the control. Only the
#'       probabilities of the other algorithms against the control will be
#'       plotted
#' @param rope The low and the high values of ROPE
#' @param hdi The HDI probability to indicate in the plot.
#'
#' @return A ggplot2 plot
#' @export
#'
#' @examples
#' \donttest{
#' m1 = bbtcomp(ll, lrope = FALSE, deal_with_ties = "davidson")
#' plot_pwin(m1, control = "rf")
#' plot_pwin(m1, selected = c("svm","lda","gbm","passive"), rope = c(0.48, 0.52))
#' plot_pwin(m1,  hdi = 0.95)
#' }
plot_pwin <- function(modout,
                    selected = NULL,
                    control = NULL,
                    rope = c(0.45, 0.55),
                    hdi = 0.89) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop(
      "Package \"bayesplot\" must be installed to use this function.",
      call. = FALSE
    )
  }
  testit::assert("modout must be a BBT model",
                 is_bbt_model(modout))

  testit::assert("selected  must be a list of alg names",
                 is.null(selected) || is.character(selected) )

  testit::assert("control must be an alg name",
                 is.null(control) || ( is.character(control) && length(control) == 1) )

  testit::assert("rope  must be pair of probabilities",
                 is.null(rope) ||
                   (is.numeric(rope) && all(rope <= 1.0) && all(rope >= 0.0) &&
                      length(rope) == 2) )


  aa <- get_pwin(modout, selected)
  nn <- ncol(aa)
  stp <- Filter(function(x) x < nn, 0.5 + cumsum(1 : nn))
  g1 <- bayesplot::mcmc_intervals(aa, point_size = 2, prob_outer = 1, prob = hdi) +
    ggplot2::geom_vline(xintercept = 0.5, color = "black") +
    ggplot2::geom_hline(yintercept = stp, alpha = 0.2)
  if (!is.null(rope)) g1 <- g1 + ggplot2::geom_vline(xintercept = rope, color = "gray")
  return(g1)
}


