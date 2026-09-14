#' Print a summary of the probabilities of one algorithm being better than another
#'
#' @param modout A BBT model
#' @param selected List of algorithm names to be printed. NULL for all
#' @param control An algorithm name, to be seen as the control. Only the
#'       probabilities of the other algorithms against the control will be
#'       printed
#' @param short Print only the short summary (mean, delta, above.50 and in.rope)
#' @param rope The low and the high values of ROPE
#' @param columns Which columns to print in the summary
#' @param hdi HDI for the low/high and delta columns in the summary
#' @param ndigits Number of digits to print the probabilities
#'
#' @return A data.frame (S3 class \code{"bbt_pwin_table"}) with \code{larger}
#'   and \code{smaller} columns identifying each compared pair of algorithms
#'   (\code{larger} is the one with the higher estimated ability), plus the
#'   requested summary columns. When printed, \code{larger}/\code{smaller} are
#'   shown combined as a single readable \code{"A > B"} column, but the two
#'   fields stay separately accessible for programmatic use, e.g.
#'   \code{tab$larger} or \code{subset(tab, larger == "svm")}.
#' @export
#'
#' @examples
#' \donttest{
#' m1 = bbtcomp(ll)
#' table_pwin(m1, control = "rf", columns = c("median","delta","in.rope"))
#' table_pwin(m1, selected = c("svm","lda","gbm","passive"), rope = c(0.48, 0.52))
#' table_pwin(m1, short = FALSE, hdi = 0.95)
#' table_pwin(m1, columns = c("mean","low", "high"), hdi = 0.90)
#' }
table_pwin <- function(
    modout,
    selected = NULL,
    control = NULL,
    short = TRUE,
    rope = c(0.45, 0.55),
    columns = c("median", "mean", "low", "high", "delta",  "above.50", "in.rope"),
    hdi = 0.89,
    ndigits = 2) {

  testit::assert("modout must be a BBT model",
                 is_bbt_model(modout))

  testit::assert("selected  must be a list of alg names",
                 is.null(selected) || is.character(selected) )

  testit::assert("control must be an alg name",
                 is.null(control) || (is.character(control) && length(control) == 1) )

  if (short && length(columns) == 7) columns <- c("mean","delta", "above.50", "in.rope")
  columns <- intersect( columns, c("median", "mean", "low", "high", "delta", "above.50", "in.rope"))
  zz <- get_pwin(modout, selected, control)
  larger <- attr(zz, "larger")
  smaller <- attr(zz, "smaller")
  n <- length(larger)

  df <- data.frame(larger = larger,
                   smaller = smaller,
                   median = numeric(n),
                   mean = numeric(n),
                   low = numeric(n),
                   high = numeric(n),
                   delta = numeric(n),
                   above.50 = numeric(n),
                   in.rope = numeric(n),
                   stringsAsFactors = FALSE)

  for (i in 1 : n) {
    aux <- zz[, i]
    hdaux <- HDInterval::hdi(aux, hdi)
    a1 <- stats::median(aux, na.rm = T)
    a2 <- as.numeric(hdaux[1])
    a3 <- as.numeric(hdaux[2])
    a4 <- as.numeric(mean(aux))
    a5 <- as.numeric(mean(aux > 0.5))
    a6 <- as.numeric(mean(aux <= rope[2] & aux >= rope[1]))
    df[i, "median"] <- round(a1, ndigits)
    df[i, "low"] <- round(a2, ndigits)
    df[i, "high"] <- round(a3, ndigits)
    df[i, "delta"] <- round(a3 - a2, ndigits)
    df[i, "mean"] <- round(a4, ndigits)
    df[i, "above.50"] <- round(a5, ndigits)
    df[i, "in.rope"] <- round(a6, ndigits)
  }
  result <- df[, c("larger", "smaller", columns)]
  class(result) <- c("bbt_pwin_table", class(result))
  return(result)
}
