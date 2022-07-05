#' Print a summary of the probabilities of an algorithms being better than another
#'
#' @param modout A BBT model
#' @param selected List of algorithms names to be printed. NULL for all
#' @param control An algorithm name, to be seen as the control. Only the
#'       probabilities of the other algorithms against the control will be
#'       printed
#' @param short Print only the short summary (above.50 and in.rope)
#' @param rope Thelow and the high values of ROPE
#' @param columns Which columns to print in the summary
#' @param hdi HDI for the low/high  and delta  columns in the summary
#' @param ndigits Number of digits to print the probabilties
#'
#' @return A data.frame
#' @export
#'
#' @examples
#' \donttest{
#' m1 = bbtcomp(ll)
#' summary_pwin(m1, control = "rf", columns = c("median","delta","in.rope"))
#' summary_pwin(m1, selected = c("svm","lda","gbm","passive"), rope = c(0.48, 0.52))
#' summary_pwin(m1, short = FALSE, hdi = 0.95)
#' }
summary_pwin <- function(
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

  if (short && length(columns) == 7) columns <- c("mean","delta", "in.rope")
  columns <- intersect( columns, c("median", "mean", "low", "high", "delta", "above.50", "in.rope"))
  zz <- get_pwin(modout, selected, control)
  nn <- colnames(zz)
  n <- length(nn)

  df <- data.frame(pair = nn,
                   median = numeric(n),
                   mean = numeric(n),
                   low = numeric(n),
                   high = numeric(n),
                   delta = numeric(n),
                   above.50 = numeric(n),
                   in.rope = numeric(n))

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
  return(df[, c("pair", columns)])
}
