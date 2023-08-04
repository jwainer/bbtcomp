#' Summary of the PPC for all parameters
#'
#' @param modout A BBT model
#'
#' @return A data.frame with the different HDI and the proportion of data that
#'         falls within the HDI
#' @export
#'
#' @examples
#' \donttest{
#' m1 = bbtcomp(ll, lrope = FALSE, deal_with_ties = "davidson")
#' m2 = bbtcomp(ll,  deal_with_ties = "spread")
#' table_ppc(m1)
#' table_ppc(m2)
#' }
table_ppc <- function(modout) {

  testit::assert("modout must be a BBT model",
                 is_bbt_model(modout))

  sums <- modout$wintable$table
  y <- sums$win1
  yrep <- aux_get_win1_rep(modout)
  out <- data.frame(hdi = c(0.5, 0.9, 0.95, 1.0), proportion = numeric(4))
  if (modout$davidson) {
    yt <- sums$ties
    ytrep <- aux_get_ties_rep(modout)
    out <- data.frame(hdi = c(0.5, 0.9, 0.95, 1.0),
                      proportion = numeric(4),
                      ties = numeric(4))
  }

  i <- 0
  for (h in c(0.5, 0.9, 0.95)) {
    i <- i + 1
    z <- apply(yrep, 2, HDInterval::hdi, h)
    out[i, 2] <- mean(y <= z[2, ] & y >= z[1, ])
    if (modout$davidson) {
      z <- apply(ytrep, 2, HDInterval::hdi, h)
      out[i, 3] <- mean(yt <= z[2, ] & yt >= z[1, ])
    }
  }
  out[4, 2] <- mean(y <= apply(yrep, 2, max) & y >= apply(yrep, 2, min))
  if (modout$davidson)
    out[4, 3] <- mean(yt <= apply(ytrep, 2, max) & yt >= apply(ytrep, 2, min))
  round(out, 2)
}
