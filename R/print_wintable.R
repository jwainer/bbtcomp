#' Prints a wintable
#'
#' @param tablex A wintable
#'
#' @return A data.frame with the wintable
#' @export
#'
#' @examples
#' \donttest{
#' w1 <- make_wintable(ll)
#' print_wintable(w1)
#'
#' w2 <- make_wintable(ll, lrope = FALSE)
#' print_wintable(w2)
#'
#' w3 <- make_wintable(ll, lrope_value = 0.3)
#' print_wintable(w3)
#'
#' w4 <- make_wintable(ll, deal_with_ties = "random", paired = FALSE)
#' print_wintable(w4)
#'
#'}
print_wintable <- function(tablex) {
  sums <- tablex$table
  names <- tablex$alg_names
  a <- data.frame(
    alg1 = names[sums$pi],
    alg2 = names[sums$pj],
    win1 = sums$win1,
    win2 = sums$win2)
  if ("ties" %in% colnames(sums)) a <- cbind(a, ties = sums$ties)
  return(a)
}
