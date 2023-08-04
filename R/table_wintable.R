#' Prints a wintable
#'
#' @param tablex A wintable
#' @param which Whether to print the wintable after the processing of ties ("pos")
#'  or before the processing ("pre") or both  ("both")
#'
#' @return A data.frame with the wintable
#' @export
#'
#' @examples
#' \donttest{
#' w1 <- make_wintable(ll)
#' table_wintable(w1)
#'
#' w2 <- make_wintable(ll, lrope = FALSE)
#' table_wintable(w2)
#'
#' w3 <- make_wintable(ll, lrope_value = 0.3)
#' table_wintable(w3)
#'
#' w4 <- make_wintable(ll, deal_with_ties = "random", paired = FALSE)
#' table_wintable(w4)
#'
#'}
table_wintable <- function(tablex, which="pos") {
  if (which=="both") {
    t1 <- table_wintable(tablex,which="pos")
    t2 <- table_wintable(tablex,which="pre")
    tt = cbind(t1, t2[,-c(1,2)])
    return(tt)
  }
  if (which=="pos")  sums <- tablex$table else  sums <- tablex$table_pre
  names <- tablex$alg_names
  a <- data.frame(
    alg1 = names[sums$pi],
    alg2 = names[sums$pj],
    win1 = sums$win1,
    win2 = sums$win2)
  if ("ties" %in% colnames(sums) & !all(sums$ties==0)) a <- cbind(a, ties = sums$ties)
  if (which=="pre") colnames(a) = c(colnames(a)[1:2], paste0(colnames(a)[-c(1,2)],"(pre)"))
  return(a)

}
