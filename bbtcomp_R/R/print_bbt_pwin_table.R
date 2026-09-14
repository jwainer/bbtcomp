#' Print method for \code{table_pwin()} results
#'
#' \code{\link{table_pwin}} stores the two algorithm names of each compared
#' pair as separate \code{larger}/\code{smaller} columns (\code{larger} being
#' the one with the higher estimated ability), so they can be used
#' programmatically (e.g. \code{tab$larger}, \code{subset(tab, larger ==
#' "svm")}) without having to parse a combined string back apart. For
#' display, this print method combines them back into a single readable
#' \code{"A > B"} \code{"pair"} column, matching how the table used to look
#' when printed.
#'
#' @param x A \code{bbt_pwin_table}, as returned by \code{\link{table_pwin}}
#' @param ... Passed on to \code{print.data.frame}
#'
#' @return \code{x}, invisibly
#' @export
#' @noRd
print.bbt_pwin_table <- function(x, ...) {
  disp <- x
  class(disp) <- "data.frame"
  disp$pair <- paste(disp$larger, ">", disp$smaller)
  disp$larger <- NULL
  disp$smaller <- NULL
  disp <- disp[, c("pair", setdiff(colnames(disp), "pair")), drop = FALSE]
  print.data.frame(disp, ...)
  invisible(x)
}
