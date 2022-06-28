#' Construct a wintable from a table of fold results for each combination of algorithm and data set
#'
#' @param tabx An array or data frame, where lines are data sets and columns are algorithms.
#' There are three basic use case.
#'     * only \code{tabx}, and there is a column that contains the name or an id for each data set. The dbcol
#'       parameter indicates this column. In this case, there can be multiple lines for the same
#'       data set, where the lines are the results on the different folds of the data set. In this case
#'       the computation of wins and losses can use the local ROPE concept. If the \code{lrope} parameter is true
#'       then local ROPE will be used. Also, in this case, the algorithm can use the paired version
#'       of the local ROPE (if all folds for this data set are the same for all algorithms). If the
#'       \code{paired} parameter is set, the algorithm will use the paired local ROPE
#'     * \code{tabx} and \code{tabsd} are given. In this case, local ROPE can be used (if \code{lrope} is TRUE) but not
#'       the paired version (\code{paired} is silently ignored). If \code{lrope}
#'       is FALSE, \code{tabsd} is silently ignored.
#'     * only \code{tabx} is given and there is no data set column (\code{dbcol} is set to 0). In this case
#'       local ROPE cannot be used (and \code{lrope} and \code{paired} are silently ignored.)
#' @param tabsd A data frame or an array with the standard deviation of the fold measures for each data set and
#'     algorithm. If given,  the parameter \code{lrope} determines whether local ROPE will or will
#'     not be used to compute wins and losses,
#' @param dbcol Column number or name were the data set indicator is located, or 0 for no such column
#' @param lrope Whether to use local ROPE in the pairwise comparisons
#' @param lrope_value The value of the local ROPE threshold (default 0.4)
#' @param paired Whether to use paired local ROPE calculations
#' @param deal_with_ties How to deal with ties in the wintable. The options are:
#'          - spread - half (rounded up) of each ties is added to each algorithm.
#'          - forget: remove the ties column.
#'          - add: add the ties to each algorithm.
#'          - random: add each tie randomly to one of the algorithms.
#'          - davidson: do nothing with the ties. TO be used with the Davidson Bayesian model.
#'
#' @return A wintable object
#' @export
#'
#' @examples
#'\donttest{
#' w1 <- make_wintable(ll)
#'
#' w2 <- make_wintable(as.matrix(ll[, -1]), dbcol = 0)
#'
#' w3 <- make_wintable(ll, lrope_value = 0.3, paired = TRUE)
#'
#' w4 <- make_wintable(ll, deal_with_ties = "random", paired = FALSE)
#'
#'}
#' @importFrom rlang .data
make_wintable <- function(tabx,
                          tabsd = NULL,
                          dbcol=1,
                          lrope=NULL,
                          lrope_value=0.4,
                          paired = TRUE,
                          deal_with_ties = c("spread", "forget", "davidson", "add", "random")
) {

  testit::assert("tabx should have colnames",
                 !is.null(colnames(tabx)) )

  testit::assert("tabx must be a data.frame or an array",
                 is.data.frame(tabx) || is.array(tabx) )

  if (length(deal_with_ties) > 1) deal_with_ties <- "s"

  names <- colnames(tabx)
  if (dbcol != 0) names <- names[-dbcol]

  ## only a mean table given
  if (is.null(tabsd) && !multiple_folds(tabx, dbcol))  {
    if (dbcol != 0) tabx = tabx[,-dbcol]
    out <- compute_differences_no_sd(tabx)
    lrope <- FALSE
    lrope_value <- 0.0

  ## mean and sd  table given
  } else if (!is.null(tabsd)) {

    if (is.null(lrope) || lrope) {
      lrope <- TRUE
      if (dbcol != 0) tabx = tabx[,-dbcol]
      out <- compute_differences_with_sd(tabx, tabsd, lrope, lrope_value)
    } else {
      # lrope == FALSE
      if (dbcol != 0) tabx = tabx[,-dbcol]
      out <- compute_differences_no_sd(tabx)
    }

  ## only tabx with folds given
  } else if ((is.null(paired) || paired) && (is.null(lrope) || lrope)) {
    # paired effect size
    out <- compute_paired_differences(tabx, dbcol, lrope_value)
    lrope <- TRUE

  } else {
    ## not paired,  compute mean and sd
    if (is.null(lrope)) lrope <- TRUE
    dbcol_name <-  colnames(tabx)[dbcol]
    if (is.array(tabx)) tabx <- tibble::as_tibble(tabx)

    meanx <- tabx %>% dplyr::group_by(.data[[dbcol_name]]) %>%
        dplyr::summarise(dplyr::across(tidyselect::vars_select_helpers$where(is.numeric),
                                     ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
        dplyr::select(-{{dbcol_name}})
    sdx <- tabx %>% dplyr::group_by(.data[[dbcol_name]]) %>%
        dplyr::summarise(dplyr::across(tidyselect::vars_select_helpers$where(is.numeric),
                                     ~stats::sd(.x, na.rm = TRUE)), .groups = "drop") %>%
        dplyr::select(-{{dbcol_name}})
    out <- compute_differences_with_sd(meanx, sdx, lrope, lrope_value)
  }

  out2 <- proc_ties(out, deal_with_ties)
  ret <- list(table = out2, alg_names = names,
              lrope = lrope,
              lrope_value = lrope_value,
              table_pre = out,
              ties_proc = deal_with_ties)
  return(ret)
}



compute_differences_with_sd <- function(meanx,
                                sdx,
                                lrope,
                                lrope_value) {

  nalg <- ncol(meanx)
  npairs <- as.integer(nalg * (nalg - 1) / 2)
  out <- data.frame(pi = numeric(npairs), pj = numeric(npairs),
                    win1 = numeric(npairs), win2 =  numeric(npairs),
                    ties = numeric(npairs))

  for (xx in allpairs(nalg)) {
    i <- xx[1]
    j <- xx[2]
    k <- xx[3]
    if (lrope)
      rope <- lrope_value * sqrt((sdx[, i, drop = TRUE]**2 + sdx[, j, drop = TRUE]**2) / 2)
    else
      rope <- 0.0
    delta <- meanx[, i, drop = TRUE] - meanx[, j, drop = TRUE]
    win1 <- sum(as.integer(delta > rope), na.rm = TRUE)
    win2 <- sum(as.integer(delta < -rope), na.rm = TRUE)
    ties <- sum(as.integer(abs(delta) <= rope), na.rm = TRUE)
    out[k, ] <- c(i, j, win1, win2, ties)
  }

  return(out)
}

compute_differences_no_sd <- function(meanx){
  nalg <- ncol(meanx)
  npairs <- as.integer(nalg * (nalg - 1) / 2)
  out <- data.frame(pi = numeric(npairs), pj = numeric(npairs),
                    win1 = numeric(npairs), win2 =  numeric(npairs),
                    ties = numeric(npairs))

  for (xx in allpairs(nalg)) {
    i <- xx[1]
    j <- xx[2]
    k <- xx[3]
    delta <- meanx[, i, drop = TRUE] - meanx[, j, drop = TRUE]
    win1 <- sum(as.integer(delta > 0.0), na.rm = TRUE)
    win2 <- sum(as.integer(delta < 0.0), na.rm = TRUE)
    ties <- sum(as.integer(delta == 0.0), na.rm = TRUE)
    out[k, ] <- c(i, j, win1, win2, ties)
  }
  return(out)
}

compute_paired_differences <- function(tab,
                                       dbcol,
                                       lrope_value) {

  aux_func <- function(x, zzz) {
    nalg <- ncol(x)
    npairs <- as.integer(nalg * (nalg - 1) / 2)
    out <- data.frame(pi = numeric(npairs), pj = numeric(npairs),
                      win1 = numeric(npairs), win2 =  numeric(npairs),
                      ties = numeric(npairs))
    for (xx in allpairs(nalg)) {
      i <- xx[1]
      j <- xx[2]
      k <- xx[3]
      delta <- x[, i, drop = TRUE] - x[, j, drop = TRUE]
      m <- mean(delta, na.rm = TRUE)
      s <- stats::sd(delta, na.rm = TRUE)
      if (is.null(s) || is.na(s))
        rope <- 0.0
      else
        rope <- lrope_value * s
      out[k, ] <- c(i, j,
                    as.integer(m > rope),
                    as.integer(m < -rope),
                    as.integer(m >= -rope && m <= rope))
    }
    return(out)
  }

  names <- colnames(tab)
  dbcol_name <-  names[dbcol]
  out <- tab %>% dplyr::group_by( .data[[dbcol_name]]  ) %>%
    dplyr::group_modify(aux_func) %>%
    dplyr::ungroup() %>%
    dplyr::select(-{{dbcol_name}}) %>%
    dplyr::group_by(.data$pi,.data$pj) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop")
  return(as.data.frame(out))
}



proc_ties <- function(tab, deal_with_ties) {
  aux1 <- function(x) {
    if (x == 0) return(c(0,0))
    a <- sample(1 : x, 1)
    b <- x - a
    return(c(a, b))
  }
  choice <- substr(deal_with_ties, 1, 1)
  if (choice == "s") {
    ties <- ceiling(tab$ties / 2)
    tab$win1 <- tab$win1 + ties
    tab$win2 <- tab$win2 + ties
    tab$ties <- 0
  } else if (choice == "a") {
    ties <- tab$ties
    tab$win1 <- tab$win1 + ties
    tab$win2 <- tab$win2 + ties
    tab$ties <- 0
  } else if (choice == "f") {
    tab$ties <- 0
  } else if (choice == "r") {
    ties <- do.call(rbind, lapply(tab$ties, aux1))
    tab$win1 <- tab$win1 + ties[, 1]
    tab$win2 <- tab$win2 + ties[, 2]
    tab$ties <- 0
  }

  return(tab)
}


multiple_folds <- function(tab, dbcol)  dbcol != 0 &&  length(unique(tab[[dbcol]])) != nrow(tab)
