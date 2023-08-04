#' Generates a BBT model from a table with the results of algorithms on
#'      different data sets
#'
#' @inheritParams make_wintable
#' @param hyper_prior The code of the hyper-prior for the sigma
#'
#'     * 0 = lognormal(0.5) - the default,
#'     * 1 = lognormal(scale),
#'     * 2 = cauchy(scale),
#'     * 3 = normal(scale)
#' @param scale The scale of the hyper prior for the sigma parameter
#' @param ... Extra parameters to be passed to the MCMC algorithm (cmdstanr)
#'
#' @return A BBT model. The model contains:
#'
#'        - `model`: a fit object from cmdstanr
#'        - `wintable`: the wintable generated from the data
#'        - `davidson` whether the Davidson model was used
#'
#' @export
#'
#' @examples
#' \donttest{
#' ss <- ll[1:80, 1:6] # first 20 data sets and 5 algorithms
#'
#' options(mc.cores = parallel::detectCores(logical = FALSE))
#'
#' m1 <- bbtcomp(ss)
#' m2 <- bbtcomp(ss, lrope = FALSE, iter_sampling = 2000, seed = 1234)
#' m3 <- bbtcomp(ss, lrope_value = 0.2, paired = FALSE)
#' m4 <- bbtcomp(ss, deal_with_ties = "davidson", hyper_prior = 3, scale = 3.0)
#'
#' }
bbtcomp <- function(tabx,
                    tabsd = NULL, #data with standard deviations
                    dbcol = 1, ## what column has the db indicator

                    ## build win table parameters
                    lrope = NULL,  ## compute local rope?
                    lrope_value = 0.4,
                    paired = TRUE,
                    deal_with_ties = c("spread", "forget", "davidson",
                                      "add", "random"), ## how to deal with ties

                    ## MCMC options
                    hyper_prior = NULL,
                    scale = NULL,
                    ...  ## for the MCMC
) {

  testit::assert("tabx must be a data frame or an array",
                 is.data.frame(tabx) || is.array(tabx) )
  testit::assert("tabxsd must be a data frame or an array",
                 is.null(tabsd) || is.data.frame(tabsd) || is.array(tabsd))
  testit::assert("dbcol must be a number or a name",
                 is.null(dbcol) || (is.numeric(dbcol) && length(dbcol) == 1) ||
                   (is.character(dbcol) && length(dbcol) == 1) )

  if (length(deal_with_ties) > 1) deal_with_ties <- "spread" #default is a list!

  # hyperprior

  if (is.null(hyper_prior)) hyper_prior <- 0
  if (is.null(scale)) scale <- 0.5


  ## deal with ties

  if (substr(deal_with_ties, 1, 1) == "d")
    use_davidson <- 1
  else
    use_davidson <- 0

  nalg <- ncol(tabx)

  # deal with no dbcol
  if (is.null(dbcol) || dbcol < 1 || dbcol > nalg) dbcol <- 0

  # compute lrope
  if (!is.null(tabsd)) lrope <- TRUE
  if (dbcol == 0 && is.null(tabsd)) lrope <- FALSE
  if (multiple_folds(tabx, dbcol) && is.null(lrope)) lrope <- TRUE


  # create algorithm names
  if (is.null(colnames(tabx))) {
    print("Creating algorithm names")
    if (dbcol == 0) {
      new_names <- paste0("alg", 1 : ncol(tabx))
      colnames(tabx) <- new_names
    } else {
      new_names <- paste0("alg", 1 : (ncol(tabx) - 1))
      colnames(tabx[, -dbcol]) <- new_names
      colnames(tabx[, dbcol]) <- "db"
    }
  }

  # wintable post processing ties
  tablex <- make_wintable(tabx,
                          tabsd,
                          dbcol,
                          lrope,
                          lrope_value,
                          paired,
                          deal_with_ties)



  mod <- mcmcbbt( tablex,
                  hyper_prior,
                  scale,
                  use_davidson,
                   ...)

  return(mod)
}


