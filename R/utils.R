allpairs <- function(n,
                     init = 1,
                     count = TRUE) {
  aa <- utils::combn(init : n, 2)
  if (count) aa <- rbind(aa, 1 : ncol(aa))
  bb <- apply(aa, 2, list)
  return(lapply(bb, function(x) x[[1]]))
}

is_bbt_model <- function(modout) is.list(modout) && !is.null(modout$model)

get_pwin <- function(modout,
                     selected = NULL,
                     control = NULL) {

  toprob <- function(a, b) a / (a + b)

  aa <- posterior::as_draws_df(modout$model$draws())
  bb <- tibble::tibble(aa) %>% dplyr::select(dplyr::starts_with("beta"))
  colnames(bb) <- modout$wintable$alg_names
  res <- colMeans(bb)
  names <- names(sort(res, decreasing = TRUE))

  if (!is.null(selected)) names <- intersect(names, selected)
  bb <- bb[, names]
  bb <- exp(bb)
  n <- ncol(bb)
  l <- nrow(bb)

  if (is.null(control) || !(control %in% names)) {
    out <- array(NA_real_, dim = c(l, round(n * (n - 1) / 2)))
    out_names <- NULL
    for (xx in allpairs(n)) {
      i <- xx[1]
      j <- xx[2]
      k <- xx[3]
      out[, k] <- toprob(bb[, i, drop = TRUE], bb[, j, drop = TRUE])
      out_names <- c(out_names, paste(names[i], ">", names[j]))
    }
  } else {
    out <- array(NA_real_, dim = c(l, n - 1))
    out_names <- NULL
    i <- which(names == control)
    k <- 1
    for (j in 1 : n)
      if (i < j) {
        out[, k] <- toprob(bb[, i, drop = TRUE], bb[, j, drop = TRUE])
        out_names <- c(out_names, paste(names[i], ">", names[j]))
        k <- k + 1
      } else if (i > j) {
        out[, k] <- toprob(bb[, j, drop = TRUE], bb[, i, drop = TRUE])
        out_names <- c(out_names, paste(names[j], ">", names[i]))
        k <- k + 1
      }
  }

  colnames(out) <- out_names
  return(out)
}


aux_get_win1_rep <- function(modout) {
  mod <- modout$model
  tabx <- modout$wintable
  nn <- tabx$alg_names
  dd <- mod$draws()
  ddim <- dim(dd)
  l <- ddim[1] * ddim[2]
  npairs <- nrow(tabx$table)
  out <- array(dim = c(l, npairs))
  for (i in 1:npairs) out[, i] <- posterior::extract_variable(dd, paste0("win1_rep[", i, "]"))
  coln <- NULL
  for (i in 1 : (length(nn) - 1))
    for (j in (i + 1) : length(nn))
      coln <- c(coln, paste(nn[i], ">", nn[j]))
  colnames(out) <- coln
  return(out)
}


aux_get_ties_rep <- function(modout) {
  mod <- modout$model
  tabx <- modout$wintable
  nn <- tabx$alg_names
  dd <- mod$draws()
  ddim <- dim(dd)
  l <- ddim[1] * ddim[2]
  npairs <- nrow(tabx$table)
  out <- array(dim = c(l, npairs))
  for (i in 1:npairs) out[, i] <- posterior::extract_variable(dd, paste0("tie_rep[", i, "]"))
  coln <- NULL
  for (i in 1 : (length(nn) - 1))
    for (j in (i + 1) : length(nn))
      coln <- c(coln, paste(nn[i], "=", nn[j]))
  colnames(out) <- coln
  return(out)
}


