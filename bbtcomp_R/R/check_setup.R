#' Check whether CmdStan is installed and discoverable by cmdstanr
#'
#' Installing \code{bbtcomp} (via CRAN, \code{remotes::install_github()},
#' etc.) installs the \code{cmdstanr} R package, but it cannot install or
#' verify CmdStan itself (the compiled Stan backend) -- CmdStan is not an
#' R package, so R's package installers have no way to know about it, let
#' alone install it. Run this once after installing bbtcomp (before
#' fitting a model) to get a clear message if CmdStan is missing, instead
#' of discovering it later as a cryptic error from deep inside cmdstanr.
#'
#' @param raise_on_error Whether to \code{stop()} with setup instructions
#'   if CmdStan can't be found, instead of just returning \code{FALSE}
#'
#' @return \code{TRUE} (invisibly) if CmdStan was found (its path is
#'   printed); \code{FALSE} if it wasn't (a message with setup
#'   instructions is printed), unless \code{raise_on_error} is
#'   \code{TRUE}, in which case this stops with those instructions
#'   instead.
#' @export
#'
#' @examples
#' \donttest{
#' check_setup()
#' }
check_setup <- function(raise_on_error = FALSE) {
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)

  if (is.null(path)) {
    msg <- paste0(
      "CmdStan could not be found.\n\n",
      "bbtcomp needs CmdStan (the Stan C++ backend) in addition to the ",
      "cmdstanr R package -- installing bbtcomp installs cmdstanr but ",
      "cannot install or check for CmdStan itself, since it isn't an R ",
      "package.\n\n",
      "To fix this, do ONE of the following:\n\n",
      "  1. If you don't have CmdStan yet, install it with:\n",
      "         cmdstanr::install_cmdstan()\n",
      "     (this downloads and compiles CmdStan; needs a C++ compiler and\n",
      "     takes several minutes). See\n",
      "     https://mc-stan.org/cmdstanr/articles/cmdstanr.html for details.\n\n",
      "  2. If you already have a CmdStan installation (for example, one\n",
      "     installed for the Python package via\n",
      "     cmdstanpy.install_cmdstan()), point cmdstanr at it instead of\n",
      "     installing a second copy:\n",
      "         cmdstanr::set_cmdstan_path(\"/path/to/cmdstan-x.y.z\")\n"
    )
    if (raise_on_error) stop(msg, call. = FALSE)
    cat(msg)
    return(invisible(FALSE))
  }

  cat("CmdStan found at:", path, "\n")
  return(invisible(TRUE))
}


#' @keywords internal
.ensure_cmdstan <- function() {
  ok <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  if (!ok) check_setup(raise_on_error = TRUE)
}
