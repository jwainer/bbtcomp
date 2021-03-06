% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcbbt.R
\name{mcmcbbt}
\alias{mcmcbbt}
\title{Run the MCMC on a BBT model}
\usage{
mcmcbbt(
  tablex,
  simplified = TRUE,
  hyper_prior = 0,
  scale = 0.5,
  use_davidson = FALSE,
  use_log_lik = FALSE,
  dir = NULL,
  ...
)
}
\arguments{
\item{tablex}{A wintable}

\item{simplified}{Whether to use the simplified model. If hyper_prior, or
scale are set, or davidson is set to TRUE, then simplified should be FALSE}

\item{hyper_prior}{The code of the hyper-prior for the sigma
* 0 = lognormal(0.5) - the default,
* 1 = lognormal(scale),
* 2 = cauchy(scale),
* 3 = normal(scale)}

\item{scale}{The scale of the hyper prior for the sigma parameter}

\item{use_davidson}{Whether to use the Davidson model to deal with ties}

\item{use_log_lik}{Whether to compute the log likelihood for each data
in order to compute the WAIC and loo of the model}

\item{dir}{The direcotry to store the cmdstanr files (compiled model and samples)}

\item{...}{Further parameters to the MCMC (cmdstanr model and fit)}
}
\value{
A fitted BBT model
}
\description{
Run the MCMC on a BBT model
}
\examples{
\donttest{
ss <- ll[1:80, 1:6] # first 20 datasets and 5 algorithms

wintable <- make_wintable(ss, lrope = FALSE)

m1 <- mcmcbbt(wintable, TRUE, hyper_prior = 0, scale = 0.5,
     use_davidson = TRUE, use_log_lik = FALSE)

m2 <- mcmcbbt(wintable, TRUE, iter_sampling = 2000, seed = 1234)
}
}
