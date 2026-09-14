"""WAIC and PSIS-LOO model-fit criteria, via arviz."""

from __future__ import annotations

from ._mcmc import BBTModel


def _idata(mod: BBTModel):
    import arviz as az
    return az.from_cmdstanpy(posterior=mod.model, log_likelihood="log_lik")


def get_waic(mod: BBTModel):
    """Compute the WAIC for a fitted BBT model (via ``arviz.waic``)."""
    import arviz as az
    return az.waic(_idata(mod))


def get_loo(mod: BBTModel):
    """Compute the PSIS-LOO for a fitted BBT model (via ``arviz.loo``)."""
    import arviz as az
    return az.loo(_idata(mod))
