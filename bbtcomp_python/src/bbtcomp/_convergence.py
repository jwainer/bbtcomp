"""Convergence diagnostics for a fitted BBT model."""

from __future__ import annotations

from ._mcmc import BBTModel


def convergence_check(mod: BBTModel) -> str:
    """Run cmdstan's sampler diagnostics on a fitted BBT model.

    Thin wrapper over ``cmdstanpy``'s ``CmdStanMCMC.diagnose()`` (the
    equivalent of the R package's call to ``cmdstan_diagnose`` via
    cmdstanr), printed and also returned as a string.
    """
    report = mod.model.diagnose()
    print(report)
    return report
