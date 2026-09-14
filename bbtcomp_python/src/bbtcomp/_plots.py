"""Plots for a fitted BBT model (matplotlib-based equivalents of the R
package's ``bayesplot``-based ``plot_pwin`` / ``plot_ppc``)."""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np

from ._mcmc import BBTModel, aux_get_ties_rep, aux_get_win1_rep, get_pwin
from ._utils import hdi, is_bbt_model


def plot_pwin(modout: BBTModel,
             selected: Optional[List[str]] = None,
             control: Optional[str] = None,
             rope: Optional[Sequence[float]] = (0.45, 0.55),
             hdi_prob: float = 0.89):
    """Forest plot of P(alg_i beats alg_j) for every requested pair.

    Mirrors the R package's ``plot_pwin`` (built on
    ``bayesplot::mcmc_intervals``): one horizontal interval per pair,
    a vertical line at 0.5, and (optionally) a shaded ROPE band.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    assert is_bbt_model(modout), "modout must be a BBT model"
    assert selected is None or all(isinstance(s, str) for s in selected)
    assert control is None or isinstance(control, str)
    assert rope is None or (len(rope) == 2 and 0.0 <= rope[0] <= 1.0 and 0.0 <= rope[1] <= 1.0)

    zz, names, _larger, _smaller = get_pwin(modout, selected, control)
    n = len(names)

    means = zz.mean(axis=0)
    bounds = np.array([hdi(zz[:, i], hdi_prob) for i in range(n)])
    # innermost interval (50%) drawn thicker, like bayesplot's default
    inner = np.array([hdi(zz[:, i], 0.5) for i in range(n)])

    fig, ax = plt.subplots(figsize=(6, 0.4 * n + 1))
    y = np.arange(n, 0, -1)

    ax.hlines(y, bounds[:, 0], bounds[:, 1], color="steelblue", linewidth=1.5)
    ax.hlines(y, inner[:, 0], inner[:, 1], color="steelblue", linewidth=4)
    ax.scatter(means, y, color="black", zorder=3, s=20)

    ax.axvline(0.5, color="black", linewidth=1)
    if rope is not None:
        ax.axvspan(rope[0], rope[1], color="gray", alpha=0.2)

    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlabel("P(row algorithm beats column algorithm)")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    return fig


def plot_ppc(modout: BBTModel, maxshow: int = 20):
    """Posterior predictive check plots: observed vs. replicated win counts.

    Mirrors the R package's ``plot_ppc`` (built on
    ``bayesplot::ppc_stat``): a small histogram per pair, of the
    posterior-predictive distribution of ``win1`` (and ``ties`` when the
    Davidson model is used), with a vertical line at the observed value.

    Parameters
    ----------
    modout : BBTModel
    maxshow : int, optional
        Cap on the number of pairs plotted (a random subset is chosen
        when there are more).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    assert is_bbt_model(modout), "modout must be a BBT model"

    sums = modout.wintable.table
    y = sums["win1"].to_numpy()
    yrep = aux_get_win1_rep(modout)
    names = list(yrep.columns)

    if modout.davidson:
        y = np.concatenate([y, sums["ties"].to_numpy()])
        ties_rep = aux_get_ties_rep(modout)
        yrep_arr = np.concatenate([yrep.to_numpy(), ties_rep.to_numpy()], axis=1)
        names = names + list(ties_rep.columns)
    else:
        yrep_arr = yrep.to_numpy()

    n = len(names)
    if maxshow is not None and maxshow < n:
        rng = np.random.default_rng()
        keep = np.sort(rng.choice(n, size=maxshow, replace=False))
        y = y[keep]
        yrep_arr = yrep_arr[:, keep]
        names = [names[i] for i in keep]
        n = len(names)

    ncols = min(4, n) or 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.2 * nrows), squeeze=False)

    for idx in range(nrows * ncols):
        ax = axes[idx // ncols][idx % ncols]
        if idx >= n:
            ax.axis("off")
            continue
        ax.hist(yrep_arr[:, idx], bins=20, color="steelblue", alpha=0.7)
        ax.axvline(y[idx], color="black", linewidth=1.5)
        ax.set_title(names[idx], fontsize=9)
        ax.set_yticks([])

    fig.tight_layout()
    return fig
