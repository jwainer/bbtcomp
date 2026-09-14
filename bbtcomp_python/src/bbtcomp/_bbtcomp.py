"""The main, end-to-end entry point: table of results -> fitted BBT model."""

from __future__ import annotations

from typing import List, Optional, Union

import pandas as pd

from ._mcmc import BBTModel, mcmcbbt
from ._wintable import make_wintable


def bbtcomp(tabx: pd.DataFrame,
           tabsd: Optional[pd.DataFrame] = None,
           dbcol: Union[int, str, None] = 0,
           lrope: Optional[bool] = None,
           lrope_value: float = 0.4,
           paired: bool = True,
           deal_with_ties: Union[str, List[str]] = "spread",
           hyper_prior: int = 0,
           scale: float = 0.5,
           output_dir: Optional[str] = None,
           **kwargs) -> BBTModel:
    """Fit a Bayesian Bradley-Terry model comparing algorithms across data sets.

    This builds the win/loss/tie table with :func:`make_wintable` and then
    fits the model with :func:`mcmcbbt`, in one call.

    Parameters
    ----------
    tabx : pd.DataFrame
        Rows are data sets, columns are algorithms (higher is better; for
        error-like metrics, negate the values first). One column may
        instead be a data-set identifier -- see ``dbcol``.
    tabsd : pd.DataFrame, optional
        Standard deviation of the fold measures, aligned with ``tabx``.
    dbcol : int, str or None
        Column position (default 0) or name identifying the data set,
        or ``None``/``-1`` if there is no such column.
    lrope : bool, optional
        Whether to use the local ROPE. Auto-detected from the data shape
        when left as ``None``.
    lrope_value : float
        The local ROPE threshold (default 0.4).
    paired : bool
        Whether to use the paired version of the local ROPE.
    deal_with_ties : str
        "spread" (default), "forget", "davidson", "add" or "random".
    hyper_prior : int
        0 = lognormal(0, 0.5) [default], 1 = lognormal(0, scale),
        2 = cauchy(0, scale), 3 = normal(0, scale).
    scale : float
        Scale of the hyper-prior for ``sigma``.
    output_dir : str, optional
        Directory for cmdstanpy's compiled model / sampler csv output.
    **kwargs
        Extra arguments forwarded to ``cmdstanpy.CmdStanModel.sample``.

    Returns
    -------
    BBTModel

    Examples
    --------
    >>> from bbtcomp import bbtcomp, load_ll
    >>> ll = load_ll()
    >>> ss = ll.iloc[0:80, 0:6]
    >>> m1 = bbtcomp(ss)                                   # doctest: +SKIP
    >>> m2 = bbtcomp(ss, lrope=False, iter_sampling=2000)   # doctest: +SKIP
    """
    use_davidson = str(deal_with_ties).lower().startswith("d")

    wintable = make_wintable(
        tabx,
        tabsd=tabsd,
        dbcol=dbcol,
        lrope=lrope,
        lrope_value=lrope_value,
        paired=paired,
        deal_with_ties=deal_with_ties,
    )

    return mcmcbbt(
        wintable,
        hyper_prior=hyper_prior,
        scale=scale,
        use_davidson=use_davidson,
        output_dir=output_dir,
        **kwargs,
    )
