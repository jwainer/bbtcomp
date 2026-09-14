"""Tabular summaries of a fitted BBT model."""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

from ._mcmc import BBTModel, aux_get_ties_rep, aux_get_win1_rep, get_pwin
from ._pwin_table import PwinTable
from ._utils import hdi, is_bbt_model
from ._wintable import Wintable

_ALL_COLUMNS = ("median", "mean", "low", "high", "delta", "above.50", "in.rope")


def table_pwin(modout: BBTModel,
              selected: Optional[List[str]] = None,
              control: Optional[str] = None,
              short: bool = True,
              rope: Sequence[float] = (0.45, 0.55),
              columns: Sequence[str] = _ALL_COLUMNS,
              hdi_prob: float = 0.89,
              ndigits: int = 2) -> pd.DataFrame:
    """Summarize P(alg_i beats alg_j) for every pair of algorithms.

    Parameters
    ----------
    modout : BBTModel
    selected : list of str, optional
        Restrict the comparison to these algorithms.
    control : str, optional
        If given, only comparisons against this algorithm are shown.
    short : bool
        If True (default) and ``columns`` was left at its default, only
        "mean", "delta", "above.50" and "in.rope" are shown.
    rope : (float, float)
        Region of practical equivalence for the "in.rope" column.
    columns : sequence of str
        Which summary columns to include.
    hdi_prob : float
        HDI mass used for the "low"/"high"/"delta" columns.
    ndigits : int
        Rounding precision.

    Returns
    -------
    PwinTable
        A ``pandas.DataFrame`` subclass with ``larger`` and ``smaller``
        columns identifying each compared pair of algorithms (``larger``
        is the one with the higher estimated ability), plus the requested
        summary columns. When printed (or shown in Jupyter),
        ``larger``/``smaller`` are displayed combined as a single readable
        ``"pair"`` column (e.g. ``"A > B"``), but the two fields stay
        separately accessible for programmatic use, e.g. ``tp["larger"]``
        or ``tp[tp["larger"] == "svm"]`` -- it's an ordinary DataFrame
        otherwise, just with a different display.
    """
    assert is_bbt_model(modout), "modout must be a BBT model"
    assert selected is None or all(isinstance(s, str) for s in selected)
    assert control is None or isinstance(control, str)

    columns = list(columns)
    if short and len(columns) == 7:
        columns = ["mean", "delta", "above.50", "in.rope"]
    columns = [c for c in columns if c in _ALL_COLUMNS]

    zz, names, larger, smaller = get_pwin(modout, selected, control)
    n = len(names)

    out = {"larger": larger, "smaller": smaller}
    medians = np.zeros(n)
    means = np.zeros(n)
    lows = np.zeros(n)
    highs = np.zeros(n)
    deltas = np.zeros(n)
    above50 = np.zeros(n)
    inrope = np.zeros(n)

    for i in range(n):
        col = zz[:, i]
        low, high = hdi(col, hdi_prob)
        medians[i] = np.median(col)
        means[i] = np.mean(col)
        lows[i] = low
        highs[i] = high
        deltas[i] = high - low
        above50[i] = np.mean(col > 0.5)
        inrope[i] = np.mean((col >= rope[0]) & (col <= rope[1]))

    out["median"] = np.round(medians, ndigits)
    out["mean"] = np.round(means, ndigits)
    out["low"] = np.round(lows, ndigits)
    out["high"] = np.round(highs, ndigits)
    out["delta"] = np.round(deltas, ndigits)
    out["above.50"] = np.round(above50, ndigits)
    out["in.rope"] = np.round(inrope, ndigits)

    df = PwinTable(out)
    return df[["larger", "smaller"] + columns]


def table_ppc(modout: BBTModel) -> pd.DataFrame:
    """Posterior predictive check summary: HDI coverage of ``win1`` (and
    ``ties`` when the Davidson model was used).

    Returns
    -------
    pd.DataFrame with columns "hdi", "proportion" (and "ties" if Davidson).
    """
    assert is_bbt_model(modout), "modout must be a BBT model"

    sums = modout.wintable.table
    y = sums["win1"].to_numpy()
    yrep = aux_get_win1_rep(modout).to_numpy()

    hdis = [0.5, 0.9, 0.95]
    proportion = []
    for h in hdis:
        bounds = np.array([hdi(yrep[:, k], h) for k in range(yrep.shape[1])])
        proportion.append(float(np.mean((y >= bounds[:, 0]) & (y <= bounds[:, 1]))))
    proportion.append(float(np.mean((y >= yrep.min(axis=0)) & (y <= yrep.max(axis=0)))))

    out = {"hdi": [0.5, 0.9, 0.95, 1.0], "proportion": np.round(proportion, 2)}

    if modout.davidson:
        yt = sums["ties"].to_numpy()
        ytrep = aux_get_ties_rep(modout).to_numpy()
        tie_proportion = []
        for h in hdis:
            bounds = np.array([hdi(ytrep[:, k], h) for k in range(ytrep.shape[1])])
            tie_proportion.append(float(np.mean((yt >= bounds[:, 0]) & (yt <= bounds[:, 1]))))
        tie_proportion.append(float(np.mean((yt >= ytrep.min(axis=0)) & (yt <= ytrep.max(axis=0)))))
        out["ties"] = np.round(tie_proportion, 2)

    return pd.DataFrame(out)


def table_wintable(tablex: Wintable, which: str = "pos") -> pd.DataFrame:
    """Render a :class:`Wintable` as a readable data frame.

    Parameters
    ----------
    tablex : Wintable
    which : str
        "pos" (default) for the table after tie-processing, "pre" for
        before, or "both" for both side by side.
    """
    if which == "both":
        t1 = table_wintable(tablex, which="pos")
        t2 = table_wintable(tablex, which="pre")
        return pd.concat([t1, t2.iloc[:, 2:]], axis=1)

    sums = tablex.table if which == "pos" else tablex.table_pre
    names = np.array(tablex.alg_names)
    out = pd.DataFrame({
        "alg1": names[sums["pi"].to_numpy()],
        "alg2": names[sums["pj"].to_numpy()],
        "win1": sums["win1"].to_numpy(),
        "win2": sums["win2"].to_numpy(),
    })
    if "ties" in sums.columns and not (sums["ties"] == 0).all():
        out["ties"] = sums["ties"].to_numpy()

    if which == "pre":
        out.columns = list(out.columns[:2]) + [f"{c}(pre)" for c in out.columns[2:]]
    return out
