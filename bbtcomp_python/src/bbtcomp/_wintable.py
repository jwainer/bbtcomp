"""Construct a wintable from raw per-fold algorithm results.

This mirrors ``R/make_wintable.R`` in the R package: given a table of
algorithm measures across data sets (optionally with several rows per
data set for cross-validation folds), compute, for every pair of
algorithms, how many times each one "won", "lost" or "tied" -- using the
local ROPE (region of practical equivalence) concept when fold-level
variability is available.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from ._utils import allpairs

_TIE_STRATEGIES = {"spread": "s", "add": "a", "forget": "f",
                   "random": "r", "davidson": "d"}


def _normalize_ties(deal_with_ties: Union[str, List[str]]) -> str:
    if isinstance(deal_with_ties, (list, tuple)):
        deal_with_ties = deal_with_ties[0]
    deal_with_ties = deal_with_ties.lower()
    if deal_with_ties in _TIE_STRATEGIES:
        return _TIE_STRATEGIES[deal_with_ties]
    if deal_with_ties in _TIE_STRATEGIES.values():
        return deal_with_ties
    raise ValueError(
        f"deal_with_ties must be one of {sorted(_TIE_STRATEGIES)} "
        f"(or their first letters), got {deal_with_ties!r}"
    )


@dataclass
class Wintable:
    """The result of :func:`make_wintable`.

    Attributes
    ----------
    table : pd.DataFrame
        Columns ``pi, pj, win1, win2, ties`` (0-indexed algorithm
        positions into ``alg_names``) after ties have been processed
        according to ``ties_proc``.
    alg_names : list of str
        Names of the algorithms being compared, in ``tabx`` column order.
    lrope : bool
        Whether the local ROPE was used to compute wins/losses/ties.
    lrope_value : float
        The local ROPE threshold used.
    table_pre : pd.DataFrame
        Same as ``table`` but *before* ties were processed.
    ties_proc : str
        The tie-handling strategy that was applied ("s", "a", "f", "r"
        or "d").
    paired : bool
        Whether the paired version of the local ROPE was used.
    """

    table: pd.DataFrame
    alg_names: List[str]
    lrope: bool
    lrope_value: float
    table_pre: pd.DataFrame
    ties_proc: str
    paired: bool = True


def _resolve_dbcol(columns: List[str], dbcol) -> Optional[str]:
    if dbcol is None or dbcol == -1:
        return None
    if isinstance(dbcol, str):
        return dbcol
    return columns[dbcol]


def multiple_folds(tabx: pd.DataFrame, dbcol_name: Optional[str]) -> bool:
    """True if the data set column has repeated values (i.e. folds)."""
    if dbcol_name is None:
        return False
    return tabx[dbcol_name].nunique() != len(tabx)


def _empty_out(npairs: int) -> pd.DataFrame:
    return pd.DataFrame({
        "pi": np.zeros(npairs, dtype=int),
        "pj": np.zeros(npairs, dtype=int),
        "win1": np.zeros(npairs, dtype=int),
        "win2": np.zeros(npairs, dtype=int),
        "ties": np.zeros(npairs, dtype=int),
    })


def compute_differences_no_sd(meanx: pd.DataFrame) -> pd.DataFrame:
    """Wins/losses/ties from means only (no ROPE, an exact tie is a tie)."""
    nalg = meanx.shape[1]
    npairs = nalg * (nalg - 1) // 2
    out = _empty_out(npairs)
    for i, j, k in allpairs(nalg):
        delta = meanx.iloc[:, i] - meanx.iloc[:, j]
        out.loc[k, "pi"] = i
        out.loc[k, "pj"] = j
        out.loc[k, "win1"] = int(np.nansum(delta > 0.0))
        out.loc[k, "win2"] = int(np.nansum(delta < 0.0))
        out.loc[k, "ties"] = int(np.nansum(delta == 0.0))
    return out


def compute_differences_with_sd(meanx: pd.DataFrame, sdx: Optional[pd.DataFrame],
                                lrope: bool, lrope_value: float) -> pd.DataFrame:
    """Wins/losses/ties using a (non-paired) local ROPE from mean+sd tables."""
    nalg = meanx.shape[1]
    npairs = nalg * (nalg - 1) // 2
    out = _empty_out(npairs)
    for i, j, k in allpairs(nalg):
        if lrope and sdx is not None:
            rope = lrope_value * np.sqrt((sdx.iloc[:, i] ** 2 + sdx.iloc[:, j] ** 2) / 2)
        else:
            rope = 0.0
        delta = meanx.iloc[:, i] - meanx.iloc[:, j]
        out.loc[k, "pi"] = i
        out.loc[k, "pj"] = j
        out.loc[k, "win1"] = int(np.nansum(delta > rope))
        out.loc[k, "win2"] = int(np.nansum(delta < -rope))
        out.loc[k, "ties"] = int(np.nansum(np.abs(delta) <= rope))
    return out


def compute_paired_differences(tab: pd.DataFrame, dbcol_name: str,
                               lrope_value: float) -> pd.DataFrame:
    """Wins/losses/ties using the *paired* local ROPE.

    For each data set (group of folds sharing the same ``dbcol_name``
    value), the per-fold differences between two algorithms are used to
    estimate a mean and standard deviation of the effect, which in turn
    defines the local ROPE for that data set. Results are then summed
    across data sets.
    """
    names = [c for c in tab.columns if c != dbcol_name]
    nalg = len(names)
    npairs = nalg * (nalg - 1) // 2
    totals = _empty_out(npairs)
    totals.loc[:, ["pi", "pj"]] = [(i, j) for i, j, _ in allpairs(nalg)]

    for _, group in tab.groupby(dbcol_name, sort=False):
        x = group[names]
        for i, j, k in allpairs(nalg):
            delta = x.iloc[:, i] - x.iloc[:, j]
            m = delta.mean()
            s = delta.std()
            rope = 0.0 if (s is None or pd.isna(s)) else lrope_value * s
            totals.loc[k, "win1"] += int(m > rope)
            totals.loc[k, "win2"] += int(m < -rope)
            totals.loc[k, "ties"] += int(rope * -1 <= m <= rope)
    return totals


def proc_ties(tab: pd.DataFrame, deal_with_ties: str) -> pd.DataFrame:
    """Redistribute the ``ties`` column according to ``deal_with_ties``.

    - ``s`` (spread): half (rounded up) of each tie goes to each algorithm.
    - ``a`` (add): the whole tie count is added to *both* algorithms.
    - ``f`` (forget): ties are dropped (set to 0), not added to either.
    - ``r`` (random): each tie is randomly assigned to one of the two.
    - ``d`` (davidson): left untouched, to be modeled explicitly.
    """
    tab = tab.copy()
    if deal_with_ties == "d":
        return tab
    if deal_with_ties == "a":
        tab["win1"] = tab["win1"] + tab["ties"]
        tab["win2"] = tab["win2"] + tab["ties"]
        tab["ties"] = 0
    elif deal_with_ties == "f":
        tab["ties"] = 0
    elif deal_with_ties == "s":
        half = np.ceil(tab["ties"] / 2).astype(int)
        tab["win1"] = tab["win1"] + half
        tab["win2"] = tab["win2"] + half
        tab["ties"] = 0
    elif deal_with_ties == "r":
        rng = np.random.default_rng()
        a = np.array([rng.integers(1, x, endpoint=True) if x > 0 else 0
                     for x in tab["ties"]])
        b = tab["ties"].to_numpy() - a
        tab["win1"] = tab["win1"] + a
        tab["win2"] = tab["win2"] + b
        tab["ties"] = 0
    return tab


def make_wintable(tabx: pd.DataFrame,
                  tabsd: Optional[pd.DataFrame] = None,
                  dbcol: Union[int, str, None] = 0,
                  lrope: Optional[bool] = None,
                  lrope_value: float = 0.4,
                  paired: bool = True,
                  deal_with_ties: Union[str, List[str]] = "spread") -> Wintable:
    """Build a :class:`Wintable` from a table of algorithm results.

    Parameters
    ----------
    tabx : pd.DataFrame
        Rows are data sets, columns are algorithms. One column may
        instead identify the data set (see ``dbcol``); in that case
        there can be multiple rows (folds) per data set.
    tabsd : pd.DataFrame, optional
        Standard deviation of the fold measures, aligned with ``tabx``
        (only used when ``tabx`` has no ``dbcol``/folds already).
    dbcol : int, str or None
        Column position (default 0, the first column) or name that
        identifies the data set. Use ``None`` (or ``-1``) when every row
        is already a distinct, single-measure data set.
    lrope : bool, optional
        Whether to use the local ROPE. If ``None``, it is turned on
        automatically whenever fold information (or ``tabsd``) is
        available.
    lrope_value : float
        The local ROPE threshold (default 0.4).
    paired : bool
        Whether to use the paired version of the local ROPE (only
        possible when ``dbcol`` is set and folds are present).
    deal_with_ties : str
        One of "spread" (default), "forget", "davidson", "add", "random"
        (or their first letters).

    Returns
    -------
    Wintable
    """
    if not isinstance(tabx, pd.DataFrame):
        raise TypeError("tabx must be a pandas DataFrame")
    if tabsd is not None and not isinstance(tabsd, pd.DataFrame):
        raise TypeError("tabsd must be a pandas DataFrame")

    ties_code = _normalize_ties(deal_with_ties)
    dbcol_name = _resolve_dbcol(list(tabx.columns), dbcol)
    names = [c for c in tabx.columns if c != dbcol_name]

    has_folds = multiple_folds(tabx, dbcol_name)

    if tabsd is None and not has_folds:
        # Only a mean table given.
        meanx = tabx[names] if dbcol_name is not None else tabx
        out = compute_differences_no_sd(meanx)
        lrope = False
        lrope_value = 0.0

    elif tabsd is not None:
        # Mean and sd tables given.
        meanx = tabx[names] if dbcol_name is not None else tabx
        if lrope is None or lrope:
            lrope = True
            out = compute_differences_with_sd(meanx, tabsd, lrope, lrope_value)
        else:
            out = compute_differences_no_sd(meanx)

    elif (paired is None or paired) and (lrope is None or lrope):
        # Only tabx with folds given -- paired local ROPE.
        out = compute_paired_differences(tabx, dbcol_name, lrope_value)
        lrope = True

    else:
        # Not paired -- compute mean and sd per data set, then use those.
        if lrope is None:
            lrope = True
        meanx = tabx.groupby(dbcol_name, sort=False)[names].mean().reset_index(drop=True)
        sdx = tabx.groupby(dbcol_name, sort=False)[names].std().reset_index(drop=True)
        out = compute_differences_with_sd(meanx, sdx, lrope, lrope_value)

    out2 = proc_ties(out, ties_code)
    return Wintable(
        table=out2.reset_index(drop=True),
        alg_names=names,
        lrope=bool(lrope),
        lrope_value=lrope_value,
        table_pre=out.reset_index(drop=True),
        ties_proc=ties_code,
        paired=bool(paired),
    )
