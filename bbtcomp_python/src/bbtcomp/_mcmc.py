"""Fit the Bayesian Bradley-Terry model with Stan (via cmdstanpy)."""

from __future__ import annotations

import importlib.resources
import os
import shutil
import tempfile
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import pandas as pd

from ._check import _ensure_cmdstan
from ._utils import allpairs
from ._wintable import Wintable

_STAN_FULL = "bbt-full.stan"


def _stan_file(name: str = _STAN_FULL) -> str:
    ref = importlib.resources.files("bbtcomp") / "stan" / name
    return str(ref)


@dataclass
class BBTModel:
    """A fitted BBT model.

    Attributes
    ----------
    model : cmdstanpy.CmdStanMCMC
        The cmdstanpy fit object (posterior samples).
    davidson : bool
        Whether the Davidson extension (explicit tie modeling) was used.
    wintable : Wintable
        The wintable the model was fit on.
    """

    model: "object"
    davidson: bool
    wintable: Wintable


def mcmcbbt(wintable: Wintable,
           hyper_prior: int = 0,
           scale: float = 0.5,
           use_davidson: bool = False,
           output_dir: Optional[str] = None,
           **kwargs) -> BBTModel:
    """Run MCMC sampling on a :class:`Wintable` to fit the BBT model.

    Parameters
    ----------
    wintable : Wintable
        As produced by :func:`bbtcomp.make_wintable`.
    hyper_prior : int
        0 = lognormal(0, 0.5) [default], 1 = lognormal(0, scale),
        2 = cauchy(0, scale), 3 = normal(0, scale).
    scale : float
        Scale of the hyper-prior for ``sigma`` (used when ``hyper_prior``
        is 1, 2 or 3).
    use_davidson : bool
        Whether to use the Davidson extension of the model, which
        explicitly models ties instead of having them pre-processed away.
    output_dir : str, optional
        Directory used by cmdstanpy to store the compiled model and the
        sampler's csv output. Defaults to ``bbtcomp.DEFAULT_OUTPUT_DIR``
        if that package-level variable has been set (Python has no
        R-style ``options()`` registry, so this plays the role of R's
        ``options(bbtcomp.dir = ...)`` -- see the README), otherwise to
        a fresh temporary directory.
    **kwargs
        Extra keyword arguments forwarded to
        ``cmdstanpy.CmdStanModel.sample`` (e.g. ``chains``, ``seed``,
        ``iter_sampling``).

    Returns
    -------
    BBTModel

    Raises
    ------
    RuntimeError
        If CmdStan cannot be found (`pip install` only installs the
        cmdstanpy Python package -- CmdStan itself, the compiled Stan
        backend, must be set up separately; see
        :func:`bbtcomp.check_setup`).
    """
    import cmdstanpy

    _ensure_cmdstan()

    tab = wintable.table
    K = len(wintable.alg_names)
    N = len(tab)

    data = dict(
        hyp=int(hyper_prior),
        scale=float(scale),
        use_davidson=int(use_davidson),
        K=K,
        N=N,
        player1=(tab["pi"].to_numpy() + 1).tolist(),
        player2=(tab["pj"].to_numpy() + 1).tolist(),
        win1=tab["win1"].astype(int).tolist(),
        win2=tab["win2"].astype(int).tolist(),
        ties=tab["ties"].astype(int).tolist(),
    )

    # bbtcomp has no R-style options() registry; a package-level default
    # (set as `bbtcomp.DEFAULT_OUTPUT_DIR = ...`) plays that role instead.
    # Imported lazily and read from the top-level package here (rather than
    # via `from . import DEFAULT_OUTPUT_DIR`) so that reassigning the
    # attribute on the package after import is actually picked up.
    import bbtcomp as _pkg

    outdir = output_dir or _pkg.DEFAULT_OUTPUT_DIR or tempfile.mkdtemp(prefix="bbtcomp_")
    outdir = os.path.expanduser(outdir)
    os.makedirs(outdir, exist_ok=True)

    # cmdstanpy.CmdStanModel has no equivalent of cmdstanr's `dir=`: it
    # compiles the binary next to the .stan file it is given. Copy the
    # bundled .stan file into outdir first so the compiled binary (and
    # its cache) lands there rather than inside the installed package,
    # which may not be writable.
    stan_src = _stan_file()
    stan_copy = os.path.join(outdir, os.path.basename(stan_src))
    if not os.path.exists(stan_copy):
        shutil.copy(stan_src, stan_copy)

    model = cmdstanpy.CmdStanModel(stan_file=stan_copy)
    fit = model.sample(data=data, output_dir=outdir, **kwargs)

    return BBTModel(model=fit, davidson=bool(use_davidson), wintable=wintable)


def get_pwin(mod: BBTModel, selected: Optional[List[str]] = None,
            control: Optional[str] = None):
    """Posterior samples of P(alg_i beats alg_j) for every requested pair.

    Returns
    -------
    (np.ndarray, list of str, list of str, list of str)
        A (draws x pairs) array; the matching "name1 > name2" labels;
        and the same two algorithm names again as separate ``larger``
        and ``smaller`` lists (``larger`` is always the one with the
        higher posterior mean ability), for callers that want them as
        two programmatically usable fields instead of parsing them back
        out of the combined label. Algorithms are ordered by decreasing
        posterior mean ability.
    """
    beta = mod.model.stan_variable("beta")  # (draws, K)
    names = np.array(mod.wintable.alg_names)
    order = np.argsort(-beta.mean(axis=0))
    ordered_names = names[order]
    w = np.exp(beta[:, order])

    if selected is not None:
        keep = [i for i, n in enumerate(ordered_names) if n in selected]
        ordered_names = ordered_names[keep]
        w = w[:, keep]

    n = len(ordered_names)

    def toprob(a, b):
        return a / (a + b)

    larger: List[str] = []
    smaller: List[str] = []

    if control is None or control not in ordered_names:
        out = np.empty((w.shape[0], n * (n - 1) // 2))
        out_names = []
        for i, j, k in allpairs(n):
            out[:, k] = toprob(w[:, i], w[:, j])
            out_names.append(f"{ordered_names[i]} > {ordered_names[j]}")
            larger.append(str(ordered_names[i]))
            smaller.append(str(ordered_names[j]))
        return out, out_names, larger, smaller

    ci = int(np.where(ordered_names == control)[0][0])
    out = np.empty((w.shape[0], n - 1))
    out_names = []
    k = 0
    for i in range(n):
        if i == ci:
            continue
        if i < ci:
            out[:, k] = toprob(w[:, i], w[:, ci])
            out_names.append(f"{ordered_names[i]} > {ordered_names[ci]}")
            larger.append(str(ordered_names[i]))
            smaller.append(str(ordered_names[ci]))
        else:
            out[:, k] = toprob(w[:, ci], w[:, i])
            out_names.append(f"{ordered_names[ci]} > {ordered_names[i]}")
            larger.append(str(ordered_names[ci]))
            smaller.append(str(ordered_names[i]))
        k += 1
    return out, out_names, larger, smaller


def _pair_labels(alg_names: List[str]) -> List[str]:
    labels = []
    for i, j, _ in allpairs(len(alg_names)):
        labels.append(f"{alg_names[i]} > {alg_names[j]}")
    return labels


def aux_get_win1_rep(mod: BBTModel) -> pd.DataFrame:
    """Posterior-predictive draws of ``win1`` for every pair."""
    rep = mod.model.stan_variable("win1_rep")  # (draws, N)
    cols = _pair_labels(mod.wintable.alg_names)
    return pd.DataFrame(rep, columns=cols)


def aux_get_ties_rep(mod: BBTModel) -> pd.DataFrame:
    """Posterior-predictive draws of ``tie_rep`` (Davidson model only)."""
    rep = mod.model.stan_variable("tie_rep")  # (draws, N)
    names = mod.wintable.alg_names
    cols = []
    for i, j, _ in allpairs(len(names)):
        cols.append(f"{names[i]} = {names[j]}")
    return pd.DataFrame(rep, columns=cols)
