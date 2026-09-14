"""Small numerical/statistical helpers shared across the package."""

from __future__ import annotations

import math
from itertools import combinations
from typing import Iterator, Tuple

import numpy as np


def allpairs(n: int) -> Iterator[Tuple[int, int, int]]:
    """Yield (i, j, k) for every 0-indexed pair i < j in range(n).

    k is a running counter (0, 1, 2, ...) giving the position of the
    pair, mirroring the R package's ``allpairs`` helper.
    """
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        yield i, j, k


def is_bbt_model(modout) -> bool:
    """Return True if ``modout`` looks like a fitted BBT model."""
    return hasattr(modout, "model") and modout.model is not None


def hdi(x: np.ndarray, credmass: float = 0.89) -> Tuple[float, float]:
    """Highest Density Interval of a sample.

    A direct port of the algorithm used by R's ``HDInterval::hdi`` (and
    the ``newhdi`` helper in the original python prototype): among all
    intervals that contain ``credmass`` fraction of the (sorted) sample,
    return the narrowest one.
    """
    x = np.sort(np.asarray(x))
    n = len(x)
    exclude = n - math.floor(n * credmass) - 1
    if exclude <= 0:
        return float(x[0]), float(x[-1])
    low_poss = x[0:exclude]
    upp_poss = x[(n - exclude):n]
    best = np.argmin(upp_poss - low_poss)
    return float(low_poss[best]), float(upp_poss[best])
