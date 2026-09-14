"""Access to the bundled ``ll`` dataset."""

from __future__ import annotations

import importlib.resources

import pandas as pd


def load_ll() -> pd.DataFrame:
    """Load the bundled ``ll`` dataset: accuracy of 17 classifiers on 132
    data sets (4 cross-validation folds each), as used in the BBT paper.

    Returns
    -------
    pd.DataFrame
        526 rows x 17 columns; the first column ``db`` is the data-set
        name, the rest are per-fold accuracy values for each algorithm.
    """
    ref = importlib.resources.files("bbtcomp") / "data" / "ll.csv"
    with importlib.resources.as_file(ref) as path:
        return pd.read_csv(path)
