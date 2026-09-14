import os

import pandas as pd
import pytest

_HERE = os.path.dirname(__file__)


@pytest.fixture
def ss() -> pd.DataFrame:
    """Fold-level results: 5 algorithms, 20 data sets, 4 folds each."""
    return pd.read_csv(os.path.join(_HERE, "ss.csv"))


@pytest.fixture
def ssmean() -> pd.DataFrame:
    """Mean-only results (with a db column), same 5 algorithms/20 data sets."""
    return pd.read_csv(os.path.join(_HERE, "ssmean.csv"))


@pytest.fixture
def ssmeanx(ssmean) -> pd.DataFrame:
    """``ssmean`` with a handful of rows dropped (matches the R test suite)."""
    drop_idx = [2, 4, 6, 7, 11, 19]  # 0-indexed equivalents of R's c(3,5,7,8,12,20)
    return ssmean.drop(index=drop_idx).reset_index(drop=True)


@pytest.fixture
def amean(ssmean) -> pd.DataFrame:
    """``ssmean`` without the db column, as a bare numeric frame."""
    return ssmean.drop(columns=["db"])


@pytest.fixture
def sssd() -> pd.DataFrame:
    """Standard deviations aligned with ``ssmean`` (with a db column)."""
    return pd.read_csv(os.path.join(_HERE, "sssd.csv"))
