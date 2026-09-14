"""End-to-end tests that actually run MCMC via cmdstanpy/cmdstan.

These mirror ``tests/testthat/test-bbtcomp.R`` from the R package. They
need a working CmdStan installation (see
https://mc-stan.org/cmdstanpy/installation.html) and are skipped
automatically if one isn't available, since compiling/installing CmdStan
is a heavy, environment-specific step outside this package's control.
"""

import pytest

from bbtcomp import bbtcomp, is_bbt_model, plot_ppc, plot_pwin, table_ppc, table_pwin


def _cmdstan_available() -> bool:
    try:
        import cmdstanpy
        cmdstanpy.cmdstan_path()
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _cmdstan_available(),
    reason="CmdStan is not installed/configured (see cmdstanpy install docs)",
)

_SAMPLE_KWARGS = dict(chains=2, iter_warmup=200, iter_sampling=200, seed=1234)


def test_bbtcomp_basic(ss):
    m1 = bbtcomp(ss, lrope=True, **_SAMPLE_KWARGS)
    assert is_bbt_model(m1)

    tp = table_pwin(m1)
    assert tp.shape == (10, 6)
    assert "larger" in tp.columns and "smaller" in tp.columns
    printed = repr(tp)
    assert " > " in printed
    assert "larger" not in printed.splitlines()[0]

    fig1 = plot_pwin(m1)
    assert fig1 is not None

    fig2 = plot_ppc(m1)
    assert fig2 is not None

    tppc = table_ppc(m1)
    assert tppc.shape == (4, 2)


@pytest.mark.parametrize("kwargs", [
    dict(lrope=True, paired=False),
    dict(lrope_value=0.2, deal_with_ties="forget"),
    dict(deal_with_ties="davidson"),
    dict(lrope=False, hyper_prior=1, scale=2.0),
])
def test_bbtcomp_variations(ss, kwargs):
    m = bbtcomp(ss, **kwargs, **_SAMPLE_KWARGS)
    assert is_bbt_model(m)
