import pandas as pd
import pandas.testing as pdt
import pytest

from bbtcomp import Wintable, make_wintable


def test_no_ties_ssmeanx_all_the_same(ssmeanx):
    """When every data set has a single row (no folds), lrope/paired/tie
    settings that don't apply should all be silently ignored and give
    the same result -- mirrors the R testthat suite's first test."""
    w1 = make_wintable(ssmeanx)
    w2 = make_wintable(ssmeanx, lrope=False)
    w3 = make_wintable(ssmeanx, lrope=True, paired=False)
    w4 = make_wintable(ssmeanx, lrope=True, paired=True)
    w5 = make_wintable(ssmeanx, lrope=False, deal_with_ties="spread")
    w6 = make_wintable(ssmeanx, lrope=False, deal_with_ties="random")
    w7 = make_wintable(ssmeanx, lrope=False, deal_with_ties="forget")
    w8 = make_wintable(ssmeanx, lrope=False, deal_with_ties="davidson")

    for w in (w2, w3, w4, w5, w7, w8):
        pdt.assert_frame_equal(w1.table, w.table)
    # w6 uses "random" tie-breaking but ssmeanx (real-valued measures) has
    # no actual ties, so it should still match.
    pdt.assert_frame_equal(w1.table, w6.table)


def test_tabmean_and_tabsd(ss, ssmean, sssd, amean):
    a = make_wintable(amean, sssd.drop(columns=["db"]), lrope=True, dbcol=None, paired=False)
    b = make_wintable(ss, lrope=True, paired=False)
    pdt.assert_frame_equal(a.table, b.table)

    a = make_wintable(ssmean.drop(columns=["db"]), sssd.drop(columns=["db"]),
                      lrope=True, dbcol=None, paired=False)
    pdt.assert_frame_equal(a.table, b.table)

    a = make_wintable(ssmean.drop(columns=["db"]), sssd.drop(columns=["db"]),
                      dbcol=None, lrope=True, lrope_value=0.3, paired=False)
    b = make_wintable(ss, lrope=True, paired=False, lrope_value=0.3)
    pdt.assert_frame_equal(a.table, b.table)

    a = make_wintable(ssmean.drop(columns=["db"]), sssd.drop(columns=["db"]),
                      dbcol=None, lrope=True, deal_with_ties="davidson", paired=False)
    b = make_wintable(ss, lrope=True, paired=False, deal_with_ties="davidson")
    pdt.assert_frame_equal(a.table, b.table)

    a = make_wintable(amean, sssd.drop(columns=["db"]),
                      dbcol=None, lrope=True, deal_with_ties="forget", paired=False)
    b = make_wintable(ss, lrope=True, paired=False, deal_with_ties="forget")
    pdt.assert_frame_equal(a.table, b.table)


@pytest.mark.parametrize("kwargs", [
    dict(lrope=False, paired=True),
    dict(lrope=True, paired=False),
    dict(deal_with_ties="spread"),
    dict(deal_with_ties="forget"),
    dict(deal_with_ties="davidson"),
    dict(deal_with_ties="random", lrope=False),
])
def test_variations_all_work(ss, kwargs):
    w = make_wintable(ss, **kwargs)
    assert isinstance(w, Wintable)
    assert set(w.table.columns) == {"pi", "pj", "win1", "win2", "ties"}
    assert w.alg_names == [c for c in ss.columns if c != "db"]
