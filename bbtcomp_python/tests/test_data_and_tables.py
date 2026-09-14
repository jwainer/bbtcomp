import pandas as pd

from bbtcomp import load_ll, make_wintable, table_wintable


def test_load_ll():
    ll = load_ll()
    assert isinstance(ll, pd.DataFrame)
    assert ll.shape[0] == 528
    assert "db" in ll.columns
    assert ll.shape[1] == 17


def test_table_wintable_pre_pos_both(ss):
    w = make_wintable(ss, lrope=True, paired=True)

    t_pos = table_wintable(w)
    assert list(t_pos.columns[:4]) == ["alg1", "alg2", "win1", "win2"]

    t_pre = table_wintable(w, which="pre")
    assert t_pre.columns[2] == "win1(pre)"

    t_both = table_wintable(w, which="both")
    assert t_both.shape[1] >= t_pos.shape[1] + 3
