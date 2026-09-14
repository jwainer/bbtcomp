"""Unit tests for the PwinTable display class used by table_pwin()."""

from bbtcomp import PwinTable


def _sample():
    return PwinTable({
        "larger": ["a", "b"],
        "smaller": ["c", "d"],
        "mean": [0.6, 0.7],
    })


def test_pwin_table_is_a_dataframe_with_larger_smaller_columns():
    df = _sample()
    assert list(df.columns) == ["larger", "smaller", "mean"]
    assert df["larger"].tolist() == ["a", "b"]
    assert df["smaller"].tolist() == ["c", "d"]


def test_pwin_table_repr_shows_combined_pair_column():
    df = _sample()
    printed = repr(df)
    assert "a > c" in printed
    assert "b > d" in printed
    # the header should show "pair", not the raw larger/smaller columns
    header = printed.splitlines()[0]
    assert "pair" in header
    assert "larger" not in header
    assert "smaller" not in header


def test_pwin_table_slicing_preserves_class_and_behavior():
    df = _sample()
    sub = df[["larger", "smaller"]]
    assert isinstance(sub, PwinTable)
    assert "a > c" in repr(sub)


def test_pwin_table_falls_back_when_larger_smaller_missing():
    df = _sample()
    sub = df[["mean"]]
    # no larger/smaller columns left: display just falls back to a plain
    # DataFrame repr rather than erroring
    printed = repr(sub)
    assert "mean" in printed
