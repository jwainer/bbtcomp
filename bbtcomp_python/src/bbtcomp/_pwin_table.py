"""A small ``pandas.DataFrame`` subclass used by :func:`bbtcomp.table_pwin`.

``table_pwin()`` needs to store the two algorithm names of each compared
pair as separate, programmatically usable fields (``larger``/``smaller``,
``larger`` being the one with the higher estimated ability) rather than a
single pre-formatted ``"A > B"`` string -- so that code can do things like
``tp["larger"]`` or ``tp[tp["larger"] == "svm"]`` without having to parse a
string back apart. At the same time, the table should still *look* the way
it always did when printed: one readable ``"pair"`` column showing
``"A > B"``.

:class:`PwinTable` is an ordinary ``pandas.DataFrame`` (so every DataFrame
method -- indexing, ``.query()``, ``.to_csv()``, and so on -- works exactly
as normal, and it still has ``larger``/``smaller`` as real columns) with
just its display (``__repr__`` and ``_repr_html_``, used by the terminal
and Jupyter respectively) overridden to show ``larger``/``smaller`` merged
into a single ``"pair"`` column.
"""

from __future__ import annotations

import pandas as pd


class PwinTable(pd.DataFrame):
    """DataFrame returned by :func:`bbtcomp.table_pwin`.

    Has ``larger`` and ``smaller`` columns (plus the requested summary
    columns) for programmatic use, but displays them combined as a single
    ``"pair"`` column (e.g. ``"A > B"``) when printed or shown in Jupyter.
    """

    # tells pandas to keep returning PwinTable (rather than a plain
    # DataFrame) from operations like slicing, .copy(), etc.
    @property
    def _constructor(self):
        return PwinTable

    def _display_frame(self) -> pd.DataFrame:
        """Plain DataFrame used for display: ``larger``/``smaller`` merged
        into a single ``"pair"`` column. Falls back to showing the table
        as-is if it's been sliced down to something that no longer has
        both columns (e.g. ``tp[["mean"]]``)."""
        if "larger" not in self.columns or "smaller" not in self.columns:
            return pd.DataFrame(self)
        disp = pd.DataFrame(self).copy()
        pair = disp["larger"].astype(str) + " > " + disp["smaller"].astype(str)
        disp = disp.drop(columns=["larger", "smaller"])
        disp.insert(0, "pair", pair)
        return disp

    def __repr__(self) -> str:
        return repr(self._display_frame())

    def _repr_html_(self):  # pragma: no cover - exercised only in Jupyter
        disp = self._display_frame()
        to_html = getattr(disp, "_repr_html_", None)
        return to_html() if to_html is not None else disp.to_html()
