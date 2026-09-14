"""bbtcomp: Bayesian Bradley-Terry models to compare multiple algorithms
on multiple data sets.

A Python port of the R package of the same name -- see
https://github.com/jwainer/bbtcomp and the accompanying paper for the
statistical background (Bayesian Bradley-Terry model, ROPE and "local
ROPE").
"""

from typing import Optional

from ._bbtcomp import bbtcomp
from ._check import check_setup
from ._convergence import convergence_check
from ._data import load_ll
from ._mcmc import BBTModel, mcmcbbt
from ._plots import plot_ppc, plot_pwin
from ._pwin_table import PwinTable
from ._tables import table_ppc, table_pwin, table_wintable
from ._utils import is_bbt_model
from ._waic_loo import get_loo, get_waic
from ._wintable import Wintable, make_wintable

__version__ = "0.4.0"

#: Package-wide default for ``mcmcbbt``'s/``bbtcomp``'s ``output_dir``
#: argument, used whenever a call doesn't pass ``output_dir`` explicitly.
#: ``None`` (the default) means "use a fresh temporary directory per
#: call". Python has no R-style ``options()`` global settings registry,
#: so this plays the role of the R package's
#: ``options(bbtcomp.dir = "~/.bbtcomp")``: set it once per session,
#: e.g. ``bbtcomp.DEFAULT_OUTPUT_DIR = "~/.bbtcomp"``, and every
#: subsequent call that doesn't pass its own ``output_dir`` reuses that
#: directory instead of recompiling the Stan model each time. See the
#: README's "Output / auxiliary directory" section.
DEFAULT_OUTPUT_DIR: Optional[str] = None

__all__ = [
    "bbtcomp",
    "check_setup",
    "DEFAULT_OUTPUT_DIR",
    "make_wintable",
    "mcmcbbt",
    "plot_pwin",
    "plot_ppc",
    "table_pwin",
    "table_ppc",
    "table_wintable",
    "convergence_check",
    "get_waic",
    "get_loo",
    "load_ll",
    "Wintable",
    "BBTModel",
    "is_bbt_model",
    "PwinTable",
]
