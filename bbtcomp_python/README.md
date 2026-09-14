# bbtcomp

A Bayesian Bradley-Terry (BBT) model to compare multiple algorithms on
multiple data sets, based on any metric (accuracy, F1, AUC, RMSE, ...).

This is the Python port of the [`bbtcomp` R
package](https://github.com/jwainer/bbtcomp); see the accompanying paper
for the statistical background (the BBT model, the region of practical
equivalence -- ROPE -- and "local ROPE").

## Installation

`bbtcomp` needs [CmdStan](https://mc-stan.org/cmdstanpy/) to fit the
model. Install `cmdstanpy` and CmdStan first, following
<https://mc-stan.org/cmdstanpy/installation.html>, then:

```bash
pip install bbtcomp
```

Optional extras:

```bash
pip install "bbtcomp[plots]"  # plot_pwin / plot_ppc (matplotlib)
pip install "bbtcomp[waic]"   # get_waic / get_loo (arviz)
```

`pip install bbtcomp` installs the `cmdstanpy` *Python package*, but it
cannot install or verify CmdStan itself (the compiled Stan backend) --
CmdStan isn't a Python package, so pip has no way to know about it. After
installing, run this once to confirm CmdStan is actually set up:

```python
from bbtcomp import check_setup

check_setup()
# CmdStan found at: /path/to/cmdstan-x.y.z
```

If CmdStan isn't found, `check_setup()` prints (or, with
`check_setup(raise_on_error=True)`, raises) clear instructions instead of
you having to decode a `cmdstanpy` traceback later. `bbtcomp()` /
`mcmcbbt()` also run this check automatically and fail fast with the same
message if CmdStan is missing.

## Data format

Data is a `pandas.DataFrame` where rows are data sets and columns are
algorithms. The column names are used as algorithm names.

One column may instead be a data-set identifier (the `dbcol` parameter,
default: the first column). When present, there can be multiple rows
per data set -- the results of different cross-validation folds -- which
lets `bbtcomp` use the "local ROPE" to turn some close wins into ties.

Data can be:

- a DataFrame with a `dbcol` identifying the data set, and potentially
  several rows (folds) per data set;
- a DataFrame or array with no `dbcol` -- one row per data set, one
  value per algorithm (usually the mean over folds);
- two DataFrames/arrays with the same column names, one with the mean
  measure per algorithm/data set and the other with the standard
  deviation across folds -- this still allows the (non-paired) local
  ROPE to be used.

Higher values must mean "better" (accuracy, AUC, F1, ...); negate
error-like metrics (RMSE, MAE, ...) first. Missing entries should be
`NaN`.

The package bundles the `ll` data set used in the paper: 17 classifiers
on 132 data sets (4 folds each).

```python
from bbtcomp import load_ll

ll = load_ll()
ss = ll.iloc[0:80, 0:6]  # first 20 data sets (4 folds), 5 algorithms
```

## Before you start

### Running chains in parallel

By default `cmdstanpy` samples chains sequentially. Pass `chains` and
`parallel_chains` (forwarded to `CmdStanModel.sample`) to run them in
parallel, e.g. `bbtcomp(ss, parallel_chains=4)`.

### Output / auxiliary directory

CmdStan compiles the Stan model to a binary and writes one `.csv` file
per chain with the samples. Pass `output_dir=...` to `bbtcomp` /
`mcmcbbt` to reuse a fixed directory across calls (so the model isn't
recompiled every time); leave it unset to use a fresh temporary
directory each call.

The R package lets you set this once per session with
`options(bbtcomp.dir = "~/.bbtcomp")` instead of passing `dir=...` to
every call. Python has no equivalent of R's `options()` -- there is no
global settings registry -- so the direct translation is a plain
package-level variable, `bbtcomp.DEFAULT_OUTPUT_DIR`, which every call
that omits `output_dir` falls back to:

```python
import bbtcomp

bbtcomp.DEFAULT_OUTPUT_DIR = "~/.bbtcomp"

y = bbtcomp.bbtcomp(ss)         # reuses ~/.bbtcomp, no output_dir needed
z = bbtcomp.bbtcomp(ss, iter_sampling=2000)  # still reuses it
```

An explicit `output_dir=...` on a given call always takes precedence
over `DEFAULT_OUTPUT_DIR`; if neither is set, a fresh temporary
directory is used, same as before.

## Basic use

```python
from bbtcomp import bbtcomp, table_pwin, table_ppc, table_wintable, plot_pwin, plot_ppc

x = bbtcomp(ss)

table_pwin(x)
#          pair  mean  delta  above.50  in.rope
# 0   lgbm > rf  0.53   0.22      0.66     0.49
# ...

table_pwin(x, control="rf")
table_pwin(x, columns=["median", "low", "high"], hdi_prob=0.95)

# table_pwin() returns a PwinTable (a pandas.DataFrame subclass): the two
# algorithm names of each pair are separate "larger"/"smaller" columns
# (the printed "pair" column above is just those two combined for display,
# e.g. "lgbm > rf" means larger="lgbm", smaller="rf"), so code can use
# them directly instead of parsing the printed string back apart:
tp = table_pwin(x)
tp["larger"]                        # -> a plain pandas Series
tp[tp["larger"] == "lgbm"]          # every pair lgbm won on average

plot_pwin(x)
plot_ppc(x)

table_wintable(x.wintable)                 # ties already processed
table_wintable(x.wintable, which="pre")    # ties shown explicitly
table_wintable(x.wintable, which="both")
```

For the mean-only format (no `dbcol`):

```python
x2 = bbtcomp(ssmean, dbcol=None)
```

For the mean + standard-deviation format:

```python
x3 = bbtcomp(ssmean, ssstd, dbcol=None)
```

### Local ROPE

When fold-level data is available, `bbtcomp` can convert some
"wins"/"losses" into ties whenever the difference between two
algorithms' means is small relative to the fold-to-fold variance (the
"local ROPE"). Two variants are supported: `paired=True` (default, folds
are the same test set for every algorithm) and `paired=False`.

```python
y = bbtcomp(ss, lrope=True, paired=True)
z = bbtcomp(ss, lrope=True, paired=False)
```

### Convergence check

```python
from bbtcomp import convergence_check

y = bbtcomp(ss, output_dir="~/.bbtcomp")  # keep the sampler csv files
convergence_check(y)
```

### Model comparison

```python
from bbtcomp import get_waic, get_loo

get_waic(y)
get_loo(y)
```

## Development

```bash
pip install -e ".[test]"
pytest
```

Tests that actually run MCMC (`tests/test_bbtcomp.py`) require a working
CmdStan installation and are skipped automatically otherwise.
