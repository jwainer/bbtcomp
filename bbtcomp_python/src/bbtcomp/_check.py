"""Environment checks: bbtcomp needs CmdStan, but ``pip install bbtcomp``
cannot verify or install it -- CmdStan is a compiled C++ program, not a
Python package, so pip has no way to declare or check that dependency.
This module gives a clear, actionable error instead of a deep cmdstanpy
traceback when CmdStan can't be found.
"""

from __future__ import annotations

_INSTALL_HINT = (
    "CmdStan could not be found.\n\n"
    "bbtcomp needs CmdStan (the Stan C++ backend) in addition to the "
    "cmdstanpy Python package -- `pip install bbtcomp` installs cmdstanpy "
    "but cannot install or check for CmdStan itself, since it isn't a "
    "Python package.\n\n"
    "To fix this, do ONE of the following:\n\n"
    "  1. If you don't have CmdStan yet, install it with:\n"
    "         python -m cmdstanpy.install_cmdstan\n"
    "     (this downloads and compiles CmdStan; needs a C++ compiler and\n"
    "     takes several minutes). See\n"
    "     https://mc-stan.org/cmdstanpy/installation.html for details.\n\n"
    "  2. If you already have a CmdStan installation (for example, one\n"
    "     installed for the R package via cmdstanr::install_cmdstan()),\n"
    "     point cmdstanpy at it instead of installing a second copy:\n"
    "         import cmdstanpy\n"
    "         cmdstanpy.set_cmdstan_path(\"/path/to/cmdstan-x.y.z\")\n"
    "     or set the CMDSTAN environment variable to that path before\n"
    "     running Python. In R, `cmdstanr::cmdstan_path()` prints the\n"
    "     path of an existing installation.\n"
)


def check_setup(raise_on_error: bool = False) -> bool:
    """Check whether CmdStan is installed and discoverable by cmdstanpy.

    Run this right after installing bbtcomp (before loading any data) to
    get a clear, actionable message if CmdStan is missing, instead of
    discovering it later as a cryptic error buried inside a model-fitting
    call.

    Parameters
    ----------
    raise_on_error : bool
        If True, raise a ``RuntimeError`` with setup instructions when
        CmdStan can't be found, instead of just returning False.

    Returns
    -------
    bool
        True if CmdStan was found (a message with its path is printed).
        False if it wasn't (a message with setup instructions is printed,
        unless ``raise_on_error`` is True, in which case it raises).

    Examples
    --------
    >>> from bbtcomp import check_setup
    >>> check_setup()  # doctest: +SKIP
    """
    try:
        import cmdstanpy
    except ImportError as exc:  # pragma: no cover - cmdstanpy is a hard dep
        msg = (
            "The `cmdstanpy` package is not installed. It should have been "
            "installed automatically with bbtcomp -- try `pip install "
            "--force-reinstall bbtcomp` or `pip install cmdstanpy`."
        )
        if raise_on_error:
            raise RuntimeError(msg) from exc
        print(msg)
        return False

    try:
        path = cmdstanpy.cmdstan_path()
    except Exception:
        if raise_on_error:
            raise RuntimeError(_INSTALL_HINT)
        print(_INSTALL_HINT)
        return False

    print(f"CmdStan found at: {path}")
    return True


def _ensure_cmdstan() -> None:
    """Internal: used by mcmcbbt() to fail fast with a clear message."""
    import cmdstanpy

    try:
        cmdstanpy.cmdstan_path()
    except Exception as exc:
        raise RuntimeError(_INSTALL_HINT) from exc
