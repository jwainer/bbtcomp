"""Tests for the CmdStan environment check (bbtcomp.check_setup)."""

import pytest

from bbtcomp import check_setup


def test_check_setup_returns_bool():
    # Whatever the environment, check_setup() should never raise on its
    # own and should return a plain bool.
    result = check_setup()
    assert isinstance(result, bool)


def test_check_setup_raise_on_error_when_missing(monkeypatch):
    import cmdstanpy

    def _boom():
        raise RuntimeError("no CmdStan installation found")

    monkeypatch.setattr(cmdstanpy, "cmdstan_path", _boom)

    assert check_setup() is False
    with pytest.raises(RuntimeError, match="CmdStan could not be found"):
        check_setup(raise_on_error=True)


def test_check_setup_true_when_present(monkeypatch):
    import cmdstanpy

    monkeypatch.setattr(cmdstanpy, "cmdstan_path", lambda: "/fake/cmdstan-1.2.3")

    assert check_setup() is True
