"""Tests for the `exception` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
# project import
from variantplaner import exception


def test_notavcferror() -> None:
    """Check exception NotAVCFError."""
    e = exception.NotAVCFError(pathlib.Path("test"))

    assert f"{e}" == "File test seems not be a valid vcf file."


def test_nogenotypeerror() -> None:
    """Check exception NoGenotypeError."""
    e = exception.NoGenotypeError()

    assert f"{e}" == "LazyFrame seems not contains genotypes."
