"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars

# project import
from variantplanner import io

DATA_DIR = pathlib.Path(__file__).parent / "data"

def test_info2expr() -> None:
    """Check info2expr"""

    expressions = io.vcf.info2expr(DATA_DIR / "no_info.vcf")

    assert expressions == list()


def test_into_lazyframe() -> None:
    """Check into_lazyframe."""
    pass
