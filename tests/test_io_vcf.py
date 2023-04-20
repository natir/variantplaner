"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
# project import
from variantplanner import io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_info2expr_no_info_vcf() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_info.vcf")

    assert expressions == []


def test_info2expr() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_genotypes.vcf")

    assert len(expressions) == 21
