"""Tests for the `manipulation` module."""

# std import
from __future__ import annotations

import pathlib

import polars

# 3rd party import
import pytest

# project import
from variantplanner import exception, io, manipulation

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_extract_variants() -> None:
    """Check extract variants."""
    truth = polars.read_parquet(DATA_DIR / "good.variants.parquet")

    df = io.vcf.into_lazyframe(DATA_DIR / "good.vcf")

    lf = manipulation.extract_variants(df.lazy())

    assert truth.frame_equal(lf.collect())


def test_extract_genotypes() -> None:
    """Check extract genotypes."""
    truth = polars.read_parquet(DATA_DIR / "good.genotypes.parquet")

    df = io.vcf.into_lazyframe(DATA_DIR / "good.vcf")

    lf = manipulation.extract_genotypes(df.lazy())

    assert truth.frame_equal(lf.collect())


def test_extract_genotypes_without_genotypes() -> None:
    """Check extract genotypes failled if vcf not containts genotypes."""
    df = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf")

    with pytest.raises(exception.NoGenotypeError):
        manipulation.extract_genotypes(df.lazy())
