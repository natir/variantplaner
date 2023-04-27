"""Tests for the `extract` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars
import polars.testing
import pytest

# project import
from variantplaner import exception, extract, io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_extract_variants() -> None:
    """Check extract variants."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.variants.parquet")

    df = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf")

    lf = extract.variants(df.lazy())

    polars.testing.assert_frame_equal(truth, lf)


def test_extract_genotypes() -> None:
    """Check extract genotypes."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")

    df = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf")

    vcf_header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    lf = extract.genotypes(df, io.vcf.format2expr(vcf_header, DATA_DIR / "no_info.vcf"))

    polars.testing.assert_frame_equal(truth, lf)


def test_extract_genotypes_without_genotypes() -> None:
    """Check extract genotypes failled if vcf not containts genotypes."""
    df = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf")
    vcf_header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")

    with pytest.raises(exception.NoGenotypeError):
        extract.genotypes(df.lazy(), io.vcf.format2expr(vcf_header, DATA_DIR / "no_info.vcf"))
