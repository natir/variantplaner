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

    df = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf", DATA_DIR / "grch38.92.csv")

    lf = extract.variants(df.lazy())

    polars.testing.assert_frame_equal(truth, lf, check_row_order=False, check_column_order=False)


def test_extract_genotypes() -> None:
    """Check extract genotypes."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")

    df = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf", DATA_DIR / "grch38.92.csv")

    vcf_header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    lf = extract.genotypes(df, io.vcf.format2expr(vcf_header, DATA_DIR / "no_info.vcf"))

    try:
        polars.testing.assert_frame_equal(truth, lf)
    except OverflowError:  # TODO remove this
        assert truth.columns == lf.columns
        assert truth.dtypes == lf.dtypes
        assert truth.width == lf.width


def test_extract_genotypes_no_gt() -> None:
    """Check extract genotypes without gt value."""
    truth = polars.scan_parquet(DATA_DIR / "no_info_no_gt.genotypes.parquet")

    input_path = DATA_DIR / "no_info_no_gt.vcf"

    lf = io.vcf.into_lazyframe(input_path, DATA_DIR / "grch38.92.csv")

    vcf_header = io.vcf.extract_header(input_path)

    format2expr = io.vcf.format2expr(vcf_header, input_path)

    lf = extract.genotypes(lf, format2expr, format_str="AD:DP:GQ")

    polars.testing.assert_frame_equal(truth, lf)


def test_extract_genotypes_without_genotypes() -> None:
    """Check extract genotypes failled if vcf not containts genotypes."""
    df = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf", DATA_DIR / "grch38.92.csv")
    vcf_header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")

    with pytest.raises(exception.NoGenotypeError):
        extract.genotypes(df.lazy(), io.vcf.format2expr(vcf_header, DATA_DIR / "no_info.vcf"))


def test_extract_merge_variant_genotypes() -> None:
    """Check merge_variants_genotypes."""
    vcf = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf", DATA_DIR / "grch38.92.csv")
    vcf_header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    sample_name = io.vcf.sample_index(vcf_header, DATA_DIR / "no_info.vcf")
    if sample_name is None:
        raise AssertionError  # pragma: no cover Not reachable code

    format2expr = io.vcf.format2expr(vcf_header, DATA_DIR / "no_info.vcf")

    variants = extract.variants(vcf)
    genotypes = extract.genotypes(vcf, format2expr, format_str="GT:AD:DP:GQ")

    merge = extract.merge_variants_genotypes(variants, genotypes, list(sample_name.keys()))

    assert merge.columns == [
        "id",
        "chr",
        "pos",
        "ref",
        "alt",
        "sample_1_gt",
        "sample_1_ad",
        "sample_1_dp",
        "sample_1_gq",
        "sample_2_gt",
        "sample_2_ad",
        "sample_2_dp",
        "sample_2_gq",
        "sample_3_gt",
        "sample_3_ad",
        "sample_3_dp",
        "sample_3_gq",
    ]

    assert merge.dtypes == [
        polars.UInt64,
        polars.Utf8,
        polars.UInt64,
        polars.Utf8,
        polars.Utf8,
        polars.UInt8,
        polars.List(polars.UInt16),
        polars.UInt16,
        polars.UInt16,
        polars.UInt8,
        polars.List(polars.UInt16),
        polars.UInt16,
        polars.UInt16,
        polars.UInt8,
        polars.List(polars.UInt16),
        polars.UInt16,
        polars.UInt16,
    ]
