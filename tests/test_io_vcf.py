"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import filecmp

# 3rd party import
import pathlib

import polars
import pytest

# project import
from variantplanner import exception, io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_info2expr_no_info_vcf() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_info.vcf")

    assert expressions == []


def test_info2expr() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_genotypes.vcf")

    assert len(expressions) == 21


def test_sample_index() -> None:
    """Check sample index."""
    truth = {"sample_1": 0, "sample_2": 1, "sample_3": 2}

    sample2idx = io.vcf.sample_index(DATA_DIR / "no_info.vcf")

    assert sample2idx
    assert len(truth) == len(sample2idx)
    assert all(v == sample2idx[k] for k, v in truth.items())


def test_sample_index_no_genotypes() -> None:
    """Check sample index."""
    assert io.vcf.sample_index(DATA_DIR / "no_genotypes.vcf") is None


def test_sample_index_exception() -> None:
    """Check sample index exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.sample_index(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.sample_index(DATA_DIR / "only_header.vcf")


def test_into_lazyframe() -> None:
    """Check into lazyframe."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf")

    polars.testing.assert_frame_equal(truth, lf)

    truth = polars.scan_parquet(DATA_DIR / "no_genotypes.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf")

    polars.testing.assert_frame_equal(truth, lf)


def test_into_lazyframe_exception() -> None:
    """Check into_lazyframe exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "only_header.vcf")


def test_build_rename_column() -> None:
    """Check build_rename_column."""
    assert io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt") == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": ".",
        "FILTER": ".",
        "INFO": {},
    }

    assert io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt", "quality", "FILTER", {"GENE": "gene_name"}) == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": "quality",
        "FILTER": "FILTER",
        "INFO": {"GENE": "gene_name"},
    }


def test_from_lazyframe(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe."""
    tmp_file = tmp_path / "output.vcf"

    lf = polars.scan_parquet(DATA_DIR / "no_info.parquet")

    io.vcf.from_lazyframe(lf, tmp_file)

    assert filecmp.cmp(tmp_file, DATA_DIR / "no_info.parquet2vcf.vcf")

    io.vcf.from_lazyframe(lf, tmp_file, io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt", "vid", "vid"))

    assert filecmp.cmp(tmp_file, DATA_DIR / "no_info.parquet2vcf.vcf")
