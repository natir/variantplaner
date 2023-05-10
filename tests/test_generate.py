"""Tests for the `struct.variants` module."""

# std import
from __future__ import annotations

import filecmp
import pathlib

# 3rd party import
import polars
import polars.testing
import pytest

# project import
from variantplaner import exception, generate

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_transmission(tmp_path: pathlib.Path) -> None:
    """Check add_origin of genotype."""
    out_path = tmp_path / "transmission.parquet"

    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)
    transmission.write_parquet(out_path)

    truth = DATA_DIR / "transmission.parquet"

    filecmp.cmp(truth, out_path)


def test_transmission_missing_ad(tmp_path: pathlib.Path) -> None:
    """Check add_origin of genotype."""
    tmp_path / "transmission.parquet"

    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    genotypes_lf = genotypes_lf.drop("ad")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    assert transmission.columns == [
        "id",
        "index_gt",
        "index_dp",
        "index_gq",
        "mother_gt",
        "mother_dp",
        "mother_gq",
        "father_gt",
        "father_dp",
        "father_gq",
        "origin",
    ]


def test_transmission_nogt(tmp_path: pathlib.Path) -> None:
    """Check add_origin of genotype."""
    tmp_path / "transmission.parquet"

    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    genotypes_lf = genotypes_lf.drop("gt")

    with pytest.raises(exception.NoGTError):
        generate.transmission_ped(genotypes_lf, pedigree_lf)


def test_transmission_missing_mother() -> None:
    """Check add_origin of genotype."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_3")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    polars.Config.set_tbl_cols(500)

    assert transmission.get_column("mother_gt").to_list() == [0, 0, 0, 0, 0]
    assert transmission.get_column("mother_ad").to_list() == [None, None, None, None, None]
    assert transmission.get_column("mother_dp").to_list() == [None, None, None, None, None]
    assert transmission.get_column("mother_gq").to_list() == [None, None, None, None, None]


def test_transmission_missing_father() -> None:
    """Check add_origin of genotype."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_2")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    polars.Config.set_tbl_cols(500)

    assert transmission.get_column("father_gt").to_list() == [0, 0, 0, 0, 0]
    assert transmission.get_column("father_ad").to_list() == [None, None, None, None, None]
    assert transmission.get_column("father_dp").to_list() == [None, None, None, None, None]
    assert transmission.get_column("father_gq").to_list() == [None, None, None, None, None]
