"""Tests for the `struct.variants` module."""

# std import
from __future__ import annotations

import filecmp
import pathlib
import sys

# 3rd party import
import polars
import polars.testing
import pytest

# project import
from variantplaner import exception, generate

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_transmission(tmp_path: pathlib.Path) -> None:
    """Check transmission computation."""
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
    """Check transmission computation no gt."""
    tmp_path / "transmission.parquet"

    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    genotypes_lf = genotypes_lf.drop("gt")

    with pytest.raises(exception.NoGTError):
        generate.transmission_ped(genotypes_lf, pedigree_lf)


def test_transmission_missing_mother() -> None:
    """Check transmission computation no mother."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_3")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    polars.Config.set_tbl_cols(500)

    assert transmission.get_column("mother_gt").to_list() == [3, 3, 3, 3, 3]
    assert transmission.get_column("mother_ad").to_list() == [None, None, None, None, None]
    assert transmission.get_column("mother_dp").to_list() == [None, None, None, None, None]
    assert transmission.get_column("mother_gq").to_list() == [None, None, None, None, None]


def test_transmission_missing_father() -> None:
    """Check transmission computation no father."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_2")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    polars.Config.set_tbl_cols(500)

    assert transmission.get_column("father_gt").to_list() == [3, 3, 3, 3, 3]
    assert transmission.get_column("father_ad").to_list() == [None, None, None, None, None]
    assert transmission.get_column("father_dp").to_list() == [None, None, None, None, None]
    assert transmission.get_column("father_gq").to_list() == [None, None, None, None, None]


def test_transmission_all_mother_gt_null() -> None:
    """Check transmission computation mother gt are null."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.with_columns(
        polars.when(polars.col("sample") == "sample_3")
        .then(
            polars.lit(None),
        )
        .otherwise(polars.col("gt"))
        .alias("gt"),
    )

    polars.Config.set_tbl_rows(sys.maxsize)

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf).sort(by="id")

    assert transmission.get_column("mother_gt").to_list() == [0, 0, 0, 0, 0]
    assert transmission.get_column("mother_ad").to_list() == [[1, 4092], [5, 4828], [0, 4560], [0, 3084], [2, 4734]]
    assert transmission.get_column("mother_dp").to_list() == [4093, 4833, 4560, 3084, 4736]
    assert transmission.get_column("mother_gq").to_list() == [99, 99, 99, 99, 99]
