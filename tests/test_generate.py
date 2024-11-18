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
from variantplaner import Pedigree, exception, generate

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


def test_transmissions_with_ped() -> None:
    """Check transmission computation with a ped file."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree = Pedigree()
    pedigree.from_path(DATA_DIR / "missing_father.ped")

    transmission = generate.transmission_ped(genotypes_lf, pedigree.lf)

    assert transmission.get_column("index_gt").to_list() == [2, 2, 2, 2, 2]

    assert transmission.get_column("mother_gt").to_list() == [1, 1, 1, 1, 1]
    assert transmission.get_column("father_gt").to_list() == [None, None, None, None, None]

    assert transmission.get_column("origin").to_list() == ['#"~', '#"~', '#"~', '#"~', '#"~']


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


def test_transmission_mother_no_variant() -> None:
    """Check transmission computation no mother."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_3")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    assert transmission.get_column("mother_gt").to_list() == [
        0,
        0,
        0,
        0,
        0,
    ]
    assert transmission.get_column("mother_ad").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("mother_dp").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("mother_gq").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]

    assert transmission.get_column("origin").to_list() == ['#!"', '#!"', '#!"', '#!"', '#!"']


def test_transmission_missing_mother() -> None:
    """Check transmission computation no mother."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")
    pedigree_lf = pedigree_lf.with_columns(mother_id=polars.lit(None))

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    assert transmission.get_column("mother_gt").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("mother_ad").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("mother_dp").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("mother_gq").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]

    assert transmission.get_column("origin").to_list() == ['#~"', '#~"', '#~"', '#~"', '#~"']


def test_transmission_nogt(tmp_path: pathlib.Path) -> None:
    """Check transmission computation no gt."""
    tmp_path / "transmission.parquet"

    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    genotypes_lf = genotypes_lf.drop("gt")

    with pytest.raises(exception.NoGTError):
        generate.transmission_ped(genotypes_lf, pedigree_lf)


def test_transmission_father_no_variant() -> None:
    """Check transmission computation no father."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
    genotypes_lf = genotypes_lf.filter(polars.col("sample") != "sample_2")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    assert transmission.get_column("father_gt").to_list() == [
        0,
        0,
        0,
        0,
        0,
    ]
    assert transmission.get_column("father_ad").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("father_dp").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("father_gq").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("origin").to_list() == [
        "##!",
        "##!",
        "##!",
        "##!",
        "##!",
    ]


def test_transmission_missing_father() -> None:
    """Check transmission computation no father."""
    genotypes_lf = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")

    pedigree_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")
    pedigree_lf = pedigree_lf.with_columns(father_id=polars.lit(None))

    transmission = generate.transmission_ped(genotypes_lf, pedigree_lf)

    assert transmission.get_column("father_gt").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("father_ad").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("father_dp").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("father_gq").to_list() == [
        None,
        None,
        None,
        None,
        None,
    ]
    assert transmission.get_column("origin").to_list() == [
        "##~",
        "##~",
        "##~",
        "##~",
        "##~",
    ]


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

    assert transmission.get_column("mother_gt").to_list() == [
        0,
        0,
        0,
        0,
        0,
    ]
    assert transmission.get_column("mother_ad").to_list() == [
        [0, 4560],
        [2, 4734],
        [5, 4828],
        [1, 4092],
        [0, 3084],
    ]
    assert transmission.get_column("mother_dp").to_list() == [
        4560,
        4736,
        4833,
        4093,
        3084,
    ]
    assert transmission.get_column("mother_gq").to_list() == [99, 99, 99, 99, 99]


def test_31c177() -> None:
    """Check missing father information in 31C177."""
    pedigree = Pedigree()
    pedigree.from_path(DATA_DIR / "31C177.ped")

    genotypes_lf = polars.scan_parquet(DATA_DIR / "only_genotype_31C177.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree.lf)

    assert transmission.get_column("index_gt").to_list() == [1]
    assert transmission.get_column("father_gt").to_list() == [1]


def test_one_line_ped() -> None:
    """Check behavior if ped file as one line."""
    pedigree = Pedigree()
    pedigree.from_path(DATA_DIR / "one_line.ped")

    genotypes_lf = polars.scan_parquet(DATA_DIR / "only_genotype_31C177.parquet")

    transmission = generate.transmission_ped(genotypes_lf, pedigree.lf)

    assert transmission.get_column("index_gt").to_list() == [1]
    assert transmission.get_column("father_gt").to_list() == [None]
    assert transmission.get_column("mother_gt").to_list() == [None]
