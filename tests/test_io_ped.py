"""Tests for the `io.ped` module."""

# std import
from __future__ import annotations

import filecmp

# 3rd party import
import pathlib

import polars
import polars.testing

# project import
from variantplaner import Pedigree

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_into_lazyframe() -> None:
    """Check into_lazyframe."""
    ped = Pedigree()
    ped.from_path(DATA_DIR / "sample.ped")

    truth = polars.read_parquet(DATA_DIR / "sample.parquet")

    polars.testing.assert_frame_equal(truth, ped.lf.collect())


def test_from_lazyframe(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe."""
    output_path = tmp_path / "sample.ped"
    ped = Pedigree()
    ped.lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    ped.to_path(output_path)

    filecmp.cmp(output_path, DATA_DIR / "sample.ped")
