"""Tests for the `io.ped` module."""

# std import
from __future__ import annotations

import filecmp

# 3rd party import
import pathlib

import polars
import polars.testing

# project import
from variantplaner import io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_into_lazyframe() -> None:
    """Check into_lazyframe."""
    ped_lf = io.ped.into_lazyframe(DATA_DIR / "sample.ped").collect()

    truth = polars.read_parquet(DATA_DIR / "sample.parquet")

    polars.testing.assert_frame_equal(truth, ped_lf)


def test_from_lazyframe(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe."""
    output_path = tmp_path / "sample.ped"
    ped_lf = polars.scan_parquet(DATA_DIR / "sample.parquet")

    io.ped.from_lazyframe(ped_lf, output_path)

    filecmp.cmp(output_path, DATA_DIR / "sample.ped")
