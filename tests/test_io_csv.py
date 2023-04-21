"""Tests for the `io.csv` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars
import polars.testing

# project import
from variantplanner import extract, io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_into_lazyframe() -> None:
    """Check io.csv.into_lazyframe."""
    vcf = extract.variants(io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf"))

    csv = io.csv.into_lazyframe(
        DATA_DIR / "no_info.csv",
        "chr",
        "pos",
        "ref",
        "alt",
        ["id"],
        separator=",",
        dtypes={"id": polars.UInt64},
    )

    polars.testing.assert_frame_equal(vcf, csv, check_column_order=False)
