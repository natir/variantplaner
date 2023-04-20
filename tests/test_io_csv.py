"""Tests for the `io.csv` module."""

# std import
from __future__ import annotations

import pathlib

import polars

# 3rd party import
# project import
from variantplanner import io, manipulation

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_into_lazyframe() -> None:
    """Check io.csv.into_lazyframe."""
    vcf = manipulation.extract_variants(io.vcf.into_lazyframe(DATA_DIR / "good.vcf")).collect()

    csv = (
        io.csv.into_lazyframe(
            DATA_DIR / "good.csv",
            "chr",
            "pos",
            "ref",
            "alt",
            ["id"],
            separator=",",
            dtypes={"id": polars.UInt64},
        )
        .select(["id", "chr", "pos", "ref", "alt"])
        .collect()
    )

    assert csv.frame_equal(vcf)
