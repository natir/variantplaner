"""Tests for the `struct.variants` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars
import polars.testing

# project import
from variantplaner import generate

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_origin() -> None:
    """Check add_origin of genotype."""
    genotypes = polars.scan_parquet("tests/data/no_info.genotypes.parquet")

    with_origin = generate.origin(genotypes, ["sample_1"], father="sample_2", mother="sample_3")
    truth = [
        4,
        4,
        4,
        4,
        4,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    assert truth == with_origin.collect().get_column("origin").to_list()

    with_origin = generate.origin(genotypes, ["sample_2"], father="sample_1", mother="sample_3")
    truth = [
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        4,
        4,
        4,
        4,
        4,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    assert truth == with_origin.collect().get_column("origin").to_list()

    with_origin = generate.origin(genotypes, ["sample_2", "sample_3"], mother="sample_3")
    truth = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
    ]

    assert truth == with_origin.collect().get_column("origin").to_list()

    with_origin = generate.origin(genotypes, ["sample_2", "sample_3"], father="sample_3")
    truth = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
    ]

    assert truth == with_origin.collect().get_column("origin").to_list()
