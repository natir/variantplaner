"""Tests for the `struct.variants` module."""

# std import
from __future__ import annotations

import filecmp
import pathlib

# 3rd party import
import polars
import polars.testing

# project import
from variantplaner import generate

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
