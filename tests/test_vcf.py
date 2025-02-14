"""Test for the `object.vcf` module."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars

# project import
from variantplaner.objects import Vcf, VcfParsingBehavior

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_vcf_filter_star() -> None:
    """Filter vcf star."""
    vcf_path = DATA_DIR / "no_info.vcf"

    obj = Vcf()
    obj.from_path(vcf_path, DATA_DIR / "grch38.92.csv", VcfParsingBehavior.KEEP_STAR)

    assert obj.lf.filter(polars.col("alt") == "*").collect().height == 0


def test_chrom2length_overide_vcf_contig() -> None:
    """Check chromosome2length option overide vcf contig information."""
    vcf_path = DATA_DIR / "all_info.vcf"

    obj = Vcf()
    obj.from_path(vcf_path, None)

    assert obj.lf.null_count().collect().get_column("id").to_list() == [1]

    obj = Vcf()
    obj.from_path(vcf_path, DATA_DIR / "grch38.92.csv")

    assert obj.lf.null_count().collect().get_column("id").to_list() == [0]
