"""Tests for the `struct.variants` module."""

# std import
from __future__ import annotations

import os
import pathlib
import random

# 3rd party import
import polars
import polars.testing

# project import
from variantplaner import struct

DATA_DIR = pathlib.Path(__file__).parent / "data"

MERGE_IDS = {
    148464134520839,
    149424059711495,
    149640955559940,
    149827786637317,
    150313117941765,
    1988453893931015,
    1988466778832900,
    1988475368767494,
    1988486106185734,
    1988503286054918,
    1988518318440453,
    1988526908375046,
    1988539793276934,
    1988576300498950,
    1988591332884486,
    1988599922819076,
    1988608512753668,
    1988623545139204,
    1988625692622852,
    1988627840106502,
    1988632135073798,
    1997451850416133,
    1997451850416146,
    216629560475653,
    21788369092616,
    22419729285149,
    22531398434822,
    28503550459909,
    2865209939855409174,
    2865214825380708362,
    2865214829675675658,
    2865214831823159301,
    2865214831823159302,
    3382250923425267716,
    3382273467708604421,
    3382274273014972426,
    3382286558768922633,
    3382286558768922646,
    6174051668804501511,
    6174052456931000325,
    6174053403971289094,
    6174053880712658949,
    6174054230752493573,
    6356557166505099270,
    6356557209454772229,
    6356558437815418886,
    6356559434247831557,
    6356561358393180165,
}


def test_chunk_by_memory() -> None:
    """Check by memory."""
    paths = [
        DATA_DIR / "no_info.csv",
        DATA_DIR / "no_info.genotypes.parquet",
        DATA_DIR / "no_info.parquet",
        DATA_DIR / "no_info.parquet2vcf.vcf",
        DATA_DIR / "no_info.parquet2vcf_genotypes.vcf",
        DATA_DIR / "no_info.tsv",
        DATA_DIR / "no_info.variants.parquet",
        DATA_DIR / "no_info.vcf",
        DATA_DIR / "no_info.vcf2parquet2vcf.vcf",
    ]

    chunks = list(struct.variants.__chunk_by_memory(paths, 10000))

    truth = [
        [
            DATA_DIR / "no_info.csv",
            DATA_DIR / "no_info.genotypes.parquet",
            DATA_DIR / "no_info.parquet",
            DATA_DIR / "no_info.parquet2vcf.vcf",
        ],
        [
            DATA_DIR / "no_info.parquet2vcf_genotypes.vcf",
            DATA_DIR / "no_info.tsv",
            DATA_DIR / "no_info.variants.parquet",
            DATA_DIR / "no_info.vcf",
        ],
        [DATA_DIR / "no_info.vcf2parquet2vcf.vcf"],
    ]

    assert chunks == truth


def test_concat_uniq(tmp_path: pathlib.Path) -> None:
    """Check concat_uniq_by_id."""
    tmp_file = tmp_path / "merge_by_id.parquet"

    struct.variants.__concat_uniq(
        [
            DATA_DIR / "no_genotypes.variants.parquet",
            DATA_DIR / "no_info.variants.parquet",
        ],
        tmp_file,
    )

    lf = polars.scan_parquet(tmp_file)

    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS


def test_merge(tmp_path: pathlib.Path) -> None:
    """Check merge."""
    tmp_file = tmp_path / "merge_parquet.parquet"

    os.environ["POLARS_MAX_THREADS"] = str(2)

    struct.variants.merge(
        [
            DATA_DIR / "no_genotypes.variants.parquet",
            DATA_DIR / "no_info.variants.parquet",
        ],
        tmp_file,
        False,
    )

    lf = polars.scan_parquet(tmp_file)

    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS

def test_merge_append(tmp_path: pathlib.Path) -> None:
    """Check merge append."""

    tmp_file = tmp_path / "merge_parquet.parquet"

    os.environ["POLARS_MAX_THREADS"] = str(2)

    struct.variants.merge(
        [
            DATA_DIR / "no_genotypes.variants.parquet",
            DATA_DIR / "no_info.variants.parquet",
        ],
        tmp_file,
        False,
    )

    lf = polars.scan_parquet(tmp_file)

    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS

    struct.variants.merge(
        [
            DATA_DIR / "sv.variants.parquet",
        ],
        tmp_file,
        True,
    )

    lf_sv = polars.scan_parquet(DATA_DIR / "sv.variants.parquet")
    sv_merge = MERGE_IDS | set(lf_sv.collect().get_column("id").to_list())

    lf = polars.scan_parquet(tmp_file)

    assert set(lf.collect().get_column("id").to_list()) == sv_merge
