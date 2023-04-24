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
from variantplanner import struct

DATA_DIR = pathlib.Path(__file__).parent / "data"

MERGE_IDS = [
    14242097786807219277,
    11061441912435479070,
    15556898136537930176,
    11711297912729185188,
    10807139693934085530,
    11754354483181031149,
    4797423062503056673,
    12330835685209145,
    3370214193463694380,
    13133999947024542022,
    800025099629731040,
    9722743600114191928,
    9554065504802534594,
    2295086632353399847,
    8068220590006330225,
    14452059200390301692,
    1062943517767600077,
    3925979271274072158,
    15044799663914522374,
    16588824010281864939,
    1485127797546295921,
    14498991407705406549,
    4455610085260365294,
    6760366560597639209,
    17246148306088235974,
    270423859963006067,
    11952270119407113547,
    10499890489289074188,
    11976298723311757531,
    721273245970054755,
    7645340256260622458,
    11033100074712141168,
    14319454638398535779,
    17224330959165379425,
    11920660801823307015,
    3585291275407140934,
    18172124183201916611,
    10212290916779149099,
    4759541468123016608,
    11650605831284591550,
    12945673213782943443,
    9664051898850350365,
    2373425740527121342,
    11903476464390241143,
    10487392163259126218,
    5356120651941363990,
    4235309537048834275,
    16247809233398557031,
]


def test_random_string() -> None:
    """Check random string."""
    random.seed(42)

    assert struct.variants.__random_string() == "OhbVrpoiVgRVIfLBcbfnoGMbJmTPSI"


def test_chunk_by_memory() -> None:
    """Check by memory."""
    paths = sorted(pathlib.Path(entry.path) for entry in os.scandir(DATA_DIR) if entry.is_file())

    chunks = list(struct.variants.__chunk_by_memory(paths, 10000))

    truth = [
        [
            DATA_DIR / "annotations.csv",
            DATA_DIR / "no_genotypes.parquet",
            DATA_DIR / "no_genotypes.variants.parquet",
            DATA_DIR / "no_genotypes.vcf",
        ],
        [
            DATA_DIR / "no_info.csv",
            DATA_DIR / "no_info.genotypes.parquet",
            DATA_DIR / "no_info.parquet",
            DATA_DIR / "no_info.parquet2vcf.vcf",
        ],
        [
            DATA_DIR / "no_info.tsv",
            DATA_DIR / "no_info.variants.parquet",
            DATA_DIR / "no_info.vcf",
            DATA_DIR / "only_header.vcf",
        ],
        [],
    ]

    assert chunks == truth


def test_concat_uniq_by_id(tmp_path: pathlib.Path) -> None:
    """Check concat_uniq_by_id."""
    tmp_file = tmp_path / "merge_by_id.parquet"

    struct.variants.concat_uniq_by_id(
        [DATA_DIR / "no_genotypes.variants.parquet", DATA_DIR / "no_info.variants.parquet"],
        tmp_file,
    )

    lf = polars.scan_parquet(tmp_file)

    assert lf.collect().get_column("id").to_list() == MERGE_IDS


def test_merge(tmp_path: pathlib.Path) -> None:
    """Check merge."""
    tmp_file = tmp_path / "merge_parquet.parquet"

    struct.variants.merge(
        [DATA_DIR / "no_genotypes.variants.parquet", DATA_DIR / "no_info.variants.parquet"],
        tmp_file,
    )

    lf = polars.scan_parquet(tmp_file)

    assert lf.collect().get_column("id").to_list() == MERGE_IDS
