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
    21788637528068,
    22420668831105,
    22531532652550,
    28503684677645,
    216629694693385,
    2865209939989626902,
    2865214825649143842,
    2865214829944111138,
    2865214831957377026,
    2865214831957377025,
    3382250923559485444,
    3382273467842822145,
    3382274273283407906,
    3382286558903140374,
    3382286559037358105,
    6356557166639316998,
    6356557209588989953,
    6356558437949636610,
    6356559434382049293,
    6356561358527397901,
    6174051668938719235,
    6174052457065218057,
    6174053404105506822,
    6174053880846876681,
    6174054230886711305,
    148464268738563,
    149424193929223,
    149641089777676,
    149827920855049,
    150313252159501,
    1988454028148743,
    1988466913050636,
    1988475502985222,
    1988486240403458,
    1988503420272646,
    1988518452658185,
    1988527042592774,
    1988539927494662,
    1988576434716678,
    1988591467102222,
    1988600057036812,
    1988608646971404,
    1988623679356940,
    1988625826840588,
    1988627974324230,
    1988632269291526,
    1997451984633865,
    1997452387287654,
}


def test_random_string() -> None:
    """Check random string."""
    random.seed(42)

    assert struct.variants.__random_string() == "OhbVrpoiVgRVIfLBcbfnoGMbJmTPSI"


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
        [DATA_DIR / "no_genotypes.variants.parquet", DATA_DIR / "no_info.variants.parquet"],
        tmp_file,
    )

    lf = polars.scan_parquet(tmp_file)

    print(set(lf.collect().get_column("id").to_list()))
    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS


def test_merge(tmp_path: pathlib.Path) -> None:
    """Check merge."""
    tmp_file = tmp_path / "merge_parquet.parquet"

    os.environ["POLARS_MAX_THREADS"] = str(2)

    struct.variants.merge(
        [DATA_DIR / "no_genotypes.variants.parquet", DATA_DIR / "no_info.variants.parquet"],
        tmp_file,
    )

    lf = polars.scan_parquet(tmp_file)

    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS
