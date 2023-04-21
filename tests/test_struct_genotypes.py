"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import os
import pathlib
import typing

# 3rd party import
import polars

# project import
from variantplanner import struct

DATA_DIR = pathlib.Path(__file__).parent / "data"


def __scantree(path: pathlib.Path) -> typing.Iterator[pathlib.Path]:
    """Recursively yield DirEntry objects for given directory."""
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from __scantree(pathlib.Path(entry.path))
        else:
            yield pathlib.Path(entry.path)


def test_parquet(tmp_path: pathlib.Path) -> None:
    """Check partition genotype parquet."""
    struct.genotypes.hive([DATA_DIR / "no_info.genotypes.parquet"], tmp_path)

    partition_paths = set(__scantree(tmp_path))

    assert partition_paths == {
        tmp_path / "gt=1/sample=sample_2/id_mod=47/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=15/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=28/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=27/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=40/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=44/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=49/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=8/0.parquet",
        tmp_path / "gt=1/sample=sample_2/id_mod=24/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=8/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=31/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=17/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=40/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=49/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=44/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=27/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=22/0.parquet",
        tmp_path / "gt=1/sample=sample_3/id_mod=9/0.parquet",
        tmp_path / "gt=2/sample=sample_1/id_mod=28/0.parquet",
        tmp_path / "gt=2/sample=sample_1/id_mod=8/0.parquet",
        tmp_path / "gt=2/sample=sample_1/id_mod=15/0.parquet",
        tmp_path / "gt=2/sample=sample_1/id_mod=24/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=30/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=34/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=47/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=28/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=29/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=23/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=15/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=49/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=24/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=8/0.parquet",
        tmp_path / "gt=2/sample=sample_3/id_mod=44/0.parquet",
    }

    value = polars.concat([polars.read_parquet(path) for path in partition_paths]).drop("id_mod")
    truth = polars.read_parquet(DATA_DIR / "no_info.genotypes.parquet")

    assert sorted(value.get_column("id").to_list()) == sorted(truth.get_column("id").to_list())
    assert sorted(value.get_column("sample").to_list()) == sorted(truth.get_column("sample").to_list())
    assert sorted(value.get_column("gt").to_list()) == sorted(truth.get_column("gt").to_list())
    assert sorted(value.get_column("ad").to_list()) == sorted(truth.get_column("ad").to_list())
    assert sorted(value.get_column("dp").to_list()) == sorted(truth.get_column("dp").to_list())
    assert sorted(value.get_column("gq").fill_null(0).to_list()) == sorted(
        truth.get_column("gq").fill_null(0).to_list(),
    )
