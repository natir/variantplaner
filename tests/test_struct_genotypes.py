"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import os
import pathlib
import typing

# 3rd party import
import polars

try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:  # pragma: no cover
    pass
else:
    cleanup_on_sigterm()


# project import
from variantplaner import struct

DATA_DIR = pathlib.Path(__file__).parent / "data"


def __scantree(path: pathlib.Path) -> typing.Iterator[pathlib.Path]:
    """Recursively yield DirEntry objects for given directory."""
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from __scantree(pathlib.Path(entry.path))
        else:
            yield pathlib.Path(entry.path)


def test_manage_group(tmp_path: pathlib.Path) -> None:
    """Check manage group."""
    prefix = tmp_path / "manage_group"

    func = struct.genotypes.__manage_group(prefix, ["sample", "gt"], "1")

    df = polars.read_parquet(DATA_DIR / "no_info.genotypes.parquet")

    df.groupby("sample", "gt").apply(func)

    partition_paths = set(__scantree(prefix))

    assert partition_paths == {
        prefix / "sample=sample_3/gt=2/1.parquet",
        prefix / "sample=sample_2/gt=1/1.parquet",
        prefix / "sample=sample_3/gt=1/1.parquet",
        prefix / "sample=sample_1/gt=2/1.parquet",
    }


def test_write_or_add(tmp_path: pathlib.Path) -> None:
    """Check write or add."""
    path = tmp_path / "write_or_add.parquet"

    assert not path.exists()

    df = polars.read_parquet(DATA_DIR / "no_info.genotypes.parquet")

    struct.genotypes.__write_or_add(df, path)

    assert path.exists()

    struct.genotypes.__write_or_add(df, path)

    truth = polars.concat([df, df])

    reality = polars.read_parquet(path)

    polars.testing.assert_frame_equal(truth, reality)


def test_hive(tmp_path: pathlib.Path) -> None:
    """Check partition genotype parquet."""
    struct.genotypes.hive(
        [DATA_DIR / "no_info.genotypes.parquet", DATA_DIR / "no_info.genotypes.parquet"],
        tmp_path,
        threads=2,
        file_per_thread=1,
    )

    partition_paths = set(__scantree(tmp_path))

    value = polars.concat([polars.read_parquet(path) for path in partition_paths]).drop("id_mod")
    truth = polars.concat(
        [
            polars.read_parquet(DATA_DIR / "no_info.genotypes.parquet"),
            polars.read_parquet(DATA_DIR / "no_info.genotypes.parquet"),
        ],
    )

    assert sorted(value.get_column("id").to_list()) == sorted(truth.get_column("id").to_list())
    assert sorted(value.get_column("sample").to_list()) == sorted(truth.get_column("sample").to_list())
    assert sorted(value.get_column("gt").to_list()) == sorted(truth.get_column("gt").to_list())
    assert sorted(value.get_column("ad").to_list()) == sorted(truth.get_column("ad").to_list())
    assert sorted(value.get_column("dp").to_list()) == sorted(truth.get_column("dp").to_list())
    assert sorted(value.get_column("gq").fill_null(0).to_list()) == sorted(
        truth.get_column("gq").fill_null(0).to_list(),
    )
