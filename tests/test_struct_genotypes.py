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
    )

    partition_paths = set(__scantree(tmp_path))

    assert partition_paths == {
        tmp_path / "id_mod=56/2.parquet",
        tmp_path / "id_mod=194/1.parquet",
        tmp_path / "id_mod=43/2.parquet",
        tmp_path / "id_mod=224/2.parquet",
        tmp_path / "id_mod=194/2.parquet",
        tmp_path / "id_mod=6/1.parquet",
        tmp_path / "id_mod=198/1.parquet",
        tmp_path / "id_mod=6/2.parquet",
        tmp_path / "id_mod=198/2.parquet",
        tmp_path / "id_mod=160/1.parquet",
        tmp_path / "id_mod=75/1.parquet",
        tmp_path / "id_mod=160/2.parquet",
        tmp_path / "id_mod=115/1.parquet",
        tmp_path / "id_mod=94/1.parquet",
        tmp_path / "id_mod=75/2.parquet",
        tmp_path / "id_mod=33/1.parquet",
        tmp_path / "id_mod=238/1.parquet",
        tmp_path / "id_mod=115/2.parquet",
        tmp_path / "id_mod=94/2.parquet",
        tmp_path / "id_mod=29/1.parquet",
        tmp_path / "id_mod=33/2.parquet",
        tmp_path / "id_mod=238/2.parquet",
        tmp_path / "id_mod=237/1.parquet",
        tmp_path / "id_mod=103/1.parquet",
        tmp_path / "id_mod=29/2.parquet",
        tmp_path / "id_mod=70/1.parquet",
        tmp_path / "id_mod=99/1.parquet",
        tmp_path / "id_mod=237/2.parquet",
        tmp_path / "id_mod=103/2.parquet",
        tmp_path / "id_mod=85/1.parquet",
        tmp_path / "id_mod=70/2.parquet",
        tmp_path / "id_mod=99/2.parquet",
        tmp_path / "id_mod=41/1.parquet",
        tmp_path / "id_mod=85/2.parquet",
        tmp_path / "id_mod=41/2.parquet",
        tmp_path / "id_mod=122/1.parquet",
        tmp_path / "id_mod=205/1.parquet",
        tmp_path / "id_mod=44/1.parquet",
        tmp_path / "id_mod=122/2.parquet",
        tmp_path / "id_mod=44/2.parquet",
        tmp_path / "id_mod=56/1.parquet",
        tmp_path / "id_mod=205/2.parquet",
        tmp_path / "id_mod=43/1.parquet",
        tmp_path / "id_mod=224/1.parquet",
    }

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
