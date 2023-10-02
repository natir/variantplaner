"""Benchmark of hive partition."""

# std import
from __future__ import annotations

import pathlib
import random
import typing

# 3rd party import
import polars
import pytest

# project import
from variantplaner import struct

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark


DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"


def __generate_genotypes(number: int = 10_000) -> polars.LazyFrame:
    """Generate a LazyFrame with random genotype."""
    samples = [f"sample_{i}" for i in range(20)]

    return polars.LazyFrame(
        {
            "id": random.choices(range(2**32), k=number),
            "sample": random.choices(samples, k=number),
            "gt": random.choices([1, 2], k=number),
            "ad": [[random.choice(range(256)), random.choice(range(256))] for _ in range(number)],
            "dp": random.choices(range(256), k=number),
            "gq": random.choices(range(100), k=number),
        },
        schema={
            "id": polars.UInt64,
            "sample": polars.Utf8,
            "gt": polars.UInt8,
            "ad": polars.List(polars.UInt16),
            "dp": polars.UInt16,
            "gq": polars.UInt16,
        },
    )


def __generate_hive(
    genotype_number: int,
    group_size: int,
    number_of_file: int = 40,
) -> typing.Callable[[pytest_benchmark.BenchmarkSession], None]:
    """Generate benchmark function."""

    @pytest.mark.benchmark(group="hive_partitioning")
    def inner(tmp_path: pathlib.Path, benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Hive genotype."""
        random.seed(42)

        paths = []
        for i in range(number_of_file):
            file_path = tmp_path / "samples" / f"{i}.parquet"
            paths.append(file_path)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            __generate_genotypes(genotype_number).sink_parquet(file_path)

        benchmark(
            lambda: struct.genotypes.hive(
                paths,
                tmp_path / f"hive_{group_size}",
                threads=8,
                file_per_thread=group_size,
            ),
        )

    inner.__doc__ = """Perform hive partitioning on some genotypes file"""

    return inner


for i in range(2, 16, 2):
    globals()[f"hive_partitioning_group_size_{i}"] = __generate_hive(100, i, 16)
