"""Benchmark of struct."""

# std import
from __future__ import annotations

import pathlib
import typing

import polars

# 3rd party import
import pytest

# project import

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark


DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"


def __merge_id() -> None:
    """Merge polars.LazyFrame by id."""
    lfs = polars.concat(
        [
            polars.scan_parquet(DATA_DIR / "no_genotypes.variants.parquet"),
            polars.scan_parquet(DATA_DIR / "no_info.variants.parquet"),
        ],
    )

    lf = lfs.unique(subset="id")

    lf.collect()


def __merge_variant() -> None:
    """Merge polars.LazyFrame by id."""
    lfs = polars.concat(
        [
            polars.scan_parquet(DATA_DIR / "no_genotypes.variants.parquet"),
            polars.scan_parquet(DATA_DIR / "no_info.variants.parquet"),
        ],
    )

    lf = lfs.unique(subset=("chr", "pos", "ref", "alt"))

    lf.collect()


@pytest.mark.benchmark(group="merge")
def by_id(
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct variants merge."""
    benchmark(__merge_id)


@pytest.mark.benchmark(group="merge")
def by_variant(
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct variants merge."""
    benchmark(__merge_variant)
