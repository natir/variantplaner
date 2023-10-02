"""Benchmark of struct."""

# std import
from __future__ import annotations

import os
import pathlib
import random
import typing

# 3rd party import
import polars
import pytest

# project import
from variantplaner import io, normalization, struct

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark


DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"
chrom2length = io.csv.chr2length_into_lazyframe(DATA_DIR / "grch38.92.csv")


def __generate_variant(common: int = 10_000, diff: int = 1000) -> polars.LazyFrame:
    """Generate a LazyFrame with many duplicate data and some not."""
    nuc = ["A", "C", "T", "G"]

    lf1 = polars.LazyFrame(
        {
            "chr": random.choices([f"chr{chr_name}" for chr_name in range(1, 25)], k=common),
            "pos": random.choices(range(1, 2**32), k=common),
            "ref": random.choices(nuc, k=common),
            "alt": random.choices(nuc, k=common),
        },
    )

    lf2 = polars.LazyFrame(
        {
            "chr": random.choices([f"chr{chr_name}" for chr_name in range(1, 25)], k=diff),
            "pos": random.choices(range(1, 2**32), k=diff),
            "ref": random.choices(nuc, k=diff),
            "alt": random.choices(nuc, k=diff),
        },
    )

    lf3 = polars.LazyFrame(
        {
            "chr": random.choices([f"chr{chr_name}" for chr_name in range(1, 25)], k=diff),
            "pos": random.choices(range(1, 2**32), k=diff),
            "ref": random.choices(nuc, k=diff),
            "alt": random.choices(nuc, k=diff),
        },
    )

    return normalization.add_variant_id(polars.concat([lf1, lf2, lf1, lf3, lf1]), chrom2length)


def __generate_variant_merge(
    common: int,
    merge_type: str,
) -> typing.Callable[[pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group=f"merge_{merge_type}")
    def inner(benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant."""
        lf = __generate_variant(common, common // 10)

        if merge_type == "id":
            benchmark(lambda: lf.unique(subset="id"))
        else:
            benchmark(lambda: lf.unique(subset=("chr", "pos", "ref", "alt")))

    if merge_type == "id":
        inner.__doc__ = f"""Merge variant on id with {common} variant and {common / 10 * 2} different."""
    else:
        inner.__doc__ = f"""Merge variant on variant with {common} variant and {common / 10 * 2} different."""

    return inner


def __generate_variant_merge_on_disk(
    threads: int,
    nb_file: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="merge_variant_on_disk")
    def inner(tmp_path: pathlib.Path, benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant on disk."""
        paths = []
        for i in range(nb_file):
            paths.append(tmp_path / f"{i}.parquet")
            __generate_variant().collect().write_parquet(tmp_path / f"{i}.parquet")

        os.environ["POLARS_MAX_THREADS"] = str(threads)
        benchmark(lambda: struct.variants.merge(paths, tmp_path / "output.parquet", memory_limit=5_000_000))

    return inner


for i in range(10, 20):
    common = 2**i
    globals()[f"variant_merge_id_{common}"] = __generate_variant_merge(common, "id")
    globals()[f"variant_merge_variant_{common}"] = __generate_variant_merge(common, "variant")


for i in range(10, 21, 2):
    globals()[f"variant_merge_on_disk_{i}"] = __generate_variant_merge_on_disk(8, i)
