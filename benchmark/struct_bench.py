"""Benchmark of struct."""

# std import
from __future__ import annotations

import pathlib
import typing
import random

# 3rd party import
import pytest
import polars


# project import

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark


DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"


def __generate_variant(common: int=10_000, diff: int=1000) -> polars.LazyFrame:
    """Generate a LazyFrame with many duplicate data and some not."""

    nuc = ["A", "C", "T", "G"]

    lf1 = polars.LazyFrame({
        "chr": random.choices(range(1, 25), k=common),
        "pos": random.choices(range(1, 2 ** 32), k=common),
        "ref": random.choices(nuc, k=common),
        "alt": random.choices(nuc, k=common),
    })

    lf2 = polars.LazyFrame({
        "chr": random.choices(range(1, 25), k=diff),
        "pos": random.choices(range(1, 2**32), k=diff),
        "ref": random.choices(nuc, k=diff),
        "alt": random.choices(nuc, k=diff),
    })

    lf3 = polars.LazyFrame({
        "chr": random.choices(range(1, 25), k=diff),
        "pos": random.choices(range(1, 2**32), k=diff),
        "ref": random.choices(nuc, k=diff),
        "alt": random.choices(nuc, k=diff),
    })

    return polars.concat([lf1, lf1, lf2, lf3])


def __generate_variant_merge(common: int, merge_type: str) -> typing.Callable[[pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group=f"merge_{merge_type}")
    def inner(benchmark: pytest_benchmark.BenchmarkSession) -> None:
        """Merge variant"""
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


for i in range(15, 25):
    common = 2 ** i
    globals()[f"variant_merge_id_{common}"] = __generate_variant_merge(common, "id")
    globals()[f"variant_merge_variant_{common}"] = __generate_variant_merge(common, "variant")
