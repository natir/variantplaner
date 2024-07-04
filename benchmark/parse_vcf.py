"""Benchmark of vcf parsing."""

# std import
from __future__ import annotations

import pathlib
import typing

# 3rd party import
import pytest

if typing.TYPE_CHECKING:  # pragma: no cover
    import polars
    import pytest_benchmark


# project import
from variantplaner import Vcf, VcfParsingBehavior

from benchmark import __generate_vcf

DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"


def __worker(input_path: pathlib.Path) -> polars.LazyFrame:
    """Benchmark worker."""
    vcf_df = Vcf()

    vcf_df.from_path(input_path, DATA_DIR / "grch38.92.csv", behavior=VcfParsingBehavior.MANAGE_SV)

    return vcf_df.variants()


def __generate_parse_vcf(
    number_of_line: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="vcf_parsing")
    def inner(
        tmp_path: pathlib.Path,
        benchmark: pytest_benchmark.BenchmarkSession,
    ) -> None:
        """Parsing a vcf."""
        input_path = tmp_path / f"{number_of_line}.vcf"

        __generate_vcf(input_path, number_of_line)

        benchmark(
            lambda: __worker(input_path).collect(),
        )

    inner.__doc__ = f"""Parsing a vcf of {number_of_line} variant"""

    return inner


for i in range(5, 15):
    number_of_line = 2**i
    globals()[f"parse_vcf_{number_of_line}"] = __generate_parse_vcf(number_of_line)
