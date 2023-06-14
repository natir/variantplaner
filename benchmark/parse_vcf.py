"""Benchmark of vcf parsing."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:
    import pathlib

# 3rd party import
import pytest

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark

# project import
from variantplaner import io

from benchmark import __generate_vcf


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
            lambda: io.vcf.into_lazyframe(input_path, extension=io.vcf.IntoLazyFrameExtension.MANAGE_SV).collect(),
        )

    inner.__doc__ = f"""Parsing a vcf of {number_of_line} variant"""

    return inner


for i in range(5, 15):
    number_of_line = 2**i
    globals()[f"parse_vcf_{number_of_line}"] = __generate_parse_vcf(number_of_line)
