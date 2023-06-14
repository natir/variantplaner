"""Benchmark of annotations extraction."""

# std import
from __future__ import annotations

import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib

# 3rd party import
import pytest

if typing.TYPE_CHECKING:  # pragma: no cover
    import polars
    import pytest_benchmark

# project import
from variantplaner import io

from benchmark import __generate_info, __generate_vcf


def __generate_annotations_extractions(
    number_of_col: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="annotations_extractions")
    def inner(
        tmp_path: pathlib.Path,
        benchmark: pytest_benchmark.BenchmarkSession,
    ) -> None:
        """Parsing a vcf."""
        input_path = tmp_path / f"{1_000}.vcf"

        info_names = sorted(__generate_info()[1].keys())
        __generate_vcf(input_path, 1_000)

        def worker() -> polars.DataFrame:
            vcf_header = io.vcf.extract_header(input_path)
            info_parser = io.vcf.info2expr(vcf_header, input_path, set(info_names[:number_of_col]))
            lf = io.vcf.into_lazyframe(input_path)

            lf = lf.with_columns(info_parser).drop(["chr", "pos", "ref", "alt", "filter", "qual", "info"])

            return lf.collect()

        benchmark(worker)

    inner.__doc__ = f"""Parsing a vcf of {number_of_col} variant"""

    return inner


for i in range(1, 28):
    globals()[f"annotations_extractions_{i}"] = __generate_annotations_extractions(i)
