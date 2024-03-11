"""Benchmark of vcf2parquet."""

# std import
from __future__ import annotations

import pathlib
import typing

# 3rd party import
import pytest
from click.testing import CliRunner

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark

# project import
from variantplaner import cli

from benchmark import __generate_format, __generate_vcf

DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"
NUMBER_OF_VARIANT = 1_000


@pytest.mark.benchmark(group="vcf2parquet")
def variants(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants and genotypes."""
    input_path = tmp_path / f"{NUMBER_OF_VARIANT}.vcf"
    variants_path = tmp_path / "variants.parquet"

    __generate_vcf(input_path, NUMBER_OF_VARIANT)

    runner = CliRunner()

    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "vcf2parquet",
                "-i",
                str(input_path),
                "-v",
                str(variants_path),
                "-c",
                str(DATA_DIR / "grch38.92.csv"),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="vcf2parquet")
def variants_annotations(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants and genotypes."""
    input_path = tmp_path / f"{NUMBER_OF_VARIANT}.vcf"

    variants_path = tmp_path / "variants.parquet"
    annotations_path = tmp_path / "annotations.parquet"

    __generate_vcf(input_path, NUMBER_OF_VARIANT)

    runner = CliRunner()

    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "vcf2parquet",
                "-i",
                str(input_path),
                "-v",
                str(variants_path),
                "-a",
                str(annotations_path),
                "-c",
                str(DATA_DIR / "grch38.92.csv"),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="vcf2parquet")
def variants_genotypes(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants and genotypes."""
    input_path = tmp_path / f"{NUMBER_OF_VARIANT}.vcf"

    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    __generate_vcf(input_path, NUMBER_OF_VARIANT)
    _, format_name2value = __generate_format()
    format_key = sorted(format_name2value.keys())

    runner = CliRunner()

    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "vcf2parquet",
                "-c",
                str(DATA_DIR / "grch38.92.csv"),
                "-i",
                str(input_path),
                "-v",
                str(variants_path),
                "-g",
                str(genotypes_path),
                "-f",
                ":".join(format_key),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="vcf2parquet")
def variants_genotypes_annotations(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants, genotypes and annotations."""
    input_path = tmp_path / f"{NUMBER_OF_VARIANT}.vcf"

    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"
    annotations_path = tmp_path / "annotations.parquet"

    __generate_vcf(input_path, NUMBER_OF_VARIANT)
    _, format_name2value = __generate_format()
    format_key = sorted(format_name2value.keys())

    runner = CliRunner()

    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "vcf2parquet",
                "-c",
                str(DATA_DIR / "grch38.92.csv"),
                "-i",
                str(input_path),
                "-v",
                str(variants_path),
                "-g",
                str(genotypes_path),
                "-a",
                str(annotations_path),
                "-f",
                ":".join(format_key),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="parquet2vcf")
def basic(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark parquet2vcf without other."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "parquet2vcf",
                "-i",
                str(DATA_DIR / "no_info.parquet"),
                "-o",
                str(variants_path),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="parquet2vcf")
def add_genotype(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark parquet2vcf without other."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "parquet2vcf",
                "-i",
                str(DATA_DIR / "no_info.variants.parquet"),
                "-o",
                str(variants_path),
                "-g",
                str(DATA_DIR / "no_info.genotypes.parquet"),
                "-F",
                "GT:AD:DP:GQ",
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="struct_variants")
def struct_variants(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct.variants without other."""
    merge_path = tmp_path / "merge.parquet"

    runner = CliRunner()
    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "struct",
                "-i",
                str(DATA_DIR / "no_genotypes.variants.parquet"),
                "-i",
                str(DATA_DIR / "no_info.variants.parquet"),
                "variants",
                "-o",
                str(merge_path),
            ],
        ),
    )

    assert result.exit_code == 0


@pytest.mark.benchmark(group="struct_genotypes")
def struct_genotypes(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct.genotypes without other."""
    prefix_path = tmp_path / "hive"

    runner = CliRunner()
    result = benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "struct",
                "-i",
                str(DATA_DIR / "no_info.genotypes.parquet"),
                "genotypes",
                "-p",
                str(prefix_path),
            ],
        ),
    )

    assert result.exit_code == 0
