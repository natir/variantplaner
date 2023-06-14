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

from benchmark import __generate_vcf

DATA_DIR = pathlib.Path(__file__).parent / "data"


@pytest.mark.benchmark(group="vcf2parquet")
def variants(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants and genotypes."""
    input_path = tmp_path / f"{10_000}.vcf"
    variants_path = tmp_path / "variants.parquet"

    __generate_vcf(input_path, 10_000)

    runner = CliRunner()

    benchmark(
        lambda: runner.invoke(
            cli.main,
            ["vcf2parquet", "-i", str(input_path), "-v", str(variants_path)],
        ),
    )


@pytest.mark.benchmark(group="vcf2parquet")
def variants_genotypes(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants and genotypes."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()

    benchmark(
        lambda: runner.invoke(
            cli.main,
            ["vcf2parquet", "-i", str(DATA_DIR / "all_info.vcf"), "-v", str(variants_path), "-g", str(genotypes_path)],
        ),
    )


@pytest.mark.benchmark(group="vcf2parquet")
def variants_genotypes_annotations(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark vcf2parquet variants, genotypes and annotations."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()

    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "vcf2parquet",
                "-i",
                str(DATA_DIR / "all_info.vcf"),
                "-v",
                str(variants_path),
                "-g",
                str(genotypes_path),
                "-a",
                str(annotations_path),
            ],
        ),
    )


@pytest.mark.benchmark(group="parquet2vcf")
def basic(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark parquet2vcf without other."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            ["parquet2vcf", "-i", str(DATA_DIR / "no_info.parquet"), "-o", str(variants_path)],
        ),
    )


@pytest.mark.benchmark(group="parquet2vcf")
def add_genotype(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark parquet2vcf without other."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    benchmark(
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


@pytest.mark.benchmark(group="struct_variants")
def struct_variants(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct.variants without other."""
    merge_path = tmp_path / "merge.parquet"

    runner = CliRunner()
    benchmark(
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


@pytest.mark.benchmark(group="struct_genotypes")
def struct_genotypes(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark struct.genotypes without other."""
    prefix_path = tmp_path / "hive"

    runner = CliRunner()
    benchmark(
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


ANNOTATIONS_COLUMNS = [
    "-i",
    "AF_ESP",
    "-i",
    " AF_EXAC",
    "-i",
    " AF_TGP",
    "-i",
    " ALLELEID",
    "-i",
    " CLNDN",
    "-i",
    " CLNDNINCL",
    "-i",
    " CLNDISDB",
    "-i",
    " CLNDISDBINCL",
    "-i",
    " CLNHGVS",
    "-i",
    " CLNREVSTAT",
    "-i",
    " CLNSIG",
    "-i",
    " CLNSIGCONF",
    "-i",
    " CLNSIGINCL",
    "-i",
    " CLNVC",
    "-i",
    " CLNVCSO",
    "-i",
    " CLNVI",
    "-i",
    " DBVARID",
    "-i",
    " GENEINFO",
    "-i",
    " MC",
    "-i",
    " ORIGIN",
    "-i",
    " RS",
]


@pytest.mark.benchmark(group="annotations_vcf")
def full(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS,
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def full_rename(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                "-r",
                "annot_id",
                *ANNOTATIONS_COLUMNS,
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def one(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:2],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def two(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:4],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def three(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:6],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def four(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:8],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def five(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:10],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def six(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:12],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def seven(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:14],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def height(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:16],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def nine(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:18],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def ten(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:20],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def eleven(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:22],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def twelve(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:24],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def thirteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:26],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def fourteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:28],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def fiveteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:30],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def sixteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:32],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def seventeen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:34],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def eighteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:36],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def nineteen(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:38],
            ],
        ),
    )


@pytest.mark.benchmark(group="annotations_vcf")
def twenty(
    tmp_path: pathlib.Path,
    benchmark: pytest_benchmark.BenchmarkSession,
) -> None:
    """Benchmark conversion of annotations vcf."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    benchmark(
        lambda: runner.invoke(
            cli.main,
            [
                "annotations",
                "-i",
                str(DATA_DIR / "no_genotypes.vcf"),
                "-o",
                str(annotations_path),
                "vcf",
                *ANNOTATIONS_COLUMNS[:40],
            ],
        ),
    )
