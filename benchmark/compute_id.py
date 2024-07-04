"""Benchmark of variants id computation."""

# std import
from __future__ import annotations

import pathlib
import typing

# 3rd party import
import polars
import pytest

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest_benchmark

# project import
from variantplaner import ContigsLength, VcfHeader
from variantplaner.normalization import add_variant_id as __default_add_id

from benchmark import __generate_vcf

DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data"
chrom2length = ContigsLength()
chrom2length.from_path(DATA_DIR / "grch38.92.csv")


def __custom_vcf_parsing(input_path: pathlib.Path) -> polars.LazyFrame:
    """Custom local vcf reading in polars LazyFrame."""
    header = VcfHeader()
    header.from_files(input_path)

    col_name = {f"column_{i}": name for (i, name) in enumerate(header.column_name(), start=1)}

    lf = polars.scan_csv(
        input_path,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        schema_overrides={"column_1": polars.Utf8, "column_2": polars.UInt64},
        ignore_errors=True,
    )

    return lf.rename(col_name)


def __hash_add_id(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Add column id of variant by hashing."""
    return lf.with_columns(
        polars.concat_str(
            [
                polars.col("chr"),
                polars.lit("-"),
                polars.col("pos"),
                polars.lit("-"),
                polars.col("ref"),
                polars.lit("-"),
                polars.col("alt"),
            ],
        ).hash(),
    )


def __rust_add_id(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Add column id of variant by use rust code."""
    real_pos_max = chrom2length.lf.select([polars.col("length").sum()]).collect().get_column("length").max()

    lf = lf.join(chrom2length.lf, left_on="chr", right_on="contig", how="left")
    lf = lf.with_columns(real_pos=polars.col("pos") + polars.col("offset"))

    return lf.with_columns(
        id=polars.col("real_pos").variant_id.compute(  # type: ignore # noqa: PGH003
            polars.col("ref"),
            polars.col("alt"),
            real_pos_max,
        ),
    )


def __generate_id_hash(
    number_of_line: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="hash_id")
    def inner(
        tmp_path: pathlib.Path,
        benchmark: pytest_benchmark.BenchmarkSession,
    ) -> None:
        """Compute hash id of variant."""
        input_path = tmp_path / f"{number_of_line}.vcf"

        __generate_vcf(input_path, number_of_line)
        lf = __custom_vcf_parsing(input_path)

        benchmark(
            lambda: __hash_add_id(lf).collect(),
        )

    inner.__doc__ = f"""Compute hash id of {number_of_line} variant"""

    return inner


def __generate_id_rust(
    number_of_line: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="rust_id")
    def inner(
        tmp_path: pathlib.Path,
        benchmark: pytest_benchmark.BenchmarkSession,
    ) -> None:
        """Compute hash id of variant."""
        input_path = tmp_path / f"{number_of_line}.vcf"

        __generate_vcf(input_path, number_of_line)
        lf = __custom_vcf_parsing(input_path)

        benchmark(
            lambda: __rust_add_id(lf).collect(),
        )

    inner.__doc__ = f"""Compute rust id of {number_of_line} variant"""

    return inner


def __generate_id_default(
    number_of_line: int,
) -> typing.Callable[[pathlib.Path, pytest_benchmark.BenchmarkSession], None]:
    @pytest.mark.benchmark(group="default_id")
    def inner(
        tmp_path: pathlib.Path,
        benchmark: pytest_benchmark.BenchmarkSession,
    ) -> None:
        """Compute hash id of variant."""
        input_path = tmp_path / f"{number_of_line}.vcf"

        __generate_vcf(input_path, number_of_line)
        lf = __custom_vcf_parsing(input_path)

        benchmark(
            lambda: __default_add_id(lf, chrom2length.lf).collect(),
        )

    inner.__doc__ = f"""Compute default id of {number_of_line} variant"""

    return inner


for i in range(5, 14):
    number_of_line = 2**i
    globals()[f"add_id_hash_{number_of_line}"] = __generate_id_hash(number_of_line)
    globals()[f"add_id_rust_{number_of_line}"] = __generate_id_rust(number_of_line)
    globals()[f"add_id_default_{number_of_line}"] = __generate_id_default(number_of_line)
