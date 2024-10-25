"""Read and write vcf file."""

# ruff: noqa: PLR5501

# std import
from __future__ import annotations

import logging
import typing

# 3rd party import
import polars

# project import


# type checking block
if typing.TYPE_CHECKING:  # pragma: no cover
    import pathlib
    import sys

    if sys.version_info >= (3, 11):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    T = typing.TypeVar("T")
    P = ParamSpec("P")

    from variantplaner import VcfHeader


MINIMAL_COL_NUMBER: int = 8
SAMPLE_COL_BEGIN: int = 9

logger = logging.getLogger("io.vcf")


if typing.TYPE_CHECKING:  # pragma: no cover
    RenameCol = typing.TypedDict(
        "RenameCol",
        {
            "#CHROM": str,
            "POS": str,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": list[tuple[str, str]],
            "FORMAT": str,
            "sample": dict[str, dict[str, str]],
        },
    )

DEFAULT_RENAME: RenameCol = {
    "#CHROM": "chr",
    "POS": "pos",
    "ID": "id",
    "REF": "ref",
    "ALT": "alt",
    "QUAL": ".",
    "FILTER": ".",
    "INFO": [],
    "FORMAT": "",
    "sample": {},
}


def build_rename_column(
    chromosome: str,
    pos: str,
    identifier: str,
    ref: str,
    alt: str,
    qual: str | None = ".",
    filter_col: str | None = ".",
    info: list[tuple[str, str]] | None = None,
    format_string: str | None = None,
    sample: dict[str, dict[str, str]] | None = None,
) -> RenameCol:
    """A helper function to generate rename column dict for [variantplaner.io.vcf.lazyframe_in_vcf][] function parameter.

    Returns:
        A rename column dictionary.
    """
    return {
        "#CHROM": chromosome,
        "POS": pos,
        "ID": identifier,
        "REF": ref,
        "ALT": alt,
        "QUAL": "." if qual is None else qual,
        "FILTER": "." if filter_col is None else filter_col,
        "INFO": [] if info is None else info,
        "FORMAT": "" if format_string is None else format_string,
        "sample": {} if sample is None else sample,
    }


def lazyframe_in_vcf(
    lf: polars.LazyFrame,
    output_path: pathlib.Path,
    /,
    vcf_header: VcfHeader | None = None,
    renaming: RenameCol = DEFAULT_RENAME,
) -> None:
    """Write [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) in vcf format.

    Warning: This function performs [polars.LazyFrame.collect][] before write vcf, this can have a significant impact on memory usage.

    Args:
        lf: LazyFrame contains information.
        output_path: Path to where vcf to write.

    Returns:
        None
    """
    select_column: list[str] = []

    lf = lf.with_columns(
        [
            polars.col(renaming["#CHROM"]).alias("#CHROM"),
            polars.col(renaming["POS"]).alias("POS"),
            polars.col(renaming["ID"]).alias("ID"),
            polars.col(renaming["REF"]).alias("REF"),
            polars.col(renaming["ALT"]).alias("ALT"),
        ],
    )

    select_column.extend(["#CHROM", "POS", "ID", "REF", "ALT"])

    if vcf_header is None:
        header = __generate_header(lf, renaming["INFO"], list(renaming["sample"].keys()), renaming["FORMAT"])
    else:
        header = "\n".join(vcf_header._header)

    if renaming["QUAL"] != ".":
        lf = lf.with_columns([polars.col(renaming["QUAL"]).alias("QUAL")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("QUAL")])

    select_column.append("QUAL")

    if renaming["FILTER"] != ".":
        lf = lf.with_columns([polars.col(renaming["FILTER"]).alias("FILTER")])
    else:
        lf = lf.with_columns([polars.lit(".").alias("FILTER")])

    select_column.append("FILTER")

    lf = (
        __rebuild_info_column(lf, renaming["INFO"])
        if renaming["INFO"]
        else lf.with_columns(polars.lit(".").alias("INFO"))
    )

    select_column.append("INFO")

    if renaming["FORMAT"]:
        lf = lf.with_columns(polars.lit(renaming["FORMAT"]).alias("FORMAT"))
        select_column.append("FORMAT")

    if renaming["FORMAT"] and renaming["sample"]:
        schema = lf.collect_schema()
        for sample_name in renaming["sample"]:
            lf = lf.with_columns(
                [
                    __lazy2format(
                        sample_name,
                        renaming["FORMAT"],
                        dict(zip(schema.names(), schema.dtypes())),
                    ).alias(sample_name),
                ],
            )
            select_column.append(sample_name)

    lf = lf.select([polars.col(col) for col in select_column])

    with open(output_path, "wb") as fh:
        fh.write(header.encode())
        fh.write(lf.collect().write_csv(separator="\t").encode())


def __rebuild_info_column(lf: polars.LazyFrame, vcfinfo2parquet_name: list[tuple[str, str]]) -> polars.LazyFrame:
    """Construct an INFO column from multiple columns of lf.

    Useful when you want serialise [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) in vcf file format.

    Args:
        lf: A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).
        vcfinfo2parquet_name: List of vcf column name and [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) column name.

    Returns:
        [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) with INFO column and remove lf column use.
    """
    lf = lf.with_columns(
        [
            polars.col(name).list.join(",").fill_null(".").alias(name)
            for name, dtype in zip(lf.columns, lf.dtypes)
            if isinstance(dtype, polars.List)
        ],
    )

    lf = lf.with_columns(
        [
            polars.col(name).cast(str).fill_null(".").alias(name)
            for name, dtype in zip(lf.columns, lf.dtypes)
            if not isinstance(dtype, polars.List)
        ],
    )

    lf = lf.with_columns(
        [
            polars.concat_str(
                [
                    polars.concat_str(
                        [
                            polars.lit(vcf_name),
                            polars.lit("="),
                            polars.col(parquet_name),
                        ],
                    )
                    for vcf_name, parquet_name in vcfinfo2parquet_name
                ],
                separator=";",
            ).alias("INFO"),
        ],
    )

    return lf.drop([p for (v, p) in vcfinfo2parquet_name])


def __generate_header(
    lf: polars.LazyFrame,
    vcfinfo2parquet_name: list[tuple[str, str]] | None = None,
    samples: list[str] | None = None,
    format_string: str | None = None,
) -> str:
    """Generate header of vcf file.

    Args:
        lf: Dataframe
        vcfinfo2parquet_name: vcf column name link to column name

    Returns:
        The header of vcf file
    """
    header = """##fileformat=VCFv4.3
##source=VariantPlanner
"""
    schema = lf.collect_schema()
    col2type = dict(zip(schema.names(), schema.dtypes()))

    if vcfinfo2parquet_name:
        for vcf_name, col_name in vcfinfo2parquet_name:
            number = "." if isinstance(col2type[col_name], polars.List) else "1"
            type_ = "String"
            type_ = "Float" if isinstance(col2type[col_name], polars.Float64) else type_
            type_ = "Integer" if isinstance(col2type[col_name], polars.Int64) else type_

            header += f'##INFO=<ID={vcf_name},Number={number},Type={type_},Description="Unknow">\n'

    if samples and format_string:
        sample = samples[0]
        for col in format_string.split(":"):
            col_name = f"{sample}_{col.lower()}"

            number = "." if isinstance(col2type[col_name], polars.List) else "1"
            type_ = "String"
            type_ = "Float" if isinstance(col2type[col_name], polars.Float64) else type_
            type_ = "Integer" if isinstance(col2type[col_name], polars.Int64) else type_

            header += f'##FORMAT=<ID={col.upper()},Number={number},Type={type_},Description="Unknow">\n'

    return header


def __lazy2format(sample_name: str, format_string: str, col2type: dict[str, polars.PolarsDataType]) -> polars.Expr:
    """Use sample name and column name to build sample column."""
    expression = []

    for col_name in format_string.split(":"):
        lazy_name = f"{sample_name}_{col_name.lower()}"
        if col_name == "GT":
            expression.append(
                polars.col(lazy_name)
                .cast(polars.Utf8)
                .fill_null("./.")
                .str.replace("1", "0/1")
                .str.replace("2", "1/1"),
            )
        elif isinstance(col2type[lazy_name], polars.List):
            expression.append(
                polars.col(lazy_name).cast(polars.List(polars.Utf8)).fill_null(["."]).list.join(","),
            )
        else:
            expression.append(
                polars.col(lazy_name).cast(polars.Utf8).fill_null("."),
            )

    return polars.concat_str(expression, separator=":")
