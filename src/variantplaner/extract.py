"""Extract information of [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) produce from raw vcf file parsing."""

# std import
from __future__ import annotations

import logging
import typing
from typing import Callable

# 3rd party import
import polars

# project import
from variantplaner.exception import NoGenotypeError

# typing import
if typing.TYPE_CHECKING:  # pragma: no cover
    import sys

    if sys.version_info >= (3, 11):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    T = typing.TypeVar("T")
    P = ParamSpec("P")

logger = logging.getLogger("manipulation")


def variants(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Extract variants only information of polars.LazyFrame.

    Args:
        lf: A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html)

    Returns:
        A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) with just variant information (id, chr, pos, ref, alt)
    """
    return lf.select(
        [
            polars.col("id"),
            polars.col("chr"),
            polars.col("pos"),
            polars.col("ref"),
            polars.col("alt"),
        ],
    )


def genotypes(
    lf: polars.LazyFrame,
    col2expr: dict[str, Callable[[polars.Expr, str], polars.Expr]],
    format_str: str = "GT:AD:DP:GQ",
) -> polars.LazyFrame:
    """Extract genotypes information of raw [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).

    Only line with format value match `format_str` are considered.

    Args:
        lf: The target [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html)
        col2expr: A dict associate column name and function to apply to create [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) column (produce by io.vcf.format2expr)
        format_str: Only variants match with this string format are considered

    Returns:
        A [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) with variant id, sample information and genotypes information

    Raises:
        NoGenotypeError: If none of the lf columns is equal to 'format'
    """
    if "format" not in lf.columns:
        raise NoGenotypeError

    lf = lf.select([*lf.columns[lf.columns.index("format") :]])

    # Clean bad variant
    lf = lf.filter(polars.col("format").str.starts_with(format_str)).select(*lf.columns[1:])

    # Found index of genotype value
    col_index = {
        key: index
        for (index, key) in enumerate(
            format_str.split(":"),
        )
    }

    # Pivot value
    genotypes = lf.melt(id_vars=["id"]).with_columns(
        [
            polars.col("id"),
            polars.col("variable").alias("sample"),
            polars.col("value").str.split(":"),
        ],
    )

    # Split genotype column in sub value
    genotypes = genotypes.with_columns(
        [polars.col("value").list.get(index).pipe(function=col2expr[col], name=col) for col, index in col_index.items()],  # type: ignore # noqa: PGH003
    )

    # Select intrusting column
    genotypes = genotypes.select(["id", "sample", *[col.lower() for col in col_index]])

    if "gt" in genotypes.columns:
        return genotypes.filter(polars.col("gt") != 0)

    return genotypes


def merge_variants_genotypes(
    variants_lf: polars.LazyFrame,
    genotypes_lf: polars.LazyFrame,
    sample_name: list[str],
) -> polars.LazyFrame:
    """Merge variants and genotypes [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).

    Args:
       variants_lf: lazyframe with variants, column: (id, chr, pos, ref, alt).
       genotypes_lf: lazyframe with genotypes, column: (id, sample, [genotype column]).

    Returns:
        A lazyframe with all data
    """
    for sample in sample_name:
        geno2sample = (
            genotypes_lf.filter(polars.col("sample") == sample)
            .rename(
                {col: f"{sample}_{col}" for col in genotypes_lf.columns[2:]},
            )
            .drop("sample")
        )
        variants_lf = variants_lf.join(geno2sample, on="id", how="outer_coalesce")

    return variants_lf
