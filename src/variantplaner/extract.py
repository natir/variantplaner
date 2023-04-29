"""Function to manage variant data."""

# std import
from __future__ import annotations

import logging
import typing
from typing import Callable

# 3rd party import
import polars

# project import
from variantplaner.exception import NoGenotypeError
from variantplaner.io.vcf import MINIMAL_COL_NUMBER

# typing import
if typing.TYPE_CHECKING:  # pragma: no cover
    import sys

    if sys.version_info >= (3, 11):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    T = typing.TypeVar("T")
    P = ParamSpec("P")


LF_COL_NB_NO_GENOTYPE = MINIMAL_COL_NUMBER + 1

logger = logging.getLogger("manipulation")


def variants(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Extract variants only information of lazyframe.

    Args:
        lf: A lazyframe

    Returns:
        A lazyframe with just variant information (id, chr, pos, ref, alt)
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
    """Extract genotypes information in lazyframe.

    Only variant with format column like 'GT:AD:DP:GQ' are support.

    Args:
        lf: A lazyframe
        col2expr: A dict associate colum name and function to apply to create specific column
        format_str: Only variants match with this string format are considere

    Returns:
        A lazyframe with genotype information (id, sample, gt, ad, db, gq)

    Raises:
        NoGenotypeError: If number of column in lf indicate no format or genotype (9)
    """
    if len(lf.columns) <= LF_COL_NB_NO_GENOTYPE:
        raise NoGenotypeError

    # Select genotype columns (side effect last columns are: format, [genotypes,…] and id)
    lf = lf.select([*lf.columns[MINIMAL_COL_NUMBER:]])

    # Clean bad variant
    lf = lf.filter(polars.col("format").str.starts_with(format_str))

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
        [polars.col("value").arr.get(index).pipe(function=col2expr[col], name=col) for col, index in col_index.items()],  # type: ignore # noqa: PGH003
    )

    # Select intrusting column
    genotypes = genotypes.select(["id", "sample", *[col.lower() for col in col_index]])

    return genotypes.filter(polars.col("gt") != 0)


def merge_variants_genotypes(
    variants_lf: polars.LazyFrame,
    genotypes_lf: polars.LazyFrame,
    sample_name: list[str],
) -> polars.LazyFrame:
    """Merge variants and genotypes lazyframe.

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
        variants_lf = variants_lf.join(geno2sample, on="id", how="outer")

    return variants_lf
