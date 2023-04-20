"""Function to manage variant data."""

# std import
from __future__ import annotations

import logging

# 3rd party import
import polars

from variantplanner.exception import NoGenotypeError
from variantplanner.io.vcf import MINIMAL_COL_NUMBER

# project import


LF_COL_NB_NO_GENOTYPE = MINIMAL_COL_NUMBER + 1

logger = logging.getLogger("manipulation")


def extract_variants(lf: polars.LazyFrame) -> polars.LazyFrame:
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


def extract_genotypes(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Extarct genotypes information in lazyframe.

    Only variant with format column like 'GT:AD:DP:GQ' are support.

    Args:
        lf: A lazyframe

    Returns:
        A lazyframe with genotype information (id, sample, gt, ad, db, gq)

    Raises:
        NoGenotypeError: If number of column in lf indicate no format or genotype (9)
    """
    if len(lf.columns) <= LF_COL_NB_NO_GENOTYPE:
        raise NoGenotypeError

    # Select genotype columns (side effect last columns are format, [genotypes,…] and id)
    lf = lf.select([*lf.columns[MINIMAL_COL_NUMBER:]])

    # Clean bad variant
    lf = lf.filter(polars.col("format").str.starts_with("GT:AD:DP:GQ"))

    # Found index of genotype value
    col_index = {
        key: index
        for (index, key) in enumerate(
            lf.select(["format"]).first().collect()["format"][0].split(":"),
        )
    }

    # Pivot value
    genotypes = lf.melt(id_vars=["id"]).with_columns(
        [
            polars.col("id"),
            polars.col("variable").alias("sample"),
        ],
    )

    # Split genotype column in sub value
    genotypes = genotypes.select(
        [
            polars.col("id"),
            polars.col("sample"),
            # gt column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["GT"])
            .str.count_match("1")
            .cast(polars.UInt8)
            .alias("gt"),
            # ad column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["AD"])
            .str.split(",")
            .arr.eval(polars.element().str.parse_int(10, strict=False).cast(polars.UInt16))
            .alias("ad"),
            # dp column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["DP"])
            .str.parse_int(10, strict=False)
            .cast(polars.UInt16)
            .alias("dp"),
            # gq column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["GQ"])
            .str.parse_int(10, strict=False)
            .cast(polars.UInt16)
            .alias("gq"),
        ],
    )

    return genotypes.filter(polars.col("gt") != 0)
