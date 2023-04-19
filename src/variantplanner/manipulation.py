"""Function to manage variant data."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


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
    """
    # Select genotype columns
    lf = lf.select(["id", "format", *lf.columns[9:]])

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
    genotypes = (
        lf.drop("format")
        .melt(id_vars=["id"])
        .with_columns(
            [
                polars.col("id"),
                polars.col("variable").alias("sample"),
                polars.col("value").str.replace_all(r"\.", "0"),
            ],
        )
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
            .arr.eval(polars.element().str.parse_int(10).cast(polars.UInt16))
            .alias("ad"),
            # dp column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["DP"])
            .str.parse_int(10)
            .cast(polars.UInt16)
            .alias("dp"),
            # gq column
            polars.col("value")
            .str.split(":")
            .arr.get(col_index["GQ"])
            .str.parse_int(10)
            .cast(polars.UInt16)
            .alias("gq"),
        ],
    )

    return genotypes.filter(polars.col("gt") == 0)
