"""Function to generate information."""

# std import
from __future__ import annotations

# 3rd party import
import polars

# project import


def origin(
    lf: polars.LazyFrame,
    childs: list[str],
    mother: str | None = None,
    father: str | None = None,
) -> polars.LazyFrame:
    """Add origin of variants in genotypes LazyFrame.

    Only origin of probands variants are compute. Other variants origin are mark as Unknow

    Origin could be:
      - Unknow  -> 0
      - De-novo -> 1
      - Mother  -> 2
      - Father  -> 3
      - Both    -> 4

    Args:
        lf: With column id and sample.
        childs: List of sample of child
        father: Name of sample of probands father
        mother: Name of sample of probands mother

    Returns:
        Same LazyFrame as input with a new column origin
    """
    for child in childs:
        variant_in = lf.groupby("id").agg(polars.col("sample")).filter(polars.col("sample").arr.contains(child))

        if father is not None and mother is not None:
            both_lf = (
                variant_in.filter(polars.col("sample").arr.contains(father) & polars.col("sample").arr.contains(mother))
                .with_columns(polars.lit(4).cast(polars.UInt8).alias("origin"))
                .drop("sample")
            )
        else:
            both_lf = polars.LazyFrame(schema={"id": polars.UInt64, "origin": polars.UInt8})

        if mother is not None:
            mother_lf = (
                variant_in.filter(polars.col("sample").arr.contains(mother))
                .with_columns(polars.lit(2).cast(polars.UInt8).alias("origin"))
                .drop("sample")
            )
        else:
            mother_lf = polars.LazyFrame(schema={"id": polars.UInt64, "origin": polars.UInt8})

        if father is not None:
            father_lf = (
                variant_in.filter(polars.col("sample").arr.contains(father))
                .with_columns(polars.lit(3).cast(polars.UInt8).alias("origin"))
                .drop("sample")
            )
        else:
            father_lf = polars.LazyFrame(schema={"id": polars.UInt64, "origin": polars.UInt8})

        denovo_lf = (
            variant_in.filter(polars.col("sample").arr.lengths() == 1 & polars.col("sample").arr.contains(child))
            .with_columns(polars.lit(3).cast(polars.UInt8).alias("origin"))
            .drop("sample")
        )

        concat_lf = polars.concat([both_lf, mother_lf, father_lf, denovo_lf])

        concat_lf = concat_lf.unique(keep="first", subset=["id"]).with_columns(polars.lit(child).alias("sample"))

        lf = concat_lf.join(lf, on=["id", "sample"], how="outer")

    lf = lf.with_columns(
        [
            polars.col("origin").fill_null(0),
        ],
    )

    return lf
