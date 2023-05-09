"""Function to generate information."""

# std import
from __future__ import annotations

import logging

# 3rd party import
import polars

# project import

logger = logging.getLogger("generate")


def transmission_ped(
    genotypes_lf: polars.LazyFrame,
    pedigree_lf: polars.LazyFrame,
) -> polars.DataFrame:
    """Compute transmission of each variants.

    **Warning**: only the first sample with two parent are consider.

    Args:
        genotypes_lf: Genotypes LazyFrame.
        pedigree_lf: Pedigree LazyFrame.


    Returns:
         DataFrame with transmission information
    """
    pedigree_lf = pedigree_lf.filter(polars.col("father_id") != "unknow").filter(polars.col("mother_id") != "unknow")

    familly_info = pedigree_lf.collect().row(0, named=True)

    return transmission(genotypes_lf, familly_info["personal_id"], familly_info["mother_id"], familly_info["father_id"])


def transmission(
    genotypes_lf: polars.LazyFrame,
    index_name: str,
    mother_name: str,
    father_name: str,
) -> polars.DataFrame:
    """Compute transmission of each variants.

    Args:
        genotypes_lf: Genotypes LazyFrame.
        index_name: Sample name of index case
        mother_name: Sample name of mother
        father_name: Sample name of father

    Returns:
         DataFrame with transmission information.
         With gt, ad, dp and gq for index, mother and father.
         If any of them isn't present value are set to polars.Null (0 for gt)
         Columns transmission contains: index_gt * 100 + mother_gt * 10 + father_gt
    """
    group_lf = genotypes_lf.groupby("id").all().collect()
    group_lf = group_lf.filter(polars.col("sample").arr.contains(index_name))

    logger.debug(f"{group_lf.row(0)}")

    # I assume sample order is all time the same but I'm not sure
    sample2index = {sample: idx for idx, sample in enumerate(group_lf.row(0)[1], start=0)}

    logger.debug(f"{sample2index}")

    transmission_lf = group_lf.with_columns(
        polars.col("gt").arr.get(sample2index[index_name]).alias("index_gt"),
        polars.col("ad").arr.get(sample2index[index_name]).alias("index_ad"),
        polars.col("dp").arr.get(sample2index[index_name]).alias("index_dp"),
        polars.col("gq").arr.get(sample2index[index_name]).alias("index_gq"),
    )

    if mother_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            polars.col("gt").arr.get(sample2index[mother_name]).alias("mother_gt"),
            polars.col("ad").arr.get(sample2index[mother_name]).alias("mother_ad"),
            polars.col("dp").arr.get(sample2index[mother_name]).alias("mother_dp"),
            polars.col("gq").arr.get(sample2index[mother_name]).alias("mother_gq"),
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            polars.lit(0).alias("mother_gt"),
            polars.lit(None).alias("mother_ad"),
            polars.lit(None).alias("mother_dp"),
            polars.lit(None).alias("mother_gq"),
        )

    if father_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            polars.col("gt").arr.get(sample2index[father_name]).alias("father_gt"),
            polars.col("ad").arr.get(sample2index[father_name]).alias("father_ad"),
            polars.col("dp").arr.get(sample2index[father_name]).alias("father_dp"),
            polars.col("gq").arr.get(sample2index[father_name]).alias("father_gq"),
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            polars.lit(0).alias("father_gt"),
            polars.lit(None).alias("father_ad"),
            polars.lit(None).alias("father_dp"),
            polars.lit(None).alias("father_gq"),
        )

    transmission_lf = transmission_lf.with_columns(
        (polars.col("index_gt") * 100 + polars.col("mother_gt") * 10 + polars.col("father_gt")).alias("origin"),
    )

    transmission_lf = transmission_lf.drop(["sample", "gt", "ad", "dp", "gq"])

    return transmission_lf
