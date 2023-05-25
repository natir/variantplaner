"""Function to generate information."""

# std import
from __future__ import annotations

import logging

# 3rd party import
import polars

# project import
from variantplaner.exception import NoGTError

logger = logging.getLogger("generate")


def transmission_ped(
    genotypes_lf: polars.LazyFrame,
    pedigree_lf: polars.LazyFrame,
) -> polars.DataFrame:
    """Compute transmission of each variants.

    **Warning**: only the first sample with two parent are consider.

    Args:
        genotypes_lf: Genotypes LazyFrame, `gt` column are required.
        pedigree_lf: Pedigree LazyFrame.

    Returns:
         DataFrame with transmission information

    Raises:
        NoGTError: if genotypes_lf not containts gt column.
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
    """Compute how each variant are transmite to index case.

    Args:
        genotypes_lf: Genotypes polars.LazyFrame, `gt` column are required.
        index_name: Sample name of index case.
        mother_name: Sample name of mother.
        father_name: Sample name of father.

    Returns:
         DataFrame with transmission information.
         With genotyping information for index, mother and father.
         If any of them isn't present value are set to polars.Null (3 for gt)
         Columns transmission contains: index_gt * 100 + mother_gt * 10 + father_gt.
         Transmission: 230 mean homozygote variant not present in father but with no information about mother

    Raises:
        NoGTError: if genotypes_lf not containts gt column.
    """
    genotypes_column = list(genotypes_lf.columns[2:])
    if "gt" not in genotypes_column:
        raise NoGTError("genotype polars.LazyFrame")

    group_lf = genotypes_lf.groupby("id").all().collect()
    group_lf = group_lf.filter(polars.col("sample").arr.contains(index_name))

    logger.debug(f"{group_lf.row(0)}")

    # I assume sample order is all time the same but I'm not sure
    sample2index = {sample: idx for idx, sample in enumerate(group_lf.row(0)[1], start=0)}

    logger.debug(f"{sample2index}")

    transmission_lf = group_lf.with_columns(
        [polars.col(col).arr.get(sample2index[index_name]).alias(f"index_{col}") for col in genotypes_column],
    )

    if mother_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            [polars.col(col).arr.get(sample2index[mother_name]).alias(f"mother_{col}") for col in genotypes_column],
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            [
                polars.lit(3).alias("mother_gt"),
                *[polars.lit(None).alias(f"mother_{col}") for col in genotypes_column if col != "gt"],
            ],
        )

    if father_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            [polars.col(col).arr.get(sample2index[father_name]).alias(f"father_{col}") for col in genotypes_column],
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            [
                polars.lit(3).alias("father_gt"),
                *[polars.lit(None).alias(f"father_{col}") for col in genotypes_column if col != "gt"],
            ],
        )

    transmission_lf = transmission_lf.with_columns(
        (polars.col("index_gt") * 100 + polars.col("mother_gt") * 10 + polars.col("father_gt"))
        .cast(polars.UInt8)
        .alias("origin"),
    )

    transmission_lf = transmission_lf.drop(["sample", *genotypes_column])

    return transmission_lf
