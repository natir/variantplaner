"""Function to generate information."""

# std import
from __future__ import annotations

import logging

# 3rd party import
import polars

# project import
from variantplaner.exception import NoGTError

logger = logging.getLogger("generate")

gt2chr = {i: chr(i + 33) for i in range(94)}


def transmission_ped(
    genotypes_lf: polars.LazyFrame,
    pedigree_lf: polars.LazyFrame,
) -> polars.DataFrame:
    """Compute transmission of each variants.

    **Warning**: only the first sample with two parent are considered.

    Args:
        genotypes_lf: Genotypes [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html), `gt` column are required.
        pedigree_lf: Pedigree [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html).

    Returns:
         DataFrame with transmission information

    Raises:
        NoGTError: If genotypes_lf not contains gt column.
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
        genotypes_lf: Genotypes [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html), `gt` column are required.
        index_name: Sample name of index case.
        mother_name: Sample name of mother.
        father_name: Sample name of father.

    Returns:
         [polars.DataFrame](https://pola-rs.github.io/polars/py-polars/html/reference/dataframe/index.html) with transmission information.
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

    samples = sorted(genotypes_lf.select("sample").unique().collect().get_column("sample").to_list())

    logger.debug(f"{samples=}")

    group_lf = genotypes_lf.group_by("id").all().collect()

    group_lf = group_lf.filter(polars.col("sample").list.contains(index_name))

    # I assume sample order is all time the same but I'm not sure
    sample2index = {sample: idx for idx, sample in enumerate(samples, start=0)}

    logger.debug(f"{index_name=} {mother_name=} {father_name=}")
    logger.debug(f"{sample2index}")

    transmission_lf = group_lf.with_columns(
        [polars.col(col).list.get(sample2index[index_name]).alias(f"index_{col}") for col in genotypes_column],
    )

    if mother_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            [polars.col(col).list.get(sample2index[mother_name]).alias(f"mother_{col}") for col in genotypes_column],
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            [polars.lit(None).alias(f"mother_{col}") for col in genotypes_column],
        )

    if father_name in sample2index:
        transmission_lf = transmission_lf.with_columns(
            [polars.col(col).list.get(sample2index[father_name]).alias(f"father_{col}") for col in genotypes_column],
        )
    else:
        transmission_lf = transmission_lf.with_columns(
            [polars.lit(None).alias(f"father_{col}") for col in genotypes_column],
        )

    polars.Config.set_tbl_width_chars(1000)
    polars.Config.set_tbl_cols(65)

    transmission_lf = transmission_lf.with_columns(
        polars.concat_str(
            polars.col("index_gt").map_dict(gt2chr, default="~", return_dtype=polars.Utf8),
            polars.col("mother_gt").fill_null(94).map_dict(gt2chr, default="~", return_dtype=polars.Utf8),
            polars.col("father_gt").fill_null(94).map_dict(gt2chr, default="~", return_dtype=polars.Utf8),
        ).alias("origin"),
    )

    return transmission_lf.drop(["sample", *genotypes_column])
