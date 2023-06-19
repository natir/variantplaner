"""Function use to normalize data."""

# std import
from __future__ import annotations

import logging

# 3rd party import
import polars

# project import


logger = logging.getLogger("normalization")


def chromosome2integer(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Convert chromosome string in number.

    Chromosome name mapping table:
      - X:  23
      - Y:  24
      - MT: 25

    If chromosome value can't be convert in number row is remove.

    Args:
        lf: [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) contains chr column

    Returns:
        [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) with chr column normalized

    """
    lf = lf.with_columns(
        [
            polars.col("chr")
            .str.replace("chr", "")
            .str.replace("X", "23")
            .str.replace("Y", "24")
            .str.replace("MT", "25")
            .str.parse_int(10, strict=False)
            .cast(polars.UInt8),
            polars.col("pos").cast(polars.UInt64),
        ],
    )

    return lf.filter(~polars.col("chr").is_null())


def add_variant_id(lf: polars.LazyFrame) -> polars.LazyFrame:
    """Add a column id of variants.


    This id is compute by 64 bit hash on chromosome, position, reference sequence and alternative sequence.


    If lf.columns contains SVTYPE and SVLEN variant with regex group in alt <([^:]+).*> match SVTYPE are replace by concatenation of SVTYPE and SVLEN first value.


    Colision risk:
        - human genome size: 3,117,275,501 bp
        - number of variant if each base have all sustitution $3,117,275,501 * 4 = 12,469,102,004$
        - probablity of colision $12,469,102,004 / 2^{64} = 6.7595137 * 10^{-10}$

    Args:
        lf: [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) contains chr, pos, ref, alt columns

    Returns:
        [polars.LazyFrame](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) with chr column normalized
    """
    if "SVTYPE" in lf.columns and "SVLEN" in lf.columns:
        lf = lf.with_columns(
            polars.when(
                polars.col("alt").str.replace("<(?<type>[^:]+).*>", "$type") == polars.col("SVTYPE"),
            )
            .then(
                polars.col("alt").str.replace(
                    ".+",
                    polars.concat_str(
                        [polars.col("SVTYPE"), polars.col("SVLEN").list.get(0)],
                        separator="-",
                    ),
                ),
            )
            .otherwise(
                polars.col("alt"),
            ),
        )

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
        )
        .hash()
        .alias("id"),
    )
