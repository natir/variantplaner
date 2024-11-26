"""Function relate to variant normalization."""

# std import
from __future__ import annotations

import pathlib

# 3rd party import
import polars

try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:  # pragma: no cover
    pass
else:
    cleanup_on_sigterm()

# project import
from variantplaner import normalization

DATA_DIR = pathlib.Path(__file__).parent / "data"


def __generate_chr2len() -> polars.DataFrame:
    """Generate a chr2len dataframe."""
    chr2len = polars.DataFrame(
        data={
            "contig": ["1", "2", "3", "22", "X"],
            "length": [10_000_000, 50_000, 120_000_500, 99_239_816, 10_000],
        },
        schema_overrides={"length": polars.UInt64},
    )

    return chr2len.with_columns(offset=polars.col("length").cum_sum() - polars.col("length"))


def __generate_variants() -> polars.DataFrame:
    """Generate a variant dataframe."""
    return polars.DataFrame(
        data={
            "chr": ["2", "1", "X", "3", "22", "2"],
            "pos": [19910, 3322992, 1292939, 399941, 11111, 16424],
            "ref": ["A", "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG", "A", "C", "T", "A"],
            "alt": ["T", "C", "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG", "G", "C", "*"],
        },
        schema_overrides={"pos": polars.UInt64},
    )


def test_id() -> None:
    """Check id generation."""
    chr2len = __generate_chr2len()
    df = __generate_variants()

    df = normalization.add_variant_id(df.lazy(), chr2len.lazy()).collect()

    assert df.get_column("id").to_list() == [
        344281486070906886,
        114177135718957217,
        13604463283740165373,
        359057238721036295,
        4468882925680590853,
        12704027001632293257,
    ]


def test_partition() -> None:
    """Check part generation."""
    chr2len = __generate_chr2len()
    df = __generate_variants()

    df = normalization.add_variant_id(df.lazy(), chr2len.lazy()).collect()

    df = normalization.add_id_part(df.lazy()).collect()

    assert df.get_column("id_part").to_list() == [9, 3, 255, 9, 124, 255]

    df = normalization.add_id_part(df.lazy(), number_of_bits=8).collect()

    assert df.get_column("id_part").to_list() == [9, 3, 255, 9, 124, 255]

    df = normalization.add_id_part(df.lazy(), number_of_bits=9).collect()

    assert df.get_column("id_part").to_list() == [19, 6, 511, 19, 248, 511]
