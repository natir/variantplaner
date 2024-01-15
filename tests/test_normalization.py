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
            "chr": ["1", "2", "3", "22", "X"],
            "length": [10_000_000, 50_000, 120_000_500, 99_239_816, 10_000],
        },
        schema_overrides={"length": polars.UInt64},
    )

    return chr2len.with_columns(offset=polars.col("length").cum_sum() - polars.col("length"))


def __generate_variants() -> polars.DataFrame:
    """Generate a variant dataframe."""
    return polars.DataFrame(
        data={
            "chr": ["2", "1", "X", "3", "22"],
            "pos": [19910, 3322992, 1292939, 399941, 11111],
            "ref": ["A", "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG", "A", "C", "T"],
            "alt": ["T", "C", "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG", "G", "C"],
        },
        schema_overrides={"pos": polars.UInt64},
    )


def test_id() -> None:
    """Check id generation."""
    chr2len = __generate_chr2len()
    df = __generate_variants()

    df = normalization.add_variant_id(df.lazy(), chr2len.lazy()).collect()

    assert df.get_column("id").to_list() == [
        344281487144648706,
        14605542477335234092,
        17793944462439896137,
        359057239794778119,
        4468882926754332681,
    ]


def test_partition() -> None:
    """Check part generation."""
    chr2len = __generate_chr2len()
    df = __generate_variants()

    df = normalization.add_variant_id(df.lazy(), chr2len.lazy()).collect()

    df = normalization.add_id_part(df.lazy()).collect()

    assert df.get_column("id_part").to_list() == [9, 255, 255, 9, 124]
