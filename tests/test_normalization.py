"""Tests for normalization."""

# std import
from __future__ import annotations

import io
import logging
import pathlib
import typing

# 3rd party import
import polars

# project import
from variantplaner import normalization

if typing.TYPE_CHECKING:  # pragma: no cover
    import pytest

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_normalization_warning(caplog: pytest.LogCaptureFixture) -> None:
    """Test normalization produce warning on not classic alt."""
    caplog.set_level(logging.WARNING)

    df = polars.read_csv(
        io.StringIO(
            """
chr,pos,vid,ref,alt,qual,filter,info
1,10146,.,AC,<DEL>,160.6,.,.
""",
        ),
    )

    df = normalization.add_variant_id(df.lazy()).collect()

    assert caplog.record_tuples == [
        ("normalization", 30, "Alternative column contains not classic sequence this can create variant collisions"),
    ]


def test_normalization_no_warning(caplog: pytest.LogCaptureFixture) -> None:
    """Test normalization produce warning on not classic alt."""
    caplog.set_level(logging.WARNING)

    df = polars.read_csv(
        io.StringIO(
            """
chr,pos,vid,ref,alt,qual,filter,info
1,10146,.,AC,A,160.6,.,.
""",
        ),
    )

    df = normalization.add_variant_id(df.lazy()).collect()

    assert caplog.record_tuples == []
