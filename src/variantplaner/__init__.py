"""VariantPlaner, a tool kit to manage many variants without many cpu and ram resource.

Convert a vcf in parquet, convert annotations in parquet, convert parquet in vcf.

But also build a file struct to get a fast variant database interrogations time.
"""

from __future__ import annotations

import base64

from variantplaner import extract, generate, normalization, struct
from variantplaner.objects import (
    Annotations,
    ContigsLength,
    Genotypes,
    Pedigree,
    Variants,
    Vcf,
    VcfHeader,
    VcfParsingBehavior,
)


def int2string(value: int) -> str:
    return base64.urlsafe_b64encode(
        value.to_bytes(
            (value.bit_length() + 7) // 8,
            byteorder="little",
            signed=True,
        )
    ).decode("utf-8")


__all__: list[str] = [
    "Annotations",
    "ContigsLength",
    "Genotypes",
    "Pedigree",
    "Variants",
    "Vcf",
    "VcfHeader",
    "VcfParsingBehavior",
    "exception",
    "extract",
    "generate",
    "io",
    "normalization",
    "struct",
]
__version__: str = "0.3.1"
