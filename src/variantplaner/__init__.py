"""VariantPlaner, a tool kit to manage many variants without many cpu and ram resource.

Convert a vcf in parquet, convert annotations in parquet, convert parquet in vcf.

But also build a file struct to get a fast variant database interrogations time.
"""

from __future__ import annotations

from variantplaner import extract, generate, normalization, struct
from variantplaner.objects import Annotations, ContigsLength, Genotypes, Variants, Vcf, VcfHeader, VcfParsingBehavior

__all__: list[str] = [
    "Annotations",
    "ContigsLength",
    "Genotypes",
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
__version__: str = "0.2.4"
