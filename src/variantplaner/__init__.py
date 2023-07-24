"""VariantPlaner, a tool kit to manage many variants without many cpu and ram resource.

Convert a vcf in parquet, convert annotations in parquet, convert parquet in vcf.

But also build a file struct to get a fast variant database interrogations time.
"""

from __future__ import annotations

from variantplaner import exception, extract, generate, io, normalization, struct

__all__: list[str] = ["exception", "extract", "generate", "io", "normalization", "struct"]
