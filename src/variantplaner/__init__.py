"""VariantPlaner.

A tool kit to manage many variant without many cpu and ram ressource.

Convert vcf in parquet, convert annotations in parquet, convert parquet in vcf.
Build a file struct to get a fast variants database interogations time.
"""

from __future__ import annotations

from variantplaner import exception, extract, generate, io, normalization, struct

__all__: list[str] = ["exception", "extract", "generate", "io", "normalization", "struct"]
