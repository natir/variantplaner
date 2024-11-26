"""Module to store variantplaner object."""

# std import
from __future__ import annotations

# 3rd party import
# project import
from variantplaner.objects.annotations import Annotations
from variantplaner.objects.contigs_length import ContigsLength
from variantplaner.objects.genotypes import Genotypes
from variantplaner.objects.pedigree import Pedigree
from variantplaner.objects.variants import Variants
from variantplaner.objects.vcf import Vcf, VcfParsingBehavior
from variantplaner.objects.vcf_header import VcfHeader

__all__: list[str] = [
    "Annotations",
    "ContigsLength",
    "Genotypes",
    "Pedigree",
    "Variants",
    "Vcf",
    "VcfHeader",
    "VcfParsingBehavior",
]
