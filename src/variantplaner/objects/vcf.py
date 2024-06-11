"""Declare Vcf object."""

# std import
from __future__ import annotations

import enum
import typing

# 3rd party import
import polars

# project import
from variantplaner import normalization
from variantplaner.exception import NoContigsLengthInformationError, NoGenotypeError, NotAVCFError, NotVcfHeaderError
from variantplaner.objects.contigs_length import ContigsLength
from variantplaner.objects.genotypes import Genotypes
from variantplaner.objects.variants import Variants
from variantplaner.objects.vcf_header import VcfHeader

# type checking block
if typing.TYPE_CHECKING:
    import collections
    import pathlib

    from variantplaner import Annotations


class VcfParsingBehavior(enum.Enum):
    """Enumeration use to control behavior of IntoLazyFrame."""

    NOTHING = 0
    """into_lazyframe not have any specific behavior"""

    MANAGE_SV = 1
    """into_lazyframe try to avoid structural variant id collision, SVTYPE/SVLEN info value must be present."""


class Vcf:
    """Object to manage lazyframe as Vcf."""

    def __init__(self):
        """Initialize a Vcf object."""
        self.lf = polars.LazyFrame(schema=Variants.minimal_schema())

        self.header = VcfHeader()

    def from_path(
        self, path: pathlib.Path, chr2len_path: pathlib.Path | None, behavior: VcfParsingBehavior = VcfParsingBehavior.NOTHING
    ) -> None:
        """Populate Vcf object with vcf file."""
        with open(path) as fh:
            try:
                self.header.from_lines(fh)
            except NotVcfHeaderError as e:
                raise NotAVCFError(path) from e

        chr2len = ContigsLength()
        if chr2len.from_vcf_header(self.header) == 0 and (chr2len_path is None or chr2len.from_path(chr2len_path) == 0):
            raise NoContigsLengthInformationError

        self.lf = polars.scan_csv(
            path,
            separator="\t",
            comment_prefix="#",
            has_header=False,
        )

        self.lf = self.lf.rename(dict(zip(self.lf.columns, self.header.column_name(self.lf.width))))
        self.lf = self.lf.cast(Vcf.schema())

        if behavior == VcfParsingBehavior.MANAGE_SV:
            self.lf = self.lf.with_columns(self.header.info_parser({"SVTYPE", "SVLEN"}))

        self.lf = normalization.add_variant_id(self.lf, chr2len.lf)

        if behavior == VcfParsingBehavior.MANAGE_SV:
            self.lf = self.lf.drop("SVTYPE", "SVLEN")

    def variants(self) -> Variants:
        """Get variants of vcf."""
        return self.lf.select(Variants.minimal_schema())

    def set_variants(self, variants: Variants) -> None:
        """Set variants of vcf."""
        self.lf = variants.lf

    def genotypes(self, format_str: str = "GT:AD:DP:GQ") -> Genotypes:
        """Get genotype of vcf."""
        if "format" not in self.lf.columns:
            raise NoGenotypeError

        lf = self.lf.select([*self.lf.columns[self.lf.columns.index("format") :]])

        # Clean bad variant
        lf = lf.filter(polars.col("format").str.starts_with(format_str)).select(*lf.columns[1:])

        # Found index of genotype value
        col_index = {
            key: index
            for (index, key) in enumerate(
                format_str.split(":"),
            )
        }

        # Pivot value
        genotypes = Genotypes()
        genotypes.lf = lf.melt(id_vars=["id"]).with_columns(
            [
                polars.col("id"),
                polars.col("variable").alias("sample"),
                polars.col("value").str.split(":"),
            ],
        )

        # Split genotype column in sub value
        col2expr = self.header.format_parser()

        genotypes.lf = genotypes.lf.with_columns(
            [polars.col("value").list.get(index).pipe(function=col2expr[col], col_name=col) for col, index in col_index.items()],
        )

        # Select intrusting column
        genotypes.lf = genotypes.lf.select(["id", "sample", *[col.lower() for col in col_index]])

        if "gt" in genotypes.lf.columns:
            genotypes.lf = genotypes.lf.filter(polars.col("gt") != 0)

        return genotypes

    def add_genotypes(self, genotypes_lf: Genotypes) -> None:
        """Add genotypes information in vcf."""
        for sample in genotypes_lf.samples_names():
            geno2sample = (
                genotypes_lf.lf.filter(polars.col("sample") == sample)
                .rename(
                    {col: f"{sample}_{col}" for col in genotypes_lf.lf.columns[2:]},
                )
                .drop("sample")
            )

            self.lf = self.lf.join(geno2sample, on="id", how="outer_coalesce")

    def annotations(self, select_info: set[str] | None = None) -> Annotations:
        """Get annotations of vcf."""
        lf = self.lf.with_columns(self.lf.header.info_parser(select_info))

        return lf.drop("chr", "pos", "ref", "alt", "format", "info")

    @classmethod
    def schema(cls) -> collections.abc.Mapping[polars.type_aliases.ColumnNameOrSelector, polars.type_aliases.PolarsDataType]:
        """Get schema of Vcf polars.LazyFrame."""
        return {
            "chr": polars.String,
            "pos": polars.UInt64,
            "vid": polars.String,
            "ref": polars.String,
            "alt": polars.String,
            "qual": polars.String,
            "filter": polars.String,
            "info": polars.String,
        }
