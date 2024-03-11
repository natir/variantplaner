"""Declare Vcf object."""

# std import
from __future__ import annotations

import re
import typing

# 3rd party import
import polars

# project import
from variantplaner.objects.csv import Csv

if typing.TYPE_CHECKING:
    import pathlib
    import sys

    from variantplaner.objects.csv import ScanCsv
    from variantplaner.objects.vcf_header import VcfHeader

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack


class ContigsLength:
    """Store contigs -> length information."""

    def __init__(self):
        """Initialise a contigs length."""
        self.lf = polars.LazyFrame(
            schema={
                "contig": polars.String,
                "length": polars.UInt64,
                "offset": polars.UInt64,
            }
        )

    def from_vcf_header(self, header: VcfHeader) -> int:
        """Fill a object with VcfHeader.

        Argument:
           header: VcfHeader

        Returns: Number of contigs line view
        """
        contigs_id = re.compile(r"ID=(?P<id>[^,]+)")
        contigs_len = re.compile(r"length=(?P<length>[^,>]+)")

        count = 0
        contigs2len: dict[str, list] = {"contig": [], "length": []}
        for contig_line in header.contigs:
            if (len_match := contigs_len.search(contig_line)) and (id_match := contigs_id.search(contig_line)):
                contigs2len["contig"].append(id_match.groupdict()["id"])
                contigs2len["length"].append(int(len_match.groupdict()["length"]))
            count += 1

        self.lf = polars.LazyFrame(contigs2len, schema={"contig": polars.String, "length": polars.UInt64})

        self.__compute_offset()

        return count

    def from_path(self, path: pathlib.Path, /, **scan_csv_args: Unpack[ScanCsv]) -> int:
        """Fill object with file point by pathlib.Path.

        Argument:
        path: path of input file

        Returns: Number of contigs line view
        """
        csv = Csv()
        csv.from_path(path, **scan_csv_args)
        self.lf = csv.lf

        self.__compute_offset()

        return self.lf.collect().shape[0]

    def __compute_offset(self) -> None:
        self.lf = self.lf.with_columns(offset=polars.col("length").cum_sum() - polars.col("length"))
        self.lf = self.lf.cast({"offset": polars.UInt64})
