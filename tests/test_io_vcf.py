"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import filecmp
import typing

if typing.TYPE_CHECKING:  # pragma: no cover
    import sys

    if sys.version_info >= (3, 11):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    T = typing.TypeVar("T")
    P = ParamSpec("P")

# 3rd party import
import pathlib

import polars
import polars.testing
import pytest

# project import
from variantplaner import exception, extract, io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_extract_header() -> None:
    """Check extract_header."""
    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")

    assert header == [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FILTER=<ID=base_qual,Description="alt median base quality">',
        '##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">',
        '##FILTER=<ID=contamination,Description="contamination">',
        '##FILTER=<ID=duplicate,Description="evidence for alt allele is overrepresented by apparent duplicates">',
        '##FILTER=<ID=fragment,Description="abs(ref - alt) median fragment length">',
        '##FILTER=<ID=germline,Description="Evidence indicates this site is germline, not somatic">',
        '##FILTER=<ID=haplotype,Description="Variant near filtered variant on same haplotype.">',
        '##FILTER=<ID=low_allele_frac,Description="Allele fraction is below specified threshold">',
        '##FILTER=<ID=map_qual,Description="ref - alt median mapping quality">',
        '##FILTER=<ID=multiallelic,Description="Site filtered because too many alt alleles pass tumor LOD">',
        '##FILTER=<ID=n_ratio,Description="Ratio of N to alt exceeds specified ratio">',
        '##FILTER=<ID=normal_artifact,Description="artifact_in_normal">',
        '##FILTER=<ID=numt_chimera,Description="NuMT variant with too many ALT reads originally from autosome">',
        '##FILTER=<ID=numt_novel,Description="Alt depth is below expected coverage of NuMT in autosome">',
        '##FILTER=<ID=orientation,Description="orientation bias detected by the orientation bias mixture model">',
        '##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">',
        '##FILTER=<ID=position,Description="median distance of alt variants from end of reads">',
        '##FILTER=<ID=slippage,Description="Site filtered due to contraction of short tandem repeat region">',
        '##FILTER=<ID=strand_bias,Description="Evidence for alt allele comes from one read direction only">',
        '##FILTER=<ID=strict_strand,Description="Evidence for alt allele is not represented in both directions">',
        '##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1\tsample_2\tsample_3",
    ]


def test_extract_header_exception() -> None:
    """Check extract_header."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.extract_header(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.extract_header(DATA_DIR / "only_header.vcf")


def test_info2expr_no_info_vcf() -> None:
    """Check info2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    expressions = io.vcf.info2expr(header, DATA_DIR / "no_info.vcf")

    assert expressions == []


def test_info2expr() -> None:
    """Check info2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")
    expressions = io.vcf.info2expr(header, DATA_DIR / "no_genotypes.vcf")

    assert len(expressions) == 21


def test_info2expr_exception() -> None:
    """Check info2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.info2expr(header[:-1], DATA_DIR / "no_genotypes.vcf")


def test_format2expr_no_format_vcf() -> None:
    """Check format2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")
    assert (
        io.vcf.format2expr(
            header,
            DATA_DIR / "no_genotypes.vcf",
        )
        == {}
    )


def test_format2expr() -> None:
    """Check format2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")

    assert set(io.vcf.format2expr(header, DATA_DIR / "no_info.vcf").keys()) == {"AD", "DP", "GQ", "GT"}


def test_format2expr_exception() -> None:
    """Check format2expr."""
    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.format2expr(header[:-1], DATA_DIR / "no_info.vcf")


def test_sample_index() -> None:
    """Check sample index."""
    truth = {"sample_1": 0, "sample_2": 1, "sample_3": 2}

    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    header.reverse()
    sample2idx = io.vcf.sample_index(header, DATA_DIR / "no_info.vcf")

    assert sample2idx
    assert len(truth) == len(sample2idx)
    assert all(v == sample2idx[k] for k, v in truth.items())


def test_sample_index_no_genotypes() -> None:
    """Check sample index."""
    header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")
    assert io.vcf.sample_index(header, DATA_DIR / "no_genotypes.vcf") is None


def test_sample_index_exception() -> None:
    """Check sample index exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.extract_header(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.sample_index([], DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.extract_header(DATA_DIR / "only_header.vcf")


def test__column_name() -> None:
    """Check __column_name."""
    header = io.vcf.extract_header(DATA_DIR / "no_info.vcf")
    assert io.vcf.__column_name(header, DATA_DIR / "no_info.vcf") == [
        "chr",
        "pos",
        "vid",
        "ref",
        "alt",
        "qual",
        "filter",
        "info",
        "format",
        "sample_1",
        "sample_2",
        "sample_3",
    ]

    with pytest.raises(exception.NotAVCFError):
        io.vcf.__column_name(header[:-1], DATA_DIR / "no_info.vcf")


def test_into_lazyframe() -> None:
    """Check into lazyframe."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf", DATA_DIR / "grch38.92.csv")

    polars.testing.assert_frame_equal(truth, lf)

    truth = polars.scan_parquet(DATA_DIR / "no_genotypes.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf", DATA_DIR / "grch38.92.csv")

    polars.testing.assert_frame_equal(truth, lf)


def test_into_lazyframe_sv() -> None:
    """Check into lazyframe work on structural variant."""
    input_path = DATA_DIR / "sv.vcf"

    lf = io.vcf.into_lazyframe(
        input_path,
        DATA_DIR / "grch38.92.csv",
        extension=io.vcf.IntoLazyFrameExtension.MANAGE_SV,
    )

    polars.testing.assert_frame_equal(lf, polars.scan_parquet(DATA_DIR / "sv.parquet"), check_row_order=False)


def test_into_lazyframe_exception() -> None:
    """Check into_lazyframe exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "no_info.tsv", DATA_DIR / "grch38.92.csv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "only_header.vcf", DATA_DIR / "grch38.92.csv")


def test_build_rename_column() -> None:
    """Check build_rename_column."""
    assert io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt") == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": ".",
        "FILTER": ".",
        "INFO": [],
        "FORMAT": "",
        "sample": {},
    }

    assert io.vcf.build_rename_column(
        "chr",
        "pos",
        "id",
        "ref",
        "alt",
        "quality",
        "FILTER",
        [("GENE", "gene_name")],
        "GT:AD:DP:GQ",
    ) == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": "quality",
        "FILTER": "FILTER",
        "INFO": [("GENE", "gene_name")],
        "FORMAT": "GT:AD:DP:GQ",
        "sample": {},
    }

    assert io.vcf.build_rename_column(
        "chr",
        "pos",
        "id",
        "ref",
        "alt",
        "quality",
        "FILTER",
        [("GENE", "gene_name")],
        "GT:AD:DP:GQ",
        {
            "sample": {
                "gt": "sample_gt",
                "ad": "sample_ad",
                "dp": "sample_dp",
                "gq": "sample_gq",
            },
        },
    ) == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": "quality",
        "FILTER": "FILTER",
        "INFO": [("GENE", "gene_name")],
        "FORMAT": "GT:AD:DP:GQ",
        "sample": {"sample": {"gt": "sample_gt", "ad": "sample_ad", "dp": "sample_dp", "gq": "sample_gq"}},
    }


def test_from_lazyframe(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe."""
    tmp_file = tmp_path / "output.vcf"

    lf = polars.scan_parquet(DATA_DIR / "no_info.parquet")

    io.vcf.from_lazyframe(lf, tmp_file)

    assert filecmp.cmp(tmp_file, DATA_DIR / "no_info.parquet2vcf.vcf")

    io.vcf.from_lazyframe(lf, tmp_file, io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt", "vid", "vid"))

    assert filecmp.cmp(tmp_file, DATA_DIR / "no_info.parquet2vcf.vcf")


def test_from_lazyframe_qual_filter_info(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe qual, filter, info."""
    output_path = tmp_path / "output.vcf"
    input_path = DATA_DIR / "no_genotypes.vcf"

    vcf_header = io.vcf.extract_header(input_path)
    info_parser = io.vcf.info2expr(vcf_header, input_path)
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf", DATA_DIR / "grch38.92.csv").with_columns(info_parser)

    io.vcf.from_lazyframe(
        lf,
        output_path,
        io.vcf.build_rename_column(
            "chr",
            "pos",
            "id",
            "ref",
            "alt",
            "qual",
            "filter",
            [(col_name, col_name) for col_name in lf.columns if col_name.isupper()],
        ),
    )

    assert filecmp.cmp(output_path, DATA_DIR / "no_genotypes.vcf2parquet2vcf.vcf")


def test_from_lazyframe_qual_filter_format(tmp_path: pathlib.Path) -> None:
    """Check from_lazyframe qual, filter, info."""
    output_path = tmp_path / "output.vcf"
    input_path = DATA_DIR / "no_info.vcf"

    format_string = "GT:AD:DP:GQ"
    vcf = io.vcf.into_lazyframe(input_path, DATA_DIR / "grch38.92.csv")
    vcf_header = io.vcf.extract_header(input_path)
    sample_name = io.vcf.sample_index(vcf_header, input_path)
    if sample_name is None:
        raise AssertionError  # pragma: no cover Not reachable code

    format2expr = io.vcf.format2expr(vcf_header, input_path)

    variants = extract.variants(vcf)
    genotypes = extract.genotypes(vcf, format2expr, format_string)

    merge = extract.merge_variants_genotypes(variants, genotypes, list(sample_name.keys()))
    sample2vcf_col2polars_col: dict[str, dict[str, str]] = {}
    for sample in sample_name:
        sample2vcf_col2polars_col[sample] = {}
        for format_col in format_string.split(":"):
            sample2vcf_col2polars_col[sample][format_col] = f"{sample}_{format_col.lower()}"

    io.vcf.from_lazyframe(
        merge,
        output_path,
        io.vcf.build_rename_column(
            "chr",
            "pos",
            "id",
            "ref",
            "alt",
            None,
            None,
            [],
            "GT:AD:DP:GQ",
            sample2vcf_col2polars_col,
        ),
    )

    assert filecmp.cmp(output_path, DATA_DIR / "no_info.vcf2parquet2vcf.vcf")


def test_add_info_column() -> None:
    """Check add_info_column."""
    vcf_header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")
    info_parser = io.vcf.info2expr(vcf_header, DATA_DIR / "no_genotypes.vcf")
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf", DATA_DIR / "grch38.92.csv").with_columns(info_parser)

    info = io.vcf.add_info_column(lf, [(col_name, col_name) for col_name in lf.columns if col_name.isupper()])

    truth = [
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2193183;CLNDN=Inborn_genetic_diseases;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.69134A>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2238986;CLNDN=Inborn_genetic_diseases;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.69581C>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2386655;CLNDN=Inborn_genetic_diseases;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.69682G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2278803;CLNDN=Inborn_genetic_diseases;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.69769T>C;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2333177;CLNDN=Inborn_genetic_diseases;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.69995G>C;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=OR4F5:79501;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1983057;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925946C>G;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1003021;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925952G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=1640863258",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1632777;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925956C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001819|synonymous_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2129477;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925961A>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1600580;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925969C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1396033;CLNDN=Inborn_genetic_diseases|not_provided;CLNDNINCL=.;CLNDISDB=MeSH:D030342,MedGen:C0950123|MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925976T>C;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1986319;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925980C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001819|synonymous_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1570515;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.925986C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001819|synonymous_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1502313;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926003C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1545352;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926010G>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001819|synonymous_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1473095;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926014G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001575|splice_donor_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2034738;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926018G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1550067;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926025G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=2155344;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926026G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1641615;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926027C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1629212;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.926029C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1631640;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.930136T>C;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
        "AF_ESP=.;AF_EXAC=.;AF_TGP=.;ALLELEID=1563949;CLNDN=not_provided;CLNDNINCL=.;CLNDISDB=MedGen:CN517202;CLNDISDBINCL=.;CLNHGVS=NC_000001.11:g.930139CCT[1];CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNSIGCONF=.;CLNSIGINCL=.;CLNVC=Microsatellite;CLNVCSO=SO:0000289;CLNVI=.;DBVARID=.;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;RS=.",
    ]

    assert info.select([polars.col("INFO")]).collect()["INFO"].to_list() == truth


def test_generate_header() -> None:
    """Check generate header."""
    vcf_header = io.vcf.extract_header(DATA_DIR / "no_genotypes.vcf")
    info_parser = io.vcf.info2expr(vcf_header, DATA_DIR / "no_genotypes.vcf")
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf", DATA_DIR / "grch38.92.csv").with_columns(info_parser)

    assert (
        io.vcf.__generate_header(lf)
        == """##fileformat=VCFv4.3
##source=VariantPlanner
"""
    )

    info = [(col_name, col_name) for col_name in lf.columns if col_name.isupper()]

    header = io.vcf.__generate_header(lf, info)
    truth = """##fileformat=VCFv4.3
##source=VariantPlanner
##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="Unknow">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="Unknow">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="Unknow">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Unknow">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Unknow">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Unknow">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Unknow">
##INFO=<ID=CLNVI,Number=.,Type=String,Description="Unknow">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="Unknow">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Unknow">
##INFO=<ID=MC,Number=.,Type=String,Description="Unknow">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Unknow">
##INFO=<ID=RS,Number=.,Type=String,Description="Unknow">
"""

    assert header == truth
