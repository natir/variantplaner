"""Tests for the `io.vcf` module."""

# std import
from __future__ import annotations

import filecmp

# 3rd party import
import pathlib

import polars
import pytest

# project import
from variantplanner import exception, io

DATA_DIR = pathlib.Path(__file__).parent / "data"


def test_info2expr_no_info_vcf() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_info.vcf")

    assert expressions == []


def test_info2expr() -> None:
    """Check info2expr."""
    expressions = io.vcf.info2expr(DATA_DIR / "no_genotypes.vcf")

    assert len(expressions) == 21


def test_sample_index() -> None:
    """Check sample index."""
    truth = {"sample_1": 0, "sample_2": 1, "sample_3": 2}

    sample2idx = io.vcf.sample_index(DATA_DIR / "no_info.vcf")

    assert sample2idx
    assert len(truth) == len(sample2idx)
    assert all(v == sample2idx[k] for k, v in truth.items())


def test_sample_index_no_genotypes() -> None:
    """Check sample index."""
    assert io.vcf.sample_index(DATA_DIR / "no_genotypes.vcf") is None


def test_sample_index_exception() -> None:
    """Check sample index exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.sample_index(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.sample_index(DATA_DIR / "only_header.vcf")


def test_into_lazyframe() -> None:
    """Check into lazyframe."""
    truth = polars.scan_parquet(DATA_DIR / "no_info.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_info.vcf")

    polars.testing.assert_frame_equal(truth, lf)

    truth = polars.scan_parquet(DATA_DIR / "no_genotypes.parquet")

    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf")

    polars.testing.assert_frame_equal(truth, lf)


def test_into_lazyframe_exception() -> None:
    """Check into_lazyframe exception."""
    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "no_info.tsv")

    with pytest.raises(exception.NotAVCFError):
        io.vcf.into_lazyframe(DATA_DIR / "only_header.vcf")


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
        "INFO": {},
    }

    assert io.vcf.build_rename_column("chr", "pos", "id", "ref", "alt", "quality", "FILTER", {"GENE": "gene_name"}) == {
        "#CHROM": "chr",
        "POS": "pos",
        "ID": "id",
        "REF": "ref",
        "ALT": "alt",
        "QUAL": "quality",
        "FILTER": "FILTER",
        "INFO": {"GENE": "gene_name"},
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

    info_parser = io.vcf.info2expr(input_path)
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf").with_columns(info_parser)

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
            {col_name: col_name for col_name in lf.columns if col_name.isupper()},
        ),
    )

    assert filecmp.cmp(output_path, DATA_DIR / "no_genotypes.vcf2parquet2vcf.vcf")


def test_add_info_column() -> None:
    """Check add_info_column."""
    info_parser = io.vcf.info2expr(DATA_DIR / "no_genotypes.vcf")
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf").with_columns(info_parser)

    info = io.vcf.add_info_column(lf, {col_name: col_name for col_name in lf.columns if col_name.isupper()})

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
    info_parser = io.vcf.info2expr(DATA_DIR / "no_genotypes.vcf")
    lf = io.vcf.into_lazyframe(DATA_DIR / "no_genotypes.vcf").with_columns(info_parser)

    assert (
        io.vcf.generate_header(lf)
        == """##fileformat=VCFv4.3
##source=VariantPlanner
"""
    )

    info = {col_name: col_name for col_name in lf.columns if col_name.isupper()}

    header = io.vcf.generate_header(lf, info)
    truth = """##fileformat=VCFv4.3
##source=VariantPlanner
##INFO=<ID=AF_ESP,Number=1,Type=String,Description="Unknow">
##INFO=<ID=AF_EXAC,Number=1,Type=String,Description="Unknow">
##INFO=<ID=AF_TGP,Number=1,Type=String,Description="Unknow">
##INFO=<ID=ALLELEID,Number=1,Type=String,Description="Unknow">
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
