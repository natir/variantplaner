"""Tests for cli."""

# std import
from __future__ import annotations

import filecmp
import pathlib

# 3rd party import
import polars
import polars.testing
from click.testing import CliRunner

try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:  # pragma: no cover
    pass
else:
    cleanup_on_sigterm()


# project import
from variantplaner import cli

DATA_DIR = pathlib.Path(__file__).parent / "data"


MERGE_IDS = {
    21788637528068,
    22420668831105,
    22531532652550,
    28503684677645,
    216629694693385,
    2865209939989626902,
    2865214825649143842,
    2865214829944111138,
    2865214831957377026,
    2865214831957377025,
    3382250923559485444,
    3382273467842822145,
    3382274273283407906,
    3382286558903140374,
    3382286559037358105,
    6356557166639316998,
    6356557209588989953,
    6356558437949636610,
    6356559434382049293,
    6356561358527397901,
    6174051668938719235,
    6174052457065218057,
    6174053404105506822,
    6174053880846876681,
    6174054230886711305,
    148464268738563,
    149424193929223,
    149641089777676,
    149827920855049,
    150313252159501,
    1988454028148743,
    1988466913050636,
    1988475502985222,
    1988486240403458,
    1988503420272646,
    1988518452658185,
    1988527042592774,
    1988539927494662,
    1988576434716678,
    1988591467102222,
    1988600057036812,
    1988608646971404,
    1988623679356940,
    1988625826840588,
    1988627974324230,
    1988632269291526,
    1997451984633865,
    1997452387287654,
}


def test_vcf2parquet(tmp_path: pathlib.Path) -> None:
    """vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "no_info.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
            "-g",
            str(genotypes_path),
        ],
    )

    assert result.exit_code == 0
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "no_info.variants.parquet"),
        polars.scan_parquet(variants_path),
        check_row_order=False,
    )
    try:
        polars.testing.assert_frame_equal(
            polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet").sort("id"),
            polars.scan_parquet(genotypes_path).sort("id"),
        )
    except OverflowError:  # pragma: no cover
        # TODO remove this
        truth = polars.scan_parquet(DATA_DIR / "no_info.genotypes.parquet")
        lf = polars.scan_parquet(genotypes_path)
        assert truth.columns == lf.columns
        assert truth.dtypes == lf.dtypes
        assert truth.width == lf.width


def test_vcf2parquet_ask_annotations(tmp_path: pathlib.Path) -> None:
    """Ask annotations vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
            "-a",
            str(annotations_path),
        ],
    )

    assert result.exit_code == 0
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "no_genotypes.variants.parquet"),
        polars.scan_parquet(variants_path),
        check_row_order=False,
    )

    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "no_genotypes.annotations.parquet").sort("id"),
        polars.scan_parquet(annotations_path).sort("id"),
        check_column_order=False,
    )


def test_vcf2parquet_not_ask_genotypes(tmp_path: pathlib.Path) -> None:
    """Not ask genotype vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "no_info.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
        ],
    )

    assert result.exit_code == 0
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "no_info.variants.parquet"),
        polars.scan_parquet(variants_path),
        check_row_order=False,
    )


def test_vcf2parquet_not_vcf(tmp_path: pathlib.Path) -> None:
    """Not a vcf input vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "no_info.tsv"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
            "-g",
            str(genotypes_path),
        ],
    )

    assert result.exit_code == 11


def test_vcf2parquet_no_genotype(tmp_path: pathlib.Path) -> None:
    """No genotype in vcf vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
            "-g",
            str(genotypes_path),
        ],
    )

    assert result.exit_code == 12


def test_vcf2parquet_sv(tmp_path: pathlib.Path) -> None:
    """vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "sv.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
        ],
    )

    assert result.exit_code == 0
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "sv.variants.parquet"),
        polars.scan_parquet(variants_path),
        check_row_order=False,
    )


def test_vcf2parquet_sv_genotypes(tmp_path: pathlib.Path) -> None:
    """vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "vcf2parquet",
            "-i",
            str(DATA_DIR / "sv.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-v",
            str(variants_path),
            "-g",
            str(genotypes_path),
            "-f",
            "GT:GQ:CN:CNQ",
        ],
    )

    assert result.exit_code == 0
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "sv.variants.parquet"),
        polars.scan_parquet(variants_path),
        check_row_order=False,
    )
    polars.testing.assert_frame_equal(
        polars.scan_parquet(DATA_DIR / "sv.genotypes.parquet"),
        polars.scan_parquet(genotypes_path),
        check_row_order=False,
    )


def test_parquet2vcf(tmp_path: pathlib.Path) -> None:
    """parquet2vcf run."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    result = runner.invoke(cli.main, ["parquet2vcf", "-i", str(DATA_DIR / "no_info.parquet"), "-o", str(variants_path)])

    assert result.exit_code == 0

    filecmp.cmp(variants_path, DATA_DIR / "no_info.parquet2vcf.vcf")


def test_parquet2vcf_add_genotype(tmp_path: pathlib.Path) -> None:
    """parquet2vcf run."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "parquet2vcf",
            "-i",
            str(DATA_DIR / "no_info.variants.parquet"),
            "-o",
            str(variants_path),
            "-g",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "-F",
            "GT:AD:DP:GQ",
        ],
    )

    assert result.exit_code == 0

    filecmp.cmp(variants_path, DATA_DIR / "no_info.parquet2vcf_genotypes.vcf")


def test_struct_variants(tmp_path: pathlib.Path) -> None:
    """Basic struct variant run."""
    merge_path = tmp_path / "merge.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "struct",
            "-i",
            str(DATA_DIR / "no_genotypes.variants.parquet"),
            "-i",
            str(DATA_DIR / "no_info.variants.parquet"),
            "variants",
            "-o",
            str(merge_path),
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(merge_path)

    assert set(lf.collect().get_column("id").to_list()) == MERGE_IDS


def test_struct_genotypes(tmp_path: pathlib.Path) -> None:
    """Basic struct genotypes run."""
    prefix_path = tmp_path / "hive"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "-vvvv",
            "struct",
            "-i",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "genotypes",
            "-p",
            str(prefix_path),
        ],
    )

    assert result.exit_code == 0


def test_struct_genotypes_threads(tmp_path: pathlib.Path) -> None:
    """Basic struct genotypes run."""
    prefix_path = tmp_path / "hive"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "-vvvv",
            "-t",
            "4",
            "struct",
            "-i",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "genotypes",
            "-p",
            str(prefix_path),
        ],
    )

    assert result.exit_code == 0


def test_annotations_vcf(tmp_path: pathlib.Path) -> None:
    """Basic annotations vcf run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "annotations",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-o",
            str(annotations_path),
            "vcf",
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(annotations_path)
    assert lf.columns == [
        "vid",
        "id",
        "AF_ESP",
        "AF_EXAC",
        "AF_TGP",
        "ALLELEID",
        "CLNDN",
        "CLNDNINCL",
        "CLNDISDB",
        "CLNDISDBINCL",
        "CLNHGVS",
        "CLNREVSTAT",
        "CLNSIG",
        "CLNSIGCONF",
        "CLNSIGINCL",
        "CLNVC",
        "CLNVCSO",
        "CLNVI",
        "DBVARID",
        "GENEINFO",
        "MC",
        "ORIGIN",
        "RS",
    ]


def test_annotations_vcf_not_vcf(tmp_path: pathlib.Path) -> None:
    """Not a vcf input Basic annotations vcf run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "annotations",
            "-i",
            str(DATA_DIR / "no_info.tsv"),
            "-o",
            str(annotations_path),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "vcf",
        ],
    )

    assert result.exit_code == 21


def test_annotations_vcf_select(tmp_path: pathlib.Path) -> None:
    """Select some annotations vcf run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "annotations",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "-o",
            str(annotations_path),
            "vcf",
            "-i",
            "AF_ESP",
            "-i",
            "CLNSIG",
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(annotations_path)
    assert lf.columns == ["vid", "id", "AF_ESP", "CLNSIG"]


def test_annotations_vcf_select_rename_id(tmp_path: pathlib.Path) -> None:
    """Select some annotations vcf run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "-vvvv",
            "annotations",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-o",
            str(annotations_path),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "vcf",
            "-r",
            "annot_id",
            "-i",
            "AF_ESP",
            "-i",
            "CLNSIG",
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(annotations_path)
    assert set(lf.columns) == {"annot_id", "id", "AF_ESP", "CLNSIG"}


def test_annotations_csv(tmp_path: pathlib.Path) -> None:
    """Basic annotations csv run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "annotations",
            "-i",
            str(DATA_DIR / "annotations.csv"),
            "-o",
            str(annotations_path),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "csv",
            "-c",
            "chr",
            "-p",
            "pos",
            "-r",
            "ref",
            "-a",
            "alt",
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(annotations_path)
    assert lf.columns == [
        "ALLELEID",
        "ENUM",
        "id",
    ]


def test_annotations_csv_select(tmp_path: pathlib.Path) -> None:
    """Select some annotations csv run."""
    annotations_path = tmp_path / "annotations.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "annotations",
            "-i",
            str(DATA_DIR / "annotations.csv"),
            "-o",
            str(annotations_path),
            "-c",
            str(DATA_DIR / "grch38.92.csv"),
            "csv",
            "-c",
            "chr",
            "-p",
            "pos",
            "-r",
            "ref",
            "-a",
            "alt",
            "-i",
            "ENUM",
        ],
    )

    assert result.exit_code == 0

    lf = polars.scan_parquet(annotations_path)
    assert lf.columns == [
        "ENUM",
        "id",
    ]


def test_metadata_json(tmp_path: pathlib.Path) -> None:
    """Run metadata json."""
    metadata_path = tmp_path / "metadata.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "metadata",
            "-i",
            str(DATA_DIR / "metadata.json"),
            "-o",
            str(metadata_path),
            "json",
        ],
    )

    assert result.exit_code == 0

    polars.testing.assert_frame_equal(
        polars.scan_parquet(metadata_path),
        polars.scan_parquet(DATA_DIR / "metadata.parquet"),
        check_row_order=False,
    )


def test_metadata_jsonl(tmp_path: pathlib.Path) -> None:
    """Run metadata json."""
    metadata_path = tmp_path / "metadata.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "metadata",
            "-i",
            str(DATA_DIR / "metadata.jsonl"),
            "-o",
            str(metadata_path),
            "json",
            "-l",
        ],
    )

    assert result.exit_code == 0

    polars.testing.assert_frame_equal(
        polars.scan_parquet(metadata_path),
        polars.scan_parquet(DATA_DIR / "metadata.parquet"),
        check_row_order=False,
    )


def test_metadata_json_select(tmp_path: pathlib.Path) -> None:
    """Run metadata json."""
    metadata_path = tmp_path / "metadata.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "metadata",
            "-i",
            str(DATA_DIR / "metadata.json"),
            "-o",
            str(metadata_path),
            "json",
            "-f",
            "sample",
            "-f",
            "link",
            "-f",
            "gender",
            "-f",
            "kindex",
        ],
    )

    assert result.exit_code == 0

    truth = polars.scan_parquet(DATA_DIR / "metadata.parquet").select(["sample", "link", "gender", "kindex"])
    value = polars.scan_parquet(metadata_path)

    polars.testing.assert_frame_equal(truth, value, check_row_order=False)


def test_metadata_csv(tmp_path: pathlib.Path) -> None:
    """Run metadata csv."""
    metadata_path = tmp_path / "metadata.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "metadata",
            "-i",
            str(DATA_DIR / "metadata.csv"),
            "-o",
            str(metadata_path),
            "csv",
        ],
    )

    assert result.exit_code == 0

    polars.testing.assert_frame_equal(
        polars.scan_parquet(metadata_path),
        polars.scan_parquet(DATA_DIR / "metadata.parquet"),
    )


def test_metadata_csv_select(tmp_path: pathlib.Path) -> None:
    """Run metadata csv."""
    metadata_path = tmp_path / "metadata.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "metadata",
            "-i",
            str(DATA_DIR / "metadata.csv"),
            "-o",
            str(metadata_path),
            "csv",
            "-c",
            "sample",
            "-c",
            "link",
            "-c",
            "gender",
            "-c",
            "kindex",
        ],
    )

    assert result.exit_code == 0

    truth = polars.scan_parquet(DATA_DIR / "metadata.parquet").select(["sample", "link", "gender", "kindex"])
    value = polars.scan_parquet(metadata_path)

    polars.testing.assert_frame_equal(truth, value, check_row_order=False)


def test_generate_transmission_ped(tmp_path: pathlib.Path) -> None:
    """Run generate origin."""
    transmission_path = tmp_path / "transmission_path.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "generate",
            "transmission",
            "-i",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "-p",
            str(DATA_DIR / "sample.ped"),
            "-t",
            str(transmission_path),
        ],
    )

    assert result.exit_code == 0

    truth = polars.read_parquet(DATA_DIR / "transmission.parquet").sort("id")
    value = polars.read_parquet(transmission_path).sort("id")

    polars.testing.assert_frame_equal(truth, value)


def test_generate_transmission_no_ped(tmp_path: pathlib.Path) -> None:
    """Run generate origin."""
    transmission_path = tmp_path / "transmission_path.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "generate",
            "transmission",
            "-i",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "-I",
            "sample_1",
            "-m",
            "sample_3",
            "-f",
            "sample_2",
            "-t",
            str(transmission_path),
        ],
    )

    assert result.exit_code == 0

    truth = polars.read_parquet(DATA_DIR / "transmission.parquet").sort("id")
    value = polars.read_parquet(transmission_path).sort("id")

    polars.testing.assert_frame_equal(truth, value)


def test_generate_transmission_nothing(tmp_path: pathlib.Path) -> None:
    """Run generate origin."""
    transmission_path = tmp_path / "transmission_path.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "generate",
            "transmission",
            "-i",
            str(DATA_DIR / "no_info.genotypes.parquet"),
            "-t",
            str(transmission_path),
        ],
    )

    assert result.exit_code == 31
