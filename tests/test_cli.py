"""Tests for cli."""

# std import
from __future__ import annotations

import filecmp
import pathlib

# 3rd party import
import polars
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


MERGE_IDS = [
    14242097786807219277,
    11061441912435479070,
    15556898136537930176,
    11711297912729185188,
    10807139693934085530,
    11754354483181031149,
    4797423062503056673,
    12330835685209145,
    3370214193463694380,
    13133999947024542022,
    800025099629731040,
    9722743600114191928,
    9554065504802534594,
    2295086632353399847,
    8068220590006330225,
    14452059200390301692,
    1062943517767600077,
    3925979271274072158,
    15044799663914522374,
    16588824010281864939,
    1485127797546295921,
    14498991407705406549,
    4455610085260365294,
    6760366560597639209,
    17246148306088235974,
    270423859963006067,
    11952270119407113547,
    10499890489289074188,
    11976298723311757531,
    721273245970054755,
    7645340256260622458,
    11033100074712141168,
    14319454638398535779,
    17224330959165379425,
    11920660801823307015,
    3585291275407140934,
    18172124183201916611,
    10212290916779149099,
    4759541468123016608,
    11650605831284591550,
    12945673213782943443,
    9664051898850350365,
    2373425740527121342,
    11903476464390241143,
    10487392163259126218,
    5356120651941363990,
    4235309537048834275,
    16247809233398557031,
]


def test_vcf2parquet(tmp_path: pathlib.Path) -> None:
    """vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["vcf2parquet", "-i", str(DATA_DIR / "no_info.vcf"), "-v", str(variants_path), "-g", str(genotypes_path)],
    )

    assert result.exit_code == 0
    assert filecmp.cmp(DATA_DIR / "no_info.variants.parquet", variants_path)
    assert filecmp.cmp(DATA_DIR / "no_info.genotypes.parquet", genotypes_path)


def test_vcf2parquet_not_ask_genotypes(tmp_path: pathlib.Path) -> None:
    """Not ask genotype vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["vcf2parquet", "-i", str(DATA_DIR / "no_info.vcf"), "-v", str(variants_path)],
    )

    assert result.exit_code == 0
    assert filecmp.cmp(DATA_DIR / "no_info.variants.parquet", variants_path)


def test_vcf2parquet_not_vcf(tmp_path: pathlib.Path) -> None:
    """Not a vcf input vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["vcf2parquet", "-i", str(DATA_DIR / "no_info.tsv"), "-v", str(variants_path), "-g", str(genotypes_path)],
    )

    assert result.exit_code == 1


def test_vcf2parquet_no_genotype(tmp_path: pathlib.Path) -> None:
    """No genotype in vcf vcf2parquet run."""
    variants_path = tmp_path / "variants.parquet"
    genotypes_path = tmp_path / "genotypes.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        ["vcf2parquet", "-i", str(DATA_DIR / "no_genotypes.vcf"), "-v", str(variants_path), "-g", str(genotypes_path)],
    )

    assert result.exit_code == 2


def test_parquet2vcf(tmp_path: pathlib.Path) -> None:
    """parquet2vcf run."""
    variants_path = tmp_path / "variants.vcf"

    runner = CliRunner()
    result = runner.invoke(cli.main, ["parquet2vcf", "-i", str(DATA_DIR / "no_info.parquet"), "-o", str(variants_path)])

    assert result.exit_code == 0

    assert filecmp.cmp(variants_path, DATA_DIR / "no_info.parquet2vcf.vcf")


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

    assert lf.collect().get_column("id").to_list() == MERGE_IDS


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
            "vcf",
        ],
    )

    assert result.exit_code == 1


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
            "annotations",
            "-i",
            str(DATA_DIR / "no_genotypes.vcf"),
            "-o",
            str(annotations_path),
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
    assert lf.columns == ["annot_id", "id", "AF_ESP", "CLNSIG"]


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

    assert filecmp.cmp(metadata_path, DATA_DIR / "metadata.parquet")


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

    polars.testing.assert_frame_equal(truth, value)


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

    assert filecmp.cmp(metadata_path, DATA_DIR / "metadata.parquet")


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

    polars.testing.assert_frame_equal(truth, value)


def test_generate_origin(tmp_path: pathlib.Path) -> None:
    """Run generate origin."""

    transmission_path = tmp_path / "transmission_path.parquet"

    runner = CliRunner()
    result = runner.invoke(
        cli.main,
        [
            "generate",
            "origin",
            "-i",
            DATA_DIR / "no_info.genotypes.parquet",
            "-p",
            DATA_DIR / "sample.ped",
            "-t",
            transmission_path,
        ]
    )

    assert result.exit_code == 0
