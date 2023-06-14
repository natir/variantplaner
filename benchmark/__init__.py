"""Tests suite for `variantplaner`."""

# std import
import os
import random
from pathlib import Path

# 3rd party import

# project import


TESTS_DIR = Path(__file__).parent
TMP_DIR = TESTS_DIR / "tmp"
FIXTURES_DIR = TESTS_DIR / "fixtures"

# Deactivate multi-threading
os.environ["POLARS_MAX_THREADS"] = str(1)


def __generate_vcf(tmp_path: Path, number_of_var: int) -> None:
    """Generate a fake large vcf."""
    random.seed(42)

    # Version
    header = "##fileformat=VCFv4.3\n"

    # Add info
    for info_type in ["Integer", "Float", "Flag", "Character", "String"]:
        for number in ["A", "R", "G", ".", "1", "2", "5"]:
            id_name = f"{info_type}_{number}"
            header += f'##INFO=<ID={id_name},Number={number},Type={info_type},Description="description">\n'

    # Add format
    for format_type in ["Integer", "Float", "Character", "String"]:
        for number in ["A", "R", "G", ".", "1", "2", "5"]:
            format_name = f"{format_type}_{number}"
            header += f'##FORMAT=<ID={format_name},Number={number},Type={info_type},Description="description">\n'

    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts_1\ts_2\n"

    chroms = [*list(range(22)), "Y", "X", "MT"]
    nuc = ["A", "C", "T", "G"]

    with open(tmp_path, "w") as fh:
        print(header, file=fh)
        for _ in range(number_of_var):
            chrom = random.choice(chroms)  # noqa: S311
            pos = random.randint(0, 2 ^ 32)  # noqa: S311
            vid = f"rid_{random.randint(100000, 1000000)}"  # noqa: S311
            ref = "".join(random.choices(nuc, k=random.randint(1, 7)))  # noqa: S311
            alt = "".join(random.choices(nuc, k=random.randint(1, 7)))  # noqa: S311
            info_str = "."
            format_str = "."
            s_1 = "."
            s_2 = "."
            print(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info_str}\t{format_str}\t{s_1}\t{s_2}", file=fh)
