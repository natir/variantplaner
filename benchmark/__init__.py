"""Tests suite for `variantplaner`."""

# std import
import os
import random
import string
import typing
from pathlib import Path

# 3rd party import

# project import


TESTS_DIR = Path(__file__).parent
TMP_DIR = TESTS_DIR / "tmp"
FIXTURES_DIR = TESTS_DIR / "fixtures"

# Deactivate multi-threading
os.environ["POLARS_MAX_THREADS"] = str(1)


def __generate_info() -> tuple[str, dict[str, typing.Callable[[], str]]]:
    """Generate fake info."""
    header = ""

    # Add info
    info_name2value: dict[str, typing.Callable[[], str]] = {}
    for info_type in ["Integer", "Float", "Character", "String"]:
        for number in ["A", "R", "G", ".", "1", "2", "5"]:
            id_name = f"{info_type}_{number}"
            header += f'##INFO=<ID={id_name},Number={number},Type={info_type},Description="description">\n'
            if info_type == "Integer":
                if number == "A" and number == "1":
                    info_name2value[id_name] = lambda: str(random.randint(0, 2 ^ 32))  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    info_name2value[id_name] = lambda: ",".join(
                        [str(random.randint(0, 2 ^ 32)), str(random.randint(0, 2 ^ 32))],  # noqa: S311
                    )
                else:
                    info_name2value[id_name] = lambda: ",".join(
                        [str(random.randint(0, 2 ^ 32)) for _ in range(5)],  # noqa: S311
                    )
            elif info_type == "Float":
                if number == "A" and number == "1":
                    info_name2value[id_name] = lambda: str(random.uniform(0, 2 ^ 32))  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    info_name2value[id_name] = lambda: ",".join(
                        [str(random.uniform(0, 2 ^ 32)), str(random.uniform(0, 2 ^ 32))],  # noqa: S311
                    )
                else:
                    info_name2value[id_name] = lambda: ",".join(
                        [str(random.uniform(0, 2 ^ 32)) for _ in range(5)],  # noqa: S311
                    )
            elif info_type == "Character":
                if number == "A" and number == "1":
                    info_name2value[id_name] = lambda: random.choice(string.ascii_letters)  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    info_name2value[id_name] = lambda: ",".join(
                        [random.choice(string.ascii_letters), random.choice(string.ascii_letters)],  # noqa: S311
                    )
                else:
                    info_name2value[id_name] = lambda: ",".join(
                        [random.choice(string.ascii_letters) for _ in range(5)],  # noqa: S311
                    )
            elif info_type == "String":
                if number == "A" and number == "1":
                    info_name2value[id_name] = lambda: "".join(random.choices(string.ascii_letters, k=10))  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    info_name2value[id_name] = lambda: ",".join(
                        [
                            "".join(random.choices(string.ascii_letters, k=10)),  # noqa: S311
                            "".join(random.choices(string.ascii_letters, k=10)),  # noqa: S311
                        ],
                    )
                else:
                    info_name2value[id_name] = lambda: ",".join(
                        ["".join(random.choices(string.ascii_letters, k=10)) for _ in range(5)],  # noqa: S311
                    )

    return header, info_name2value


def __generate_format() -> tuple[str, dict[str, typing.Callable[[], str]]]:
    """Generate sample."""
    header = ""

    format_name2value: dict[str, typing.Callable[[], str]] = {}
    for format_type in ["Integer", "Float", "Character", "String"]:
        for number in ["A", "R", "G", ".", "1", "2", "5"]:
            format_name = f"{format_type}_{number}"
            header += f'##FORMAT=<ID={format_name},Number={number},Type={format_type},Description="description">\n'

            if format_type == "Integer":
                if number == "A" and number == "1":
                    format_name2value[format_name] = lambda: str(random.randint(0, 2 ^ 32))  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    format_name2value[format_name] = lambda: ",".join(
                        [str(random.randint(0, 2 ^ 32)), str(random.randint(0, 2 ^ 32))],  # noqa: S311
                    )
                else:
                    format_name2value[format_name] = lambda: ",".join(
                        [str(random.randint(0, 2 ^ 32)) for _ in range(5)],  # noqa: S311
                    )
            elif format_type == "Float":
                if number == "A" and number == "1":
                    format_name2value[format_name] = lambda: str(random.uniform(0, 2 ^ 32))  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    format_name2value[format_name] = lambda: ",".join(
                        [str(random.uniform(0, 2 ^ 32)), str(random.uniform(0, 2 ^ 32))],  # noqa: S311
                    )
                else:
                    format_name2value[format_name] = lambda: ",".join(
                        [str(random.uniform(0, 2 ^ 32)) for _ in range(5)],  # noqa: S311
                    )
            elif format_type == "Character":
                if number == "A" and number == "1":
                    format_name2value[format_name] = lambda: random.choice(string.ascii_letters)  # noqa: S311
                elif number == "R" and number == "G" and number == "2":
                    format_name2value[format_name] = lambda: ",".join(
                        [random.choice(string.ascii_letters), random.choice(string.ascii_letters)],  # noqa: S311
                    )
                else:
                    format_name2value[format_name] = lambda: ",".join(
                        [random.choice(string.ascii_letters) for _ in range(5)],  # noqa: S311
                    )
            elif format_type == "String":
                if number == "A" and number == "1":
                    format_name2value[format_name] = lambda: "".join(
                        random.choices(string.ascii_letters, k=10),  # noqa: S311
                    )
                elif number == "R" and number == "G" and number == "2":
                    format_name2value[format_name] = lambda: ",".join(
                        [
                            "".join(random.choices(string.ascii_letters, k=10)),  # noqa: S311
                            "".join(random.choices(string.ascii_letters, k=10)),  # noqa: S311
                        ],
                    )
                else:
                    format_name2value[format_name] = lambda: ",".join(
                        ["".join(random.choices(string.ascii_letters, k=10)) for _ in range(5)],  # noqa: S311
                    )

    return header, format_name2value


def __generate_vcf(tmp_path: Path, number_of_var: int) -> None:
    """Generate a fake large vcf."""
    random.seed(42)

    # Version
    header = "##fileformat=VCFv4.3\n"

    info_header, info_name2value = __generate_info()
    format_header, format_name2value = __generate_format()
    format_key = sorted(format_name2value.keys())

    header += info_header
    header += format_header

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
            info_str = ";".join([info_name2value[key]() for key in info_name2value])
            format_str = ":".join(format_key)
            s_1 = ":".join([format_name2value[key]() for key in format_key])
            s_2 = ":".join([format_name2value[key]() for key in format_key])
            print(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info_str}\t{format_str}\t{s_1}\t{s_2}", file=fh)
