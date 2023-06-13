"""Script to generate the benchmark plots."""

from __future__ import annotations

import json

# std import
import os
import pathlib
import typing
from textwrap import dedent

# 3rd party import
import altair
import polars
from jinja2 import StrictUndefined
from jinja2.sandbox import SandboxedEnvironment

if typing.TYPE_CHECKING:  # pragma: no cover
    import pandas

# project import


def last_benchark() -> pathlib.Path:
    """Found last benchmark data."""
    change_entry = []
    for entry in os.scandir(".benchmarks"):
        if entry.is_dir():
            for subentry in os.scandir(entry.path):
                if subentry.is_file():
                    change_entry.append((subentry.stat().st_mtime, subentry.path))

    return pathlib.Path(sorted(change_entry)[-1][1])


def nothing(data: pandas.DataFrame, name: str) -> tuple[pandas.DataFrame, str]:
    """Do nothing."""
    return data, name


def annotations_vcf_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update dataframe and group name for annotations vcf group."""
    rename = {
        "one": 1,
        "two": 2,
        "three": 3,
        "four": 4,
        "five": 5,
        "six": 6,
        "seven": 7,
        "height": 8,
        "nine": 9,
        "ten": 10,
        "eleven": 11,
        "twelve": 12,
        "thirteen": 13,
        "fourteen": 14,
        "fiveteen": 15,
        "sixteen": 16,
        "seventeen": 17,
        "eighteen": 18,
        "nineteen": 19,
        "twenty": 20,
        "full": 21,
        "full_rename": 22,
    }

    data = data.sort_values(by=["method"], key=lambda col: col.apply(lambda x: rename[x]))

    return data, "Effect of number of columns use in annotations conversion"


def vcf2parquet_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update dataframe and group name for vcf2parquet group."""
    data = data.sort_values(by=["method"], key=lambda col: col.apply(lambda x: len(x)))

    return data, "Vcf2Parquet mode evaluation"


def create_plot(
    data: pandas.DataFrame,
    name: str,
    name2func: typing.Mapping[str, typing.Callable[[pandas.DataFrame, str], tuple[pandas.DataFrame, str]]],
) -> altair.Chart:
    """Create plot."""
    df, title = name2func[name](data.to_pandas(), name)

    df["iqr0"] = df["median"] - df["iqr"]
    df["iqr1"] = df["median"] + df["iqr"]

    line = (
        altair.Chart(df)
        .mark_line()
        .encode(
            x=altair.X("method", sort=None),
            y=altair.Y("median", sort=None),
        )
    )

    band = (
        altair.Chart(df)
        .mark_errorband()
        .encode(
            altair.Y(
                "iqr1:Q",
                scale=altair.Scale(zero=False),
                title="Times in (s)",
            ),
            altair.Y2("iqr0:Q"),
            altair.X("method", sort=None),
        )
        .properties(
            title=title,
        )
    )

    return band + line


def render_plot() -> str:
    """Generate benchmark plot."""
    with open(last_benchark()) as fh:
        json_data = json.load(fh)["benchmarks"]

    bench_data: typing.Mapping[str, list] = {"benchmark": [], "method": [], "median": [], "iqr": []}
    for obj in json_data:
        bench_data["benchmark"].append(obj["group"])
        bench_data["method"].append(obj["name"])
        bench_data["median"].append(obj["stats"]["median"])
        bench_data["iqr"].append(obj["stats"]["iqr"])

    df = polars.DataFrame(bench_data)

    name2func = {
        name: globals()[f"{name}_func"] if f"{name}_func" in globals() else nothing
        for name in df.get_column("benchmark").unique().to_list()
    }

    bench2plot = {}
    for name, data in df.groupby("benchmark"):
        bench2plot[name] = create_plot(data, str(name), name2func)

    template_data = {"bench2plot": bench2plot}

    template_text = dedent(
        """
        <script src="https://cdn.jsdelivr.net/npm/vega@5.25.0"></script>
        <script src="https://cdn.jsdelivr.net/npm/vega-lite@5.9.1"></script>
        <script src="https://cdn.jsdelivr.net/npm/vega-embed@6.22.1"></script>

        <div id="vis"></div>

        {% for bench in bench2plot|sort() -%}
        <div id="{{ bench }}"></div>
        <script>
        vegaEmbed(
          "#{{ bench }}",
          {{ bench2plot[bench].to_json() }},
        )
        </script>
        <br>
        {% endfor %}
        """,
    )

    jinja_env = SandboxedEnvironment(undefined=StrictUndefined)
    return jinja_env.from_string(template_text).render(**template_data)


print(render_plot())
