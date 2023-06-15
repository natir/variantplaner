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


def generate_python_version() -> typing.Iterator[tuple[str, str]]:
    """Generate python version available."""
    for entry in os.scandir(".benchmarks"):
        if entry.is_dir() and "Python" in entry.name:
            yield entry.name, entry.name.split("-")[2]


def last_benchark(versiondir: str) -> pathlib.Path:
    """Found last benchmark data."""
    change_entry = []
    for entry in os.scandir(".benchmarks"):
        if entry.is_dir() and entry.name == versiondir:
            for subentry in os.scandir(entry.path):
                if subentry.is_file():
                    change_entry.append((subentry.stat().st_mtime, subentry.path))

    return pathlib.Path(sorted(change_entry)[-1][1])


def nothing(data: pandas.DataFrame, name: str) -> tuple[pandas.DataFrame, str]:
    """Do nothing."""
    return data, name


def annotations_extractions_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update dataframe and group name for annotations extractions group."""
    data["method"] = data["method"].apply(lambda x: int(x.split("_")[2]))
    data = data.sort_values(by=["method"])

    return data, "Effect of number of columns use in annotations conversion"


def vcf2parquet_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update dataframe and group name for vcf2parquet group."""
    data = data.sort_values(by=["method"], key=lambda col: col.apply(lambda x: len(x)))

    return data, "Vcf2Parquet mode evaluation"


def vcf_parsing_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update dataframe and group name for vcf2parquet run group."""
    data["method"] = data["method"].apply(lambda x: int(x.split("_")[2]))
    data = data.sort_values(by=["method"])

    return data, "Vcf parsing scaling"


def merge_id_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update variant merge by id."""
    data["method"] = data["method"].apply(lambda x: int(x.split("_")[3]))
    data = data.sort_values(by=["method"])

    return data, "Merging time by id in number of variant"


def merge_variant_func(data: pandas.DataFrame, _name: str) -> tuple[pandas.DataFrame, str]:
    """Update variant merge by variant."""
    data["method"] = data["method"].apply(lambda x: int(x.split("_")[3]))
    data = data.sort_values(by=["method"])

    return data, "Merging time by variant in number of variant"


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
        .mark_line(point=True)
        .encode(
            x=altair.X("method", sort=None),
            y=altair.Y("median", sort=None),
            color="version",
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
            color="version",
        )
        .properties(
            title=title,
        )
    )

    return band + line


def render_plot() -> str:
    """Generate benchmark plot."""
    bench_data: typing.Mapping[str, list] = {"version": [], "benchmark": [], "method": [], "median": [], "iqr": []}
    for dirname, version in generate_python_version():
        with open(last_benchark(dirname)) as fh:
            json_data = json.load(fh)["benchmarks"]

            for obj in json_data:
                bench_data["version"].append(version)
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
