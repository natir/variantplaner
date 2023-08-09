# VariantPlaner

[![ci](https://github.com/natir/variantplaner/actions/workflows/ci.yml/badge.svg)](https://github.com/natir/variantplaner/actions/workflows/ci.yml)
[![doc](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://natir.github.io/variantplaner/)
[![pypi version](https://img.shields.io/pypi/v/variantplaner.svg)](https://pypi.org/project/variantplaner/)

A toolkit to manage many variants from many samples, with limited resources.

`variantplaner` was initially built for a project concerning the human genome, but is now evolving to support any organism. This feature is planned for version 0.2.0.

## Installation

With `pip`:

```bash
pip install git+https://github.com/natir/variantplaner.git@0.1.0#egg=variantplaner
```

With [`pipx`](https://github.com/pipxproject/pipx):

```bash
python -m pip install --user pipx
pipx install git+https://github.com/natir/variantplaner.git@0.1.0#egg=variantplaner
```

## Usage

This section presents basic usage. For a more complete exemple checkout our [usage page](https://natir.github.io/variantplaner/usages/).

!!! warning
    `variantplaner` doesn't support compressed VCFs. This is a downstream trouble we are aware of and sorry about.

### Extract data from one vcf into several parquet files

With `variantplaner`, you can parse an input VCF and save the relevant data into several parquet files.

```bash
variantplaner vcf2parquet -i input.vcf -v variants.parquet -g genotypes.parquet -a annotations.parquet
```

- `-g` option isn't mandatory. If not set you will lose genotyping information, and if `GT` field is present in the input VCF then only heterozygote or homozygote variants will be kept.
- `-a` option isn't mandatory. If not set you will lose "INFO" fields information.


Genotypes encoding:

| `gt` field in parquet file | Meaning                                         |
| -------------------------- | ----------------------------------------------- |
| 0                          | variant not present                             |
| 1                          | heterozygote                                    |
| 2                          | homozygote                                      |
| 3                          | no information (only used in transmission file) |

### Convert parquet files back to vcf

```bash
variantplaner parquet2vcf -i variants.parquet -g genotypes.parquet -o output.vcf
```

`-g` option isn't mandatory if not set the information isn't added.
This options has many options that control the behavior of this subcommand, we apologize for this complexity.

### Structuration of data

#### Merge variants

!!! danger
    This command can have huge memory and disk usage

```bash
variantplaner struct -i variants/1.parquet -i variants/2.parquet -i variants/3.parquet … -i variants/n.parquet variants -o variants.parquet
```

???+ tip
    By default temporary files are written to /tmp, but you can set your `TMPDIR`, `TEMP` or `TMP` environment variables to change this behavior.

This command uses the [divide-and-conquer algorithm](https://en.wikipedia.org/wiki/Divide-and-conquer_algorithm) to perform variants merging. The `-b|--bytes-memory-limit` option controls the size (in bytes) of each file chunk. Empirically RAM usage will be ten times this limit.

#### Partitioning genotypes

???+ danger
    This command can have huge disk usage

```bash
variantplaner struct -i genotypes/1.parquet -i genotypes/2.parquet -i genotypes/3.parquet … -i genotypes/n.parquet genotypes -p partition_prefix/
```

### Annotations

You can export annotation fields from VCF or CSV/TSV files into parquet files.
To do so, variantplaner provides the `annotations` subcommand.

Command:
```bash
variantplaner annotations -i $INPUT_FILE -o $OUTPUT_PARQUET $INPUT_TYPE [OPTIONS...]
```

Where:

- `-i|--input-path` is the input file (required)
- `-o|--output-path` is the output parquet file (required)
- `$INPUT_TYPE` is whether VCF or CSV (see below for the different value types)

Following OPTIONS are input type-specific (see below).

#### VCF format

If you wish to export `CLNDN` and `AF_ESP` fields from `annotations.vcf` into `clinvar.parquet`, you can run the following command:

```bash
variantplaner annotations -i annotations.vcf -o clinvar.parquet vcf -r annot_id --info CLNDN --info AF_ESP
```

`clinvar.parquet` will contain `id` of variant as well as all the info fields you've selected with the `info` option.
If not set, all the info columns will end up in the output file.

Options:

- `-r|--rename-id`: Can be used to rename vcf id column name (default is `vid`).
- `-i|--info`: Lets you select the info fields you wish to output. If not set, this will export them all.
- `vcf`: If the input file type is VCF

!!! tip
    Mind the `vcf` argument, as the following options depend on the input file type.

#### CSV or TSV format

```bash
variantplaner annotations -i annotations.tsv -o annotations.parquet csv -c chr -p pos -r ref -a alt -s$'\t' --info CLNDN --info AF_ESP
```

Unlike the VCF format, `variantplaner` has no way to tell which columns in the CSV/TSV file correspond to the relevant fields of a variant file.
This is why you need to specify the column names in the options (requires a header).

Options:

- `-c|--chromosome`: Name of chromosome column
- `-p|--position`: Name of position column
- `-r|--reference`: Name of reference column
- `-a|--alternative`: Name of alternative column
- `-i|--info`: Lets you select the info fields you'd like to output. If not set, this will export them all.
- `-s|--separator`: A single byte character to use as a delimiter in the input file (defaults to `,`)

### Metadata

#### JSON format

```bash
variantplaner metadata -i metadata.json -o metadata.parquet json -f sample -f link -f kindex
```

#### Csv format

```bash
variantplaner metadata -i metadata.csv -o metadata.parquet csv -c sample -c link -c kindex
```

### Generate

#### Variants transmission

If you study germline variants it's useful to calculate the familial origin of variants.

```bash
variantplaner generate transmission -i genotypes.parquet -I index_sample_name -m mother_sample_name -f father_sample_name -t transmission.parquet
```

`genotypes.parquet` file with variants of all family. This file must contains `gt` and `samples` columns.

In `transmission.parquet` each line contains an index sample variants, index, mother, father genotypes sample information and also column origin.

Origin column contains a number with 3 digit:
```
231
││└ father genotype
│└─ mother genotype
└── index genotype
```

You can also use pedigree file:

```bash
variantplaner generate transmission -i genotypes.parquet -p family.ped -t transmission.parquet
```

???+ danger
	This command could have important RAM usage (propotionaly to number of sample index variants)

## Contribution

All contributions are welcome, [see our "How to contribute" page.](https://natir.github.io/variantplaner/contributing/)
