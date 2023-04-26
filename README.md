# VariantPlanner

[![ci](https://github.com/natir/variantplanner/actions/workflows/ci.yml/badge.svg)](https://github.com/natir/variantplanner/actions/workflows/ci.yml)
[![doc](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://natir.github.io/variantplanner/)
[![pypi version](https://img.shields.io/pypi/v/variantplanner.svg)](https://pypi.org/project/variantplanner/)


A tool kit to manage many variant on desktop computer

## Installation

With `pip`:
```bash
pip install git+https://github.com/natir/variantplanner.git#egg=variantplanner
```

With [`pipx`](https://github.com/pipxproject/pipx):
```bash
python -m pip install --user pipx
pipx install git+https://github.com/natir/variantplanner.git#egg=variantplanner
```

## Usage

### Convert vcf in parquet

```
variantplanner vcf2parquet -i input.vcf -v variants.parquet -g genotypes.parquet
```

Convert multiple vcf in parquet

```
mkdir -p variants genotypes
for path in $(ls vcf/*.vcf)
do
vcf_basename=$(basename ${path} .vcf)
variantplanner vcf2parquet -i vcf/${vcf_basename}.vcf -v variants/${vcf_basename}.parquet -g genotypes/${vcf_basename}.parquet
done
```

or use gnu parallel

```
find tests/data/ -type f -name *.vcf -exec basename {} .vcf \; | parallel variantplanner vcf2parquet -i vcf/{}.vcf -v variants/{}.parquet -g genotypes/{}.parquet
```

### Structuration of data

#### Merge variants

**Warning**: this command could have huge memory and disk usage

```
variantplanner struct -i variants/1.parquet -i variants/2.parquet -i variants/3.parquet … -i variants/n.parquet variants -o variants.parquet
```

By default temporary file are write in /tmp you can use TMPDIR, TEMP or TMP to change this behavior.

This command use divide and conquer algorithm to perform merge of variants option `-b|--bytes-memory-limit` control bytes size of chunk of files.

#### Hive genotypes

**Warning**: this command could have huge disk usage

```
variantplanner struct -i genotypes/1.parquet -i genotypes/2.parquet -i genotypes/3.parquet … -i genotypes/n.parquet genotypes -p hive_prefix/
```

### Annotations

#### Vcf format

```
variantplanner annotations -i annotations.vcf -o annotations.parquet vcf -r annot_id --info CLNDN --info AF_ESP
```

`clinvar.parquet` containts `id` of variant and info field select if you didn't set `info` option all info column are include.

Option `-r|--rename-id` could be use to rename vcf id column name (default name is `vid`).

#### Csv format

```
variantplanner annotations -i annotations.tsv -o annotations.parquet csv -c chr -p pos -r ref -a alt -s$'\t' --info CLNDN --info AF_ESP
```

It's work same as vcf sub command but you must specify chromosome (`-c`), position (`-p`), reference (`-r`) and alternative (`-a`), you can change separate with option `-s`

### Metadata

#### Json format

```
variantplanner metadata -i metadata.json -o metadata.parquet json -f sample -f link -f kindex
```

#### Csv format

```
variantplanner metadata -i metadata.csv -o metadata.parquet csv -c sample -c link -c kindex
```

## Contribution

All contributions are welcome, [see our "How to contribute" page.](https://natir.github.io/variantplanner/contributing/)
