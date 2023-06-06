# VariantPlaner

[![ci](https://github.com/natir/variantplaner/actions/workflows/ci.yml/badge.svg)](https://github.com/natir/variantplaner/actions/workflows/ci.yml)
[![doc](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://natir.github.io/variantplaner/)
[![pypi version](https://img.shields.io/pypi/v/variantplaner.svg)](https://pypi.org/project/variantplaner/)

A tool kit to manage many variant without many cpu and ram ressource.

## Installation

With `pip`:
```bash
pip install git+https://github.com/natir/variantplaner.git@0.1.0-alpha#egg=variantplaner
```

With [`pipx`](https://github.com/pipxproject/pipx):
```bash
python -m pip install --user pipx
pipx install git+https://github.com/natir/variantplaner.git@0.1.0-alpha#egg=variantplaner
```

## Usage

This section present basic usage for a more complete exemple check our [usage page](https://natir.github.io/variantplaner/usages/)

**WARNING**: variantplaner support only not compressed file, sorry it's a downstream trouble.

### Convert vcf in parquet

```
variantplaner vcf2parquet -i input.vcf -v variants.parquet -g genotypes.parquet -a annotations.parquet
```

`-g` option isn't mandatory if you didn't set it you lose genotyping information.
`-a` option isn't mandatory if you didn't set it you lose "INFO" fields information.


Genotyping encoding:

| `gt` parquet value | translation |
| --- | --- |
| 0  | variants not present |
| 1  | heterozygote |
| 2  | homozygote |
| 3  | no information (only use in transmission file) |

### Convert parquet in vcf

```
variantplaner parquet2vcf -i variants.parquet -g genotypes.parquet -o output.vcf
```

`-g` option isn't mandatory if you didn't set it information isn't add.
This options have many options to control behavior of this subcommand, sorry for this complexity.

### Structuration of data

#### Merge variants

**Warning**: this command could have huge memory and disk usage

```
variantplaner struct -i variants/1.parquet -i variants/2.parquet -i variants/3.parquet … -i variants/n.parquet variants -o variants.parquet
```

By default temporary file are write in /tmp you can use TMPDIR, TEMP or TMP to change this behavior.

This command use divide and conquer algorithm to perform merge of variants option `-b|--bytes-memory-limit` control bytes size of chunk of files. Empirically ram usage is ten times bytes memory limit value.

#### Partitioning genotypes

**Warning**: this command could have huge disk usage

```
variantplaner struct -i genotypes/1.parquet -i genotypes/2.parquet -i genotypes/3.parquet … -i genotypes/n.parquet genotypes -p partition_prefix/
```

### Annotations

#### Vcf format

```
variantplaner annotations -i annotations.vcf -o annotations.parquet vcf -r annot_id --info CLNDN --info AF_ESP
```

`clinvar.parquet` containts `id` of variant and info field select if you didn't set `info` option all info column are include.

Option `-r|--rename-id` could be use to rename vcf id column name (default name is `vid`).

#### Csv format

```
variantplaner annotations -i annotations.tsv -o annotations.parquet csv -c chr -p pos -r ref -a alt -s$'\t' --info CLNDN --info AF_ESP
```

It's work same as vcf sub command but you must specify chromosome (`-c`), position (`-p`), reference (`-r`) and alternative (`-a`), you can change separate with option `-s`

### Metadata

#### Json format

```
variantplaner metadata -i metadata.json -o metadata.parquet json -f sample -f link -f kindex
```

#### Csv format

```
variantplaner metadata -i metadata.csv -o metadata.parquet csv -c sample -c link -c kindex
```

### Generate

#### Variants transmission

It is sometimes useful to calculate the familial origin of variants.

```
variantplaner generate transmission -i genotypes.parquet -I index_sample_name -m mother_sample_name -f father_sample_name -t transmission.parquet
```

`genotypes.parquet` file with variant of all family this file must contains `gt` and `samples` columns.

In `transmission.parquet` each line contains an index sample variants, index, mother, father genotypes sample information and also column origin.

Origin column contains a number with 3 digit:
```
231
││└ father genotype
│└─ mother genotype
└── index genotype
```

You can also use pedigree file:
```
variantplaner generate transmission -i genotypes.parquet -p family.ped -t transmission.parquet
```

**Warning**: this command could have important RAM usage (propotionaly to number of sample index variants)


## Contribution

All contributions are welcome, [see our "How to contribute" page.](https://natir.github.io/variantplaner/contributing/)
