# VariantPlaner

variantplaner is a set of tools that converts a large set of vcf's into an interoperable data structure in an efficient way.

To show the capabilities of the variant planner we will use a small example.

## Setup

This tutorial assume you are on unix like system, you have python setup and you [install variantplaner](https://github.com/natir/variantplaner#installation)


Requirements list:

- curl
- gunzip

Optional:

- [pqrs](https://github.com/manojkarthick/pqrs)
- gnu-parallel


## Download data

```bash
mkdir -p vp_tuto/vcf/
cd vp_tuto
URI_ROOT="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"
curl ${URI_ROOT}/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG001.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG002.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG003.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG004.vcf
curl ${URI_ROOT}/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG005.vcf
curl ${URI_ROOT}/ChineseTrio/HG006_NA24694_father/latest/GRCh38/HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG006.vcf
curl ${URI_ROOT}/ChineseTrio/HG007_NA24695_mother/latest/GRCh38/HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG007.vcf
```

## Variant planner presentation

variantplaner is a python module with command line tools, it is composed of several subcommands (they will be detailed later) but has two global options, one for parallelization and another for the level of verbosity you want.
```
Usage: variantplaner [OPTIONS] COMMAND [ARGS]...

  Run VariantPlanner.

Options:
  -t, --threads INTEGER RANGE  Number of threads usable  [default: 1; x>=0]
  -v, --verbose                Verbosity level  [0<=x<=4]
  --help                       Show this message and exit.

Commands:
  annotations  Convert an annotation variation file in parquet.
  metadata     Convert an metadata file in parquet.
  parquet2vcf  Convert parquet in vcf.
  struct       Struct operation on parquet file.
  vcf2parquet  Convert a vcf in parquet.
```


## Vcf2Parquet

First step is convert vcf data in [parquet](https://en.wikipedia.org/wiki/Apache_Parquet) it's a column oriented format with better performance.

We split vcf in two part on for variant information and another for genotype information.

```bash
mkdir -p variants genotypes
```

```
for vcf_path in $(ls vcf/*.vcf)
do
vcf_basename=$(basename ${vcf_path} .vcf)
variantplaner -t 4 vcf2parquet -i ${vcf_path} \
-v variants/${vcf_basename}.parquet \
-g genotypes/${vcf_basename}.parquet
done
```

/// details | gnu-parallel method
```bash
find vcf -type f -name *.vcf -exec basename {} .vcf \; | \
parallel variantplaner -t 2 vcf2parquet -i vcf/{}.vcf -v variants/{}.parquet \
-g genotypes/{}.parquet -f GT:PS:DP:ADALL:AD:GQ
```
///

Parquet variants file contains 5 column:

- chr: Chromosome name, X -> 22, Y -> 23, MT -> 24
- pos: Position of variant
- ref: Reference sequence
- alt: Alternative sequence
- id: An hash of other value colision isn't check but highly improbable [check api documentation](/variantplaner/reference/variantplaner/normalization/#variantplaner.normalization.add_variant_id)


/// details | variants parquet file content
You can inspect content of parquet file generate with pqrs
```bash
pqrs head variants/HG001.parquet
{id: 17886044532216650390, chr: 1, pos: 783006, ref: "A", alt: "G"}
{id: 7513336577790240873, chr: 1, pos: 783175, ref: "T", alt: "C"}
{id: 17987040642944149052, chr: 1, pos: 784860, ref: "T", alt: "C"}
{id: 10342734968077036194, chr: 1, pos: 785417, ref: "G", alt: "A"}
{id: 890514037559296207, chr: 1, pos: 797392, ref: "G", alt: "A"}
```
///

Parquet genotypes file contains column:

- id: Same as variant id
- gt: vcf GT value 1 -> heterozygote 2 -> homozygote (phasing information is lost)
- ps: Phase set in which this variant falls
- dp: vcf DP coverage of the variant for this sample
- adall: Net allele depths across all datasets
- ad: vcf AD per allele reads depth
- gq: vcf GQ quality of variant for this sample

/// details | genotypes parquet file content
You can inspect content of parquet file generate with pqrs
```bash
pqrs head genotypes/HG001.parquet
{id: 17886044532216650390, sample: "HG001", gt: 2, ps: null, dp: 652, adall: [16, 234], ad: [0, 82], gq: 312}
{id: 7513336577790240873, sample: "HG001", gt: 2, ps: null, dp: 639, adall: [0, 218], ad: [0, 84], gq: 194}
{id: 17987040642944149052, sample: "HG001", gt: 2, ps: null, dp: 901, adall: [105, 406], ad: [0, 74], gq: 301}
{id: 10342734968077036194, sample: "HG001", gt: 2, ps: null, dp: 820, adall: [125, 383], ad: [0, 70], gq: 339}
{id: 890514037559296207, sample: "HG001", gt: 1, ps: null, dp: 760, adall: [161, 142], ad: [25, 37], gq: 147}
```
///

## Structuration of data

### Merge all variant

We can now aggregate all variant present in our dataset to perform this operation we use divide to conquer merge methode by generate temporary file. By default file are write in `/tmp` but you can control where these files are write by set `TMPDIR`, `TEMP` or `TMP` directory.

```bash
input=$(find variants -type f -name *.parquet -exec echo "-i" {} \; | tr '\n' ' ')
variantplaner -t 8 struct ${input} variants -o variants.parquet
```

File `variants.parquet` containt all uniq variants present in dataset.

### Genotypes structuration

```bash
input=$(find variants -type f -name *.parquet -exec echo "-i" {} \; | tr '\n' ' ')
variantplaner -t 8 struct ${input} genotypes -p hive/
```

All genotypes information are split in an [hive structure](https://duckdb.org/docs/data/partitioning/hive_partitioning) to optimize request on data.

## Add annotations

To work on your variant you probably need and annotations.

### Snpeff annotations

First convert your unique variants in parquet format (`variants.parquet`) in vcf:
```bash
variantplaner -t 8 parquet2vcf -i variants.parquet -o variants.vcf
```

`parquet2vcf` subcommand have many more option but we didn't need it now.

Next annotate this `variants.vcf` [with snpeff](https://pcingola.github.io/SnpEff/), we assume you generate a file call `variants.snpeff.vcf`.

To convert annotated vcf in parquet, keep 'ANN' info column and rename vcf id column in snpeff\_id you can run:
```bash
variantplaner -t 8 annotations -i variants.snpeff.vcf -o snpeff.parquet vcf -i ANN -r snpeff_id
```

If you didn't set any value of option `-i` in vcf subsubcommand all info column are keep.

### Clinvar annotations

Download last clinvar version:

```bash
mkdir annotations
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz \
| gunzip - > annotations/clinvar.vcf
```

Convert clinvar vcf file in parquet file:

```bash
variantplaner annotations -i annotations/clinvar.vcf -o annotations/clinvar.parquet vcf -r clinvar_id
```

Parquet file produce contains many columns:

- clinvar_id: content of vcf id column if option `-r` is not set column name is `vid`
- id: variantplaner id
- All INFO filed

Annotations subcommand try to make match between vcf info type and parquet type.

/// details | `annotations/clinvar.parquet` file content
```bash
pqrs head annotations/clinvar.parquet
{clinvar_id: 2205837, id: 11650605831284591550, AF_ESP: null, AF_EXAC: null, AF_TGP: null, ALLELEID: 2193183, CLNDN: ["Inborn_genetic_diseases"], CLNDNINCL: null, CLNDISDB: ["MeSH:D030342", "MedGen:C0950123"], CLNDISDBINCL: null, CLNHGVS: ["NC_000001.11:g.69134A>G"], CLNREVSTAT: ["criteria_provided", "_single_submitter"], CLNSIG: ["Likely_benign"], CLNSIGCONF: null, CLNSIGINCL: null, CLNVC: "single_nucleotide_variant", CLNVCSO: "SO:0001483", CLNVI: null, DBVARID: null, GENEINFO: "OR4F5:79501", MC: ["SO:0001583|missense_variant"], ORIGIN: ["1"], RS: null}
{clinvar_id: 2252161, id: 2295086632353399847, AF_ESP: null, AF_EXAC: null, AF_TGP: null, ALLELEID: 2238986, CLNDN: ["Inborn_genetic_diseases"], CLNDNINCL: null, CLNDISDB: ["MeSH:D030342", "MedGen:C0950123"], CLNDISDBINCL: null, CLNHGVS: ["NC_000001.11:g.69581C>G"], CLNREVSTAT: ["criteria_provided", "_single_submitter"], CLNSIG: ["Uncertain_significance"], CLNSIGCONF: null, CLNSIGINCL: null, CLNVC: "single_nucleotide_variant", CLNVCSO: "SO:0001483", CLNVI: null, DBVARID: null, GENEINFO: "OR4F5:79501", MC: ["SO:0001583|missense_variant"], ORIGIN: ["1"], RS: null}
{clinvar_id: 2396347, id: 11033100074712141168, AF_ESP: null, AF_EXAC: null, AF_TGP: null, ALLELEID: 2386655, CLNDN: ["Inborn_genetic_diseases"], CLNDNINCL: null, CLNDISDB: ["MeSH:D030342", "MedGen:C0950123"], CLNDISDBINCL: null, CLNHGVS: ["NC_000001.11:g.69682G>A"], CLNREVSTAT: ["criteria_provided", "_single_submitter"], CLNSIG: ["Uncertain_significance"], CLNSIGCONF: null, CLNSIGINCL: null, CLNVC: "single_nucleotide_variant", CLNVCSO: "SO:0001483", CLNVI: null, DBVARID: null, GENEINFO: "OR4F5:79501", MC: ["SO:0001583|missense_variant"], ORIGIN: ["1"], RS: null}
{clinvar_id: 2288999, id: 10487392163259126218, AF_ESP: null, AF_EXAC: null, AF_TGP: null, ALLELEID: 2278803, CLNDN: ["Inborn_genetic_diseases"], CLNDNINCL: null, CLNDISDB: ["MeSH:D030342", "MedGen:C0950123"], CLNDISDBINCL: null, CLNHGVS: ["NC_000001.11:g.69769T>C"], CLNREVSTAT: ["criteria_provided", "_single_submitter"], CLNSIG: ["Uncertain_significance"], CLNSIGCONF: null, CLNSIGINCL: null, CLNVC: "single_nucleotide_variant", CLNVCSO: "SO:0001483", CLNVI: null, DBVARID: null, GENEINFO: "OR4F5:79501", MC: ["SO:0001583|missense_variant"], ORIGIN: ["1"], RS: null}
{clinvar_id: 2351346, id: 5356120651941363990, AF_ESP: null, AF_EXAC: null, AF_TGP: null, ALLELEID: 2333177, CLNDN: ["Inborn_genetic_diseases"], CLNDNINCL: null, CLNDISDB: ["MeSH:D030342", "MedGen:C0950123"], CLNDISDBINCL: null, CLNHGVS: ["NC_000001.11:g.69995G>C"], CLNREVSTAT: ["criteria_provided", "_single_submitter"], CLNSIG: ["Uncertain_significance"], CLNSIGCONF: null, CLNSIGINCL: null, CLNVC: "single_nucleotide_variant", CLNVCSO: "SO:0001483", CLNVI: null, DBVARID: null, GENEINFO: "OR4F5:79501", MC: ["SO:0001583|missense_variant"], ORIGIN: ["1"], RS: null}
```
///

With option of subcommand vcf `-i` you can select which column are include in parquet file
For example command:
```bash
variantplaner annotations -i annotations/clinvar.vcf -o annotations/clinvar.parquet vcf -r clinvar_id -i ALLELEID -i CLNDN
```

Produce a `annotations/clinvar.parquet` with columns:

- clinvar_id
- id
- ALLELEID
- CLNDN


/// details | `annotations/clinvar.parquet` file content
```bash
➜ pqrs head annotations/clinvar.parquet
{clinvar_id: 2205837, id: 11650605831284591550, ALLELEID: 2193183, CLNDN: ["Inborn_genetic_diseases"]}
{clinvar_id: 2252161, id: 2295086632353399847, ALLELEID: 2238986, CLNDN: ["Inborn_genetic_diseases"]}
{clinvar_id: 2396347, id: 11033100074712141168, ALLELEID: 2386655, CLNDN: ["Inborn_genetic_diseases"]}
{clinvar_id: 2288999, id: 10487392163259126218, ALLELEID: 2278803, CLNDN: ["Inborn_genetic_diseases"]}
{clinvar_id: 2351346, id: 5356120651941363990, ALLELEID: 2333177, CLNDN: ["Inborn_genetic_diseases"]}
```
///


## Querying
