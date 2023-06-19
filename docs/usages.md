# VariantPlaner

variantplaner is a set of tools that converts a large set of vcf's into an interoperable data structure in an efficient way.

To show the capabilities of the variant planner we will use a small example.

The purpose of this short tutorial is to present:

- how to convert vcf into a more suitable format
- how data can be restructured for querying
- how to integrate variant annotation databases
- how these different files can be used to obtain interesting biological information

This tutorial suggests an organization of files, but you're under no obligation to follow it `variantplaner` is quite flexible in its organization.

## Setup

This tutorial assume you are on unix like system, you have python setup and you [install variantplaner](https://github.com/natir/variantplaner#installation)


Requirements list:

- curl
- gunzip

Optional:

- [pqrs](https://github.com/manojkarthick/pqrs)
- gnu-parallel

Quering dataset:

- [polars-cli](https://crates.io/crates/polars-cli)
- [duckdb](https://duckdb.org/)


## Download data

```bash
mkdir -p vp_tuto/vcf/
cd vp_tuto
URI_ROOT="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"
curl ${URI_ROOT}/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG001.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG002.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG003.vcf
curl ${URI_ROOT}/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | gunzip - > vcf/HG004.vcf
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


## vcf2parquet

First step is convert vcf data in [parquet](https://en.wikipedia.org/wiki/Apache_Parquet) it's a column oriented format with better performance than indexed vcf.

We split vcf in two part on for variant information and another for genotype information.

```bash
mkdir -p variants genotypes/samples/
```

```
for vcf_path in $(ls vcf/*.vcf)
do
    sample_name=$(basename ${vcf_path} .vcf)
    variantplaner -t 4 vcf2parquet -i ${vcf_path} \
    -v variants/${sample_name}.parquet \
    -g genotypes/samples/${sample_name}.parquet \
    -f GT:PS:DP:ADALL:AD:GQ
done
```

We iterate over all vcf, variants are store in `variants/{sample_name}.parquet`, genotype information are store in `variants/{sample_name}.parquet`. Only genotypes with a format value equal to `-f` parameter value are retained.

/// details | gnu-parallel method
```bash
find vcf -type f -name *.vcf -exec basename {} .vcf \; | \
parallel variantplaner -t 2 vcf2parquet -i vcf/{}.vcf -v variants/{}.parquet \
-g genotypes/samples/{}.parquet -f GT:PS:DP:ADALL:AD:GQ
```
///

Parquet variants file contains 5 column:

- chr: Chromosome name, X -> 22, Y -> 23, MT -> 24 (**warning**: this behavior could change in next version)
- pos: Position of variant
- ref: Reference sequence
- alt: Alternative sequence
- id: An hash of other value colision isn't check but highly improbable [check api documentation](/variantplaner/reference/variantplaner/normalization/#variantplaner.normalization.add_variant_id)


/// details | variants parquet file content
You can inspect content of parquet file generate with pqrs
```bash
pqrs head variants/samples/HG001.parquet
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
input=$(ls variants/ | xargs -I {} -x echo "-i variants/"{} | tr '\n' ' ')
variantplaner -t 8 struct $(echo $input) variants -o variants.parquet
```

File `variants.parquet` containt all uniq variants present in dataset.

### Genotypes structuration

### By samples

This structurations data is already down in vcf2parquet step check content of `genotypes/samples`:
```bash
➜ ls genotypes/samples
HG001.parquet  HG002.parquet  HG003.parquet  HG004.parquet  HG006.parquet  HG007.parquet
```

### By variants

Here we'll organize the genotypes information by variants to make it easier to find samples where a variant is present or not.

```bash
mkdir -p genotypes/variants/
input=$(ls genotypes/samples/ | xargs -I {} -x echo "-i genotypes/samples/"{} | tr '\n' ' ')
variantplaner -t 8 struct $(echo $input) genotypes -p genotypes/variants/
```

All genotypes information are split in [hive like structure](https://duckdb.org/docs/data/partitioning/hive_partitioning) to optimize request on data.

### Compute transmission mode

If you are working with families, `variantplaner` can calculate the modes of transmission of the variants.

For these step we need concatenate all genotypes of a AshkenazimTrio in one parquet sample.

```bash
pqrs merge -i genotypes/samples/HG002.parquet genotypes/samples/HG003.parquet genotypes/samples/HG004.parquet -o genotypes/samples/AshkenazimTrio.parquet
mkdir -p genotypes/transmissions/
variantplaner generate transmission -i genotypes/samples/AshkenazimTrio.parquet -I HG002 -m HG003 -f HG004 -t genotypes/transmissions/AshkenazimTrio.parquet
```

`-I` parameter is use for index sample, `-m` parameter is use for mother sample, `-f` parameter is use for father sample only the index sample is mandatory if mother sample or father sample isn't present command work, you could also use a pedigree file with parameter `-p`.

/// detail | transmission parquet file content
```
{id: 17618281883055118663, index_gt: 2, index_ps: null, index_dp: 890, index_adall: [0, 266], index_ad: [0, 415], index_gq: 421, mother_gt: 2, mother_ps: null, mother_dp: 806, mother_adall: [2, 255], mother_ad: [2, 418], mother_gq: 409, father_gt: 1, father_ps: null, father_dp: 885, father_adall: [142, 138], father_ad: [225, 214], father_gq: 547, origin: 221}
{id: 16188392565314286737, index_gt: 1, index_ps: null, index_dp: 1101, index_adall: [175, 146], index_ad: [301, 260], index_gq: 658, mother_gt: 1, mother_ps: null, mother_dp: 944, mother_adall: [125, 138], mother_ad: [250, 260], mother_gq: 555, father_gt: 3, father_ps: null, father_dp: null, father_adall: null, father_ad: null, father_gq: null, origin: 113}
{id: 16258731660316021801, index_gt: 2, index_ps: null, index_dp: 1259, index_adall: [1, 404], index_ad: [65, 637], index_gq: 563, mother_gt: 1, mother_ps: null, mother_dp: 946, mother_adall: [160, 157], mother_ad: [259, 258], mother_gq: 561, father_gt: 2, father_ps: null, father_dp: 1150, father_adall: [0, 353], father_ad: [68, 624], father_gq: 486, origin: 212}
{id: 15289378089603723238, index_gt: 1, index_ps: null, index_dp: 1037, index_adall: [254, 211], index_ad: [36, 34], index_gq: 512, mother_gt: 2, mother_ps: null, mother_dp: 795, mother_adall: [82, 358], mother_ad: [0, 95], mother_gq: 461, father_gt: 1, father_ps: null, father_dp: 831, father_adall: [202, 194], father_ad: [36, 46], father_gq: 546, origin: 121}
{id: 11858617203510458545, index_gt: 1, index_ps: null, index_dp: 1271, index_adall: [224, 187], index_ad: [349, 307], index_gq: 658, mother_gt: 1, mother_ps: null, mother_dp: 1121, mother_adall: [198, 163], mother_ad: [318, 279], mother_gq: 559, father_gt: 3, father_ps: null, father_dp: null, father_adall: null, father_ad: null, father_gq: null, origin: 113}
```
///

Parquet transmissions file contains column all genotypes information with suffix `_index`, `_mother` and `_father` plus a `origin` column

Origin column contains a number with 3 digit:

```
231
││└ father genotype
│└─ mother genotype
└── index genotype
```

In this example case variants is homozygotes in index, mother information is missing, variants is heterozygotes in father.


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
mkdir -p annotations
variantplaner -t 8 annotations -i variants.snpeff.vcf -o annotations/snpeff.parquet vcf -i ANN -r snpeff_id
```

If you didn't set any value of option `-i` in vcf subsubcommand all info column are keep.

### Clinvar annotations

Download last clinvar version:

```bash
mkdir -p annotations
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
variantplaner annotations -i annotations/clinvar.vcf -o annotations/clinvar.parquet vcf -r clinvar_id -i ALLELEID -i CLNDN -i AF_ESP -i GENEINFO
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

You can use any tool or software library supporting the parquet format to use the files generated by `variantplaner`.

We show you how to use files with [polars-cli](https://crates.io/crates/polars-cli) and [duckdb](https://duckdb.org/).

### polars-cli

#### Count variants

```sql
〉select count(*) from read_parquet('variants.parquet');
┌─────────┐
│ count   │
│ ---     │
│ u32     │
╞═════════╡
│ 7852699 │
└─────────┘
```


/// details | check result with `pqrs`
We can check we have same result with pqrs

```bash
➜ pqrs rowcount variants.parquet
File Name: variants.parquet: 7852699 rows
```
///

#### Filter variants from annotations:

Get all variant with a AF_ESP upper than 0.9999

```sql
〉select chr, pos, ref, alt, AF_ESP from read_parquet('variants.parquet') as v left join read_parquet('annotations/clinvar.parquet') as c on c.id=v.id where AF_ESP>0.9999;
┌─────┬──────────┬─────┬─────┬─────────┐
│ chr ┆ pos      ┆ ref ┆ alt ┆ AF_ESP  │
│ --- ┆ ---      ┆ --- ┆ --- ┆ ---     │
│ u8  ┆ u64      ┆ str ┆ str ┆ f64     │
╞═════╪══════════╪═════╪═════╪═════════╡
│ 10  ┆ 16901372 ┆ G   ┆ C   ┆ 0.99992 │
│ 11  ┆ 78121030 ┆ T   ┆ A   ┆ 0.99992 │
└─────┴──────────┴─────┴─────┴─────────┘
```

#### Get sample have variant

Get all variant and sample with GENEINFO equal to 'SAMD11:148398'

```sql
〉select distinct chr, pos, ref, alt, sample from read_parquet('variants.parquet') as v left join read_parquet('genotypes/samples/*') as g on v.id=g.id left join read_parquet('annotations/clinvar.parquet') as a on v.id=a.id WHERE GENEINFO='SAMD11:148398';
┌─────┬────────┬─────┬─────┬────────┐
│ chr ┆ pos    ┆ ref ┆ alt ┆ sample │
│ --- ┆ ---    ┆ --- ┆ --- ┆ ---    │
│ u8  ┆ u64    ┆ str ┆ str ┆ str    │
╞═════╪════════╪═════╪═════╪════════╡
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG003  │
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG001  │
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG007  │
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG004  │
│ …   ┆ …      ┆ …   ┆ …   ┆ …      │
│ 1   ┆ 942934 ┆ G   ┆ C   ┆ HG003  │
│ 1   ┆ 943937 ┆ C   ┆ T   ┆ HG004  │
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG002  │
│ 1   ┆ 942451 ┆ T   ┆ C   ┆ HG006  │
└─────┴────────┴─────┴─────┴────────┘
```

### duckdb


#### Count variants

```sql
D select count(*) from read_parquet('variants.parquet');
┌──────────────┐
│ count_star() │
│    int64     │
├──────────────┤
│      7852699 │
└──────────────┘
```


/// details | check result with `pqrs`
We can check we have same result with pqrs

```bash
➜ pqrs rowcount variants.parquet
File Name: variants.parquet: 7852699 rows
```
///

#### Filter variants from annotations:

Get all variant with a AF_ESP upper than 0.9999

```sql
D select chr, pos, ref, alt, AF_ESP from read_parquet('variants.parquet') as v left join read_parquet('annotations/clinvar.parquet') as c on c.id=v.id where AF_ESP>0.9999;
┌───────┬──────────┬─────────┬─────────┬─────────┐
│  chr  │   pos    │   ref   │   alt   │ AF_ESP  │
│ uint8 │  uint64  │ varchar │ varchar │ double  │
├───────┼──────────┼─────────┼─────────┼─────────┤
│    10 │ 16901372 │ G       │ C       │ 0.99992 │
│    11 │ 78121030 │ T       │ A       │ 0.99992 │
└───────┴──────────┴─────────┴─────────┴─────────┘
```

#### Get sample have variant

Get all variant and sample with GENEINFO equal to 'SAMD11:148398'

```sql
D select distinct chr, pos, ref, alt, sample from read_parquet('variants.parquet') as v left join read_parquet('genotypes/samples/*') as g on v.id=g.id left join read_parquet('annotations/clinvar.parquet') as a on v.id=a.id WHERE GENEINFO='SAMD11:148398';
┌───────┬────────┬─────────┬─────────┬─────────┐
│  chr  │  pos   │   ref   │   alt   │ sample  │
│ uint8 │ uint64 │ varchar │ varchar │ varchar │
├───────┼────────┼─────────┼─────────┼─────────┤
│     1 │ 942451 │ T       │ C       │ HG002   │
│     1 │ 942934 │ G       │ C       │ HG002   │
│     1 │ 942451 │ T       │ C       │ HG003   │
│     1 │ 942934 │ G       │ C       │ HG003   │
│     1 │ 942451 │ T       │ C       │ HG007   │
│     1 │ 943937 │ C       │ T       │ HG007   │
│     1 │ 942451 │ T       │ C       │ HG001   │
│     1 │ 942451 │ T       │ C       │ HG004   │
│     1 │ 943937 │ C       │ T       │ HG004   │
│     1 │ 942451 │ T       │ C       │ HG006   │
├───────┴────────┴─────────┴─────────┴─────────┤
│ 10 rows                            5 columns │
└──────────────────────────────────────────────┘
```
