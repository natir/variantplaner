# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->
## [0.3.0](https://github.com/natir/variantplaner/releases/tag/0.3.0) - 2024-09-30

<small>[Compare with 0.2.4](https://github.com/natir/variantplaner/compare/0.2.4...0.3.0)</small>

### Features

- move ped code from io to object ([b6638ca](https://github.com/natir/variantplaner/commit/b6638ca7a562316a25c18e54459a74926582357f) by Pierre Marijon).
- parquet2vcf can extract only variant of one chromosome ([f29c6ac](https://github.com/natir/variantplaner/commit/f29c6ac2ea048714729515422858087195b31417) by Pierre Marijon).
- improve alt '*' management ([755a5fa](https://github.com/natir/variantplaner/commit/755a5fad7b822463495c3bae391542872d0a0ac2) by Pierre Marijon).
- partition support append mode ([1b6603f](https://github.com/natir/variantplaner/commit/1b6603ffd24b2f2d81bdf8198ae4378ef9e6a78e) by Pierre Marijon).
- variant merge support append ([d12d357](https://github.com/natir/variantplaner/commit/d12d3574a91e0ff5dc33d7171097a9a5f4faf8cd) by Pierre Marijon).
- Add python 3.12 in ci test and support ([c51c932](https://github.com/natir/variantplaner/commit/c51c932180b1d7bba2e00c5720cc7bd4bf987e33) by Pierre Marijon).
- Better variant hash we only need store ref length ([8a3db5e](https://github.com/natir/variantplaner/commit/8a3db5eab8c91b727a1e91f82b15c52c963c1a8c) by Pierre Marijon).
- add documentation on how interogate genotype variants partitions ([b843bca](https://github.com/natir/variantplaner/commit/b843bca05aa463703e2f708450b48656512dbcaa) by Pierre Marijon).
- New cli ([c2cb033](https://github.com/natir/variantplaner/commit/c2cb033ce00f548f7d6f1041d5381adfc0debcf5) by Pierre Marijon).

### Bug Fixes

- add output file in struct operation only if it's exist ([9a868ad](https://github.com/natir/variantplaner/commit/9a868ad58d12210fc46ccf35a0503e6bceae28ba) by Pierre Marijon).
- not failled if SVTYPE or SVLEN column isn't present ([32a8dad](https://github.com/natir/variantplaner/commit/32a8dad0d5f5654d06d39e79cc03b53cb2c91e13) by Pierre Marijon).
- correct benchmark script run ([8ebdca0](https://github.com/natir/variantplaner/commit/8ebdca08f1f7e7871771050f7df43aa96be825d6) by Pierre Marijon).

### Code Refactoring

- Move many variant operation in object. ([df36fa0](https://github.com/natir/variantplaner/commit/df36fa0db9e933030b4b72f92645f5717b74597e) by Pierre Marijon).

## [0.2.4](https://github.com/natir/variantplaner/releases/tag/0.2.4) - 2023-12-21

<small>[Compare with 0.2.4](https://github.com/natir/variantplaner/compare/0.2.3...0.2.4)</small>

### Features

- add variantplaner logo.
- add method to compute partition value of id.

### Bug Fixes

- fix: #41 vcf spec indicate Integer must be store in 32bits.
- update to polars 0.20.

## [0.2.3](https://github.com/natir/variantplaner/releases/tag/0.2.2) - 2023-11-21

<small>[Compare with 0.2.2](https://github.com/natir/variantplaner/compare/0.2.2..0.2.3)</small>

### Features:

- Use ruff format in place of black
- Minimal polars version is 0.19.15
- SV cannot have small variant id
- Type of id is include in part of id (all long variant are store in same place)

### Fix:

- Usage intgerate chromosome2length option

## [0.2.2](https://github.com/natir/variantplaner/releases/tag/0.2.2) - 2023-10-19

<small>[Compare with first 0.2.0](https://github.com/natir/variantplaner/compare/0.2.1...0.2.2)</small>

### Fix:

- Correct typo in readme
- Correct change in variantplaner_rs/Cargo.lock

## [0.2.1](https://github.com/natir/variantplaner/releases/tag/0.2.1) - 2023-10-09

<small>[Compare with first 0.2.0](https://github.com/natir/variantplaner/compare/0.2.0...0.2.1)</small>

### Fix:

- Hot fix some lazy trouble

## [0.2.0](https://github.com/natir/variantplaner/releases/tag/0.2.0) - 2023-10-02

<small>[Compare with first 0.1.0](https://github.com/natir/variantplaner/compare/0.1.0...0.2.0)</small>

### Features:

- Replace hash variant id by a uniq variant id
- Improve variant transmission system to support genome with ploidie lower than 92
- Genotype partitioning now use position part of unique id to get more local request

### Fixes:

- Documentation fix
- Test coverage improvement

## [0.1.0](https://github.com/natir/variantplaner/releases/tag/0.1.0) - 2023-07-25

<small>[Compare with first commit](https://github.com/natir/variantplaner/compare/265a95ea26746b7aa796c3df6cee2451a608dd49...0.1.0)</small>
