# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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

<!-- insertion marker -->
## [0.1.0](https://github.com/natir/variantplaner/releases/tag/0.1.0) - 2023-07-25

<small>[Compare with first commit](https://github.com/natir/variantplaner/compare/265a95ea26746b7aa796c3df6cee2451a608dd49...0.1.0)</small>
