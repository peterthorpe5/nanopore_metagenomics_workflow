# Quality report: v0.1.0

Date: 4 August 2026

Scope: standalone Snakemake workflow for host depletion and Kraken2, Metabuli
and optional KmerSutra classification of real Oxford Nanopore FASTQ data.

## Automated gate

- 69 unit and integration tests passed.
- Branch-aware coverage: 96% overall; required minimum: 90%.
- Ruff lint: passed.
- Python compilation: passed.
- Shell syntax validation: passed.
- Genuine Snakemake DAG dry-run with synthetic FASTQ and resources: passed.
- Python wheel and source distribution build: passed.
- Packaged Snakefile presence: verified.

## Safety and recovery coverage

Tests cover invalid manifests and configuration, `/tmp` capacity and filesystem
checks, resource staging validation, binary compressed output, atomic result
publication, completion tokens, interrupted tasks, per-sample KmerSutra
timeouts and failures, compact failure evidence, and restart behaviour.

The tests use synthetic data and fake executable adapters. Scientific
acceptance testing against the production host reference and classifier
databases remains required when the first real samples arrive.
