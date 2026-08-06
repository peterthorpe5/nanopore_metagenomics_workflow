# Quality report: v0.2.0

Date: 6 August 2026

Scope: standalone Snakemake workflow for raw or already host-removed Oxford
Nanopore FASTQ data, using Kraken2, Metabuli, controlled-reference minimap2 and
optional failure-tolerant KmerSutra classification.

## Release-specific behaviour

- `inputs.read_state: host_removed` bypasses host mapping and records that host
  depletion was not performed.
- The 11 supplied MRC FASTQ.GZ files have a ready Dundee sample sheet.
- The ready configuration uses the supplied Kraken2 and Metabuli databases and
  the exact KmerSutra v0.46 raw-ONT LOD-balanced panel.
- A run-local minimap2 MMI is built from
  `shared_bact_viral_plasmo_refs.cleaned.fa`; the existing ambiguously named MMI
  is not assumed to match.
- Minimap2 retains controlled-reference PAF evidence and reports unique and
  ambiguous best taxa at MAPQ 15 and 500 bp alignment length by default.
- Every heavy stage runs beneath job-local `/tmp`, publishes declared results
  atomically and removes its workspace on exit.

## Automated gate

- 81 unit and integration tests passed.
- Branch-aware coverage: 95% overall; required minimum: 90%.
- Ruff lint and formatting checks: passed.
- Python compilation: passed.
- Shell syntax validation: passed.
- Genuine Snakemake 9 DAG dry-run with synthetic host-removed FASTQ and all
  four classifier branches: passed.
- Python wheel and source distribution build: passed.
- Independent release extraction and complete retest: passed.

## Safety and recovery coverage

Tests cover invalid manifests and configuration, read-state routing, malformed
FASTQ, `/tmp` capacity and filesystem checks, reference and database staging,
FASTA/PAF parsing, ambiguous minimap2 best hits, atomic publication, completion
tokens, interrupted tasks, KmerSutra timeouts and failures, compact failure
evidence and restart behaviour.

Tests use synthetic data and fake classifier executables. The prepared Dundee
configuration still requires a preflight validation against the actual paths,
installed tool versions and production database inventories before execution.
