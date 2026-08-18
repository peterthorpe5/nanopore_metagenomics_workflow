# Quality report: v0.4.3

## Scope

Version 0.4.3 generalises the independent truth comparison beyond Plasmodium,
retains only exact species-level rows for formal concordance, and adds the
ATCC MSA-1003 HiFi reference-composition recovery preset.

## Validation

- 142 unit/integration tests passed; one Snakemake dry-run test was skipped
  because Snakemake was unavailable in the validation runtime.
- Python byte-compilation completed for `src/` and `tests/`.
- The recovery shell passed `bash -n` syntax validation.
- The ATCC preset test confirms that the truth table and minimap2 requirement
  contain the same 20 current taxonomy names.

## Behavioural safeguards

- Kraken2 rank `S` and Metabuli rank `species` are the only native classifier
  rows admitted to exact species concordance.
- Strain, species-group, serogroup and descendant rows remain available in the
  descriptive evidence table but cannot inflate the additional-species count.
- Failed or unavailable classifiers remain unavailable rather than becoming
  biological non-detections.
- The recovery does not rerun the successful Kraken2 or Metabuli branches.
