# Quality report: v0.3.0

Release date: 2026-08-06

## Scope

This release adds classifier-independent failure handling, offline HTML reports,
high-memory minimap2 settings, the Dundee four-week KmerSutra allocation, and a
safe new-dataset launcher and runbook.

## Automated validation

- 103 unit and integration tests pass under Python 3.12.
- The suite includes a genuine Snakemake 9 DAG dry run.
- The last branch-aware coverage run measured 95%, above the configured 90%
  release threshold.
- Ruff lint passes.
- Ruff formatting verification passes for all 22 Python files.
- Python byte-code compilation passes.
- Bash syntax validation passes for every script beneath `scripts/`.
- A Python 3 wheel builds successfully with no unresolved package dependency.

## Failure and reporting checks

The test suite exercises:

- independent `continue` policies for Kraken2, Metabuli, minimap2 and KmerSutra;
- partial success within bundled multi-sample classifier jobs;
- missing executable and missing database/reference handling;
- KmerSutra and minimap2 pipeline timeouts;
- all-classifier failure with a truthful final failed-status report;
- recovery reporting after a non-zero workflow-controller exit;
- malformed classifier tables becoming warnings rather than report loss;
- atomic output publication and `/tmp` cleanup;
- HTML escaping, balanced report sections and absence of external assets;
- final, comparison, per-classifier and per-sample HTML creation; and
- safe new-dataset initialisation with overwrite refusal.

## Dundee resource assertions

- The supplied minimap2 index and classification rules request 160,000 MB RAM.
- The supplied KmerSutra launcher overrides only `classify_kmersutra` to runtime
  40,320 minutes and `slurm_qos=4week`.
- KmerSutra's internal timeout is 40,100 minutes, leaving 220 minutes for status
  publication and scratch cleanup.

## Scientific safeguards

- Native evidence is retained for every method.
- Cross-classifier overlap is labelled descriptive rather than an accuracy
  estimate.
- No classifier is silently promoted to a truth set.
- Missing evidence is visible by sample and method.
- HTML pages are self-contained, portable and usable without internet access.

