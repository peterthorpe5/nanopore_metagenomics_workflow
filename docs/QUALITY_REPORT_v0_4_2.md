# Quality report: v0.4.2

Release date: 13 August 2026

## Corrected failures

- Optional PCR evaluation now formats integer and floating-point evidence
  safely, including the zero produced for PCR-negative rows with no expected
  species.
- Non-finite or invalid evidence counts are normalised defensively to zero.
- A generic selective retry can submit only named classifier arrays plus fresh
  aggregation.
- Selective retry accepts an intentional method-resource change while retaining
  stable run, configuration-path and repository identity checks.
- An interrupted selective submission resumes its stored small plan rather
  than expanding into the complete DAG.
- The prepared Dundee benchmark raises Kraken2 from 96,000 MB to 409,600 MB
  after all 11 original tasks ended `OUT_OF_MEMORY`.

## Reusability boundary

PCR truth remains a separate optional reporting input. The reusable workflow
does not derive classifier calls, database contents, thresholds or minimap2
reference requirements from PCR. PCR-free future datasets are covered by the
genuine Snakemake dry-run and report tests.

The 409,600 MB Kraken2 allocation belongs only to
`dundee_real_reads_v2_20260810.yaml`; the generic template retains its
configurable 96,000 MB starting value.

## Automated release gate

- 140 unit and integration tests passed under Python 3.12.
- A genuine Snakemake 9.25.1 dry run passed for a PCR-free future dataset.
- Branch-aware coverage was 91%, above the configured 90% threshold.
- Ruff lint and formatting checks passed for all 24 Python files.
- Python byte-code compilation passed.
- Bash syntax validation passed for every shell script under `scripts/` and
  `workflow/`.
- The source distribution and wheel built successfully as version 0.4.2.
- `git diff --check` passed.

## Selective-retry tests

The test suite verifies that:

- missing classification-ready FASTQs stop a retry before submission;
- unknown, duplicated or disabled retry methods are rejected;
- a minimap2-only retry requires its retained validated index;
- a Kraken2-only retry contains exactly the Kraken2 array and aggregation;
- the updated 409,600 MB value reaches the generated `sbatch` command;
- the earlier complete journal is archived only after active-job checks;
- Metabuli, minimap2 and KmerSutra are absent from the Kraken2 retry journal;
  and
- an interrupted retry resumes only the missing job.
