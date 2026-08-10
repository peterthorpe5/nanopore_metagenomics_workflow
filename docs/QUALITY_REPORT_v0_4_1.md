# Quality report: v0.4.1

Release date: 10 August 2026

## Scope

This release provides a reusable, detached and journalled Slurm DAG containing
independent per-sample arrays. It supports optional exact PCR truth validation
for benchmark datasets, restart-safe per-sample completion records, and
minimap2 protections against multipart-index count inflation. Routine real
samples can run without a PCR truth table.

The detached Slurm backend accepts host-removed manifests. Raw-read host
depletion remains available through the Snakemake backend and is rejected
early by detached planning rather than submitted as one bundled job.

The prepared Dundee run uses the reduced focused reference at:

```text
/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta
```

It does not build or use the much larger shared bacterial, viral and
Plasmodium reference.

## Automated release gate

- 134 unit and integration tests pass under Python 3.12.
- The suite includes a genuine Snakemake 9 dry run of the packaged DAG.
- Branch-aware coverage is 91%, above the configured 90% threshold.
- Ruff lint and formatting verification pass for all 24 Python files.
- Python byte-code compilation passes.
- Bash syntax validation passes for every shell script beneath `scripts/` and
  `workflow/`.
- Both the source distribution and Python wheel build successfully.

## Minimap2 regression checks

The suite verifies that:

- prepared configuration selects the exact reduced masked FASTA and leaves
  `genome_config` and `index` empty;
- legacy header forms including `Plasmodium inui`, `Plas_inui` and `P.inui`
  are resolved to the same canonical species;
- the masked FASTA must expose every species explicitly configured under
  `minimap2.required_species` before index building;
- the run creates a checksum-bound index and accepts exactly one minimap2 index
  part;
- the reference and index have hard size limits;
- PAF targets must occur in the configured FASTA;
- repeated non-consecutive query blocks are rejected as a multipart signature;
  and
- mapped query groups cannot exceed the validated number of input reads.

## Workflow and failure checks

The suite also exercises:

- eight-job prepared planning: preflight, host-removed input acceptance, index
  construction, four independent classifier arrays and aggregation;
- KmerSutra independence from Kraken2, Metabuli and minimap2;
- `4week` QoS on KmerSutra only;
- durable journalling after every accepted `sbatch` response;
- canonical Mac-to-GitHub-to-GPFS checkout identity validation during both
  safe planning and real submission;
- interrupted submission resume and deliberate new-attempt safeguards;
- per-sample timeout, scheduler-failure, missing-resource and invalid-output
  states;
- partial reports after classifier failure;
- atomic publication and restart from checksum-bound completion records; and
- output retention rules that avoid uncompressed large per-read tables.

## PCR and reporting checks

- PCR truth is optional for ordinary real-sample runs and remains separate from
  the sample manifest.
- PCR-free runs retain classifier reports, mark PCR as `not_configured`, and do
  not make an accuracy claim.
- When PCR is configured, every positive comparison species must also be
  declared explicitly under `minimap2.required_species`.
- The release manifest has 11 unique samples and a matching 11-row truth table.
- Ten samples enter the primary PCR comparison.
- `MRC1123_WuH001WB` is retained for classification but remains explicitly
  excluded because its PCR interpretation was not supplied.
- All ten included samples expect `Plasmodium inui`; six also expect
  `Plasmodium cynomolgi`.
- A missing, failed, timed-out or malformed classifier result remains
  unavailable and is never converted into a biological negative.
- Reports retain method-native tables, normalised evidence, exact denominators,
  offline HTML and provenance checksums.
