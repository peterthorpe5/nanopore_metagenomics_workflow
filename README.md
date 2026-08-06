# Nanopore real-data classification workflow

This standalone Snakemake package classifies real Oxford Nanopore FASTQ or
FASTQ.GZ samples. It can run as a production workflow and provide the
real-sample classifier evidence used in the KmerSutra study, without coupling
the workflow implementation to either paper repository.

The workflow supports two declared input states:

- `raw`: remove reads mapping to a configured host before classification;
- `host_removed`: validate and stage reads that were depleted previously,
  without performing host mapping again.

For each logical sample it:

1. validates the sample sheet, references, databases and software;
2. prepares raw or already host-removed reads according to `inputs.read_state`;
3. classifies the prepared reads independently with Kraken2 and Metabuli;
4. maps them to a controlled reference with minimap2;
5. optionally screens them with KmerSutra; and
6. writes separate TSV summaries, provenance and checksums.

KmerSutra is optional and failure-tolerant by default. A KmerSutra error or
time limit is recorded per sample, but it does not invalidate Kraken2,
Metabuli or minimap2 results and does not prevent final aggregation.

## Execution and storage model

The Dundee template uses `/tmp`, which has been verified on Barton compute
nodes as a writable, job-specific XFS filesystem with approximately 7 TB free.
Each heavy rule:

- creates a bounded workspace beneath `/tmp`;
- stages its reference or classifier database once per job;
- stages the required FASTQ inputs;
- performs all tool work locally;
- validates outputs;
- atomically `rsync`s declared results and logs to the persistent run root; and
- removes its local workspace on exit.

Kraken2, Metabuli and minimap2 process all samples sequentially within their
respective staged jobs. A global Snakemake resource permits only one large
database/reference staging job at a time. KmerSutra runs as one serial,
restartable job per sample because it is much slower and more likely to fail.

Prepared FASTQ files are Snakemake temporary outputs by default and disappear
after all four classifier branches finish. Kraken2 and Metabuli per-read
assignments and minimap2 PAF evidence are retained only as gzip-compressed
declared outputs. This is configurable.

## Inputs

Copy the templates:

```bash
cp config/real_data.template.yaml config/real_data.yaml
cp config/samples.template.tsv config/samples.tsv
```

The sample sheet is tab-separated. Required columns are `sample_id` and
`fastq`. Optional columns are `run_id`, `barcode` and `description`. A sample
may occupy several rows when Nanopore output is split into several FASTQ
chunks. The row order defines the input-part order.

Sample IDs may contain letters, numbers, dots, underscores and hyphens. Every
FASTQ path must be unique and may end in `.fastq`, `.fastq.gz`, `.fq` or
`.fq.gz`.

Set `inputs.read_state` deliberately:

```yaml
inputs:
  samples: /absolute/path/to/samples.tsv
  read_state: host_removed
```

When the value is `host_removed`, `host.reference` and `host.index` may be
blank. The workflow validates, merges if necessary and re-compresses the
declared inputs in `/tmp`; it records `host_depletion_performed=false` and does
not invoke minimap2 for host removal.

The controlled classification reference under `minimap2.reference` is
independent of the optional host reference. If its matching MMI cannot be
proven, leave `minimap2.index` blank. The workflow builds one run-local index
from the declared FASTA and reuses it.

## Installation

From the repository root:

```bash
conda env create --file environment/environment.yml
conda activate nanopore_realdata_workflow
python -m pip install --editable '.[test]'
```

If the environment already exists, update it after pulling a new release:

```bash
conda env update \
    --name nanopore_realdata_workflow \
    --file environment/environment.yml \
    --prune
conda activate nanopore_realdata_workflow
python -m pip install --editable '.[test]'
```

KmerSutra is not bundled into this independent package. Install the current
KmerSutra checkout into the same environment when that branch is enabled:

```bash
KMERSUTRA_REPO=/absolute/path/to/kmersutra
python -m pip install --editable "${KMERSUTRA_REPO}"
kmersutra-screen --help >/dev/null
```

To run without KmerSutra, set `kmersutra.enabled: false` and optionally leave
`databases.kmersutra_panel` blank. The other three classifiers still run.

## Validate and dry-run

```bash
nanopore-realdata \
    --action validate \
    --config config/real_data.yaml

bash scripts/run_workflow.sh \
    --config config/real_data.yaml \
    --dry-run \
    --verbose
```

## Prepared Dundee real-read benchmark

Release v0.2.0 includes a sample sheet and configuration for the 11 MRC FASTQ
files under:

```text
/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/real_reads
```

The prepared configuration declares those reads as `host_removed` and uses:

- the Plasmodium benchmarking Kraken2 database;
- the matching Metabuli database;
- `shared_bact_viral_plasmo_refs.cleaned.fa` for controlled minimap2 mapping;
- the KmerSutra v0.46 raw-ONT LOD-balanced panel supplied for this run.

The existing `shared_bact_viral_plasmo_refs.mmi` is not selected because its
filename associates it with the uncleaned FASTA and no checksum sidecar proves
that it was built from the cleaned reference. The workflow safely builds
`00_preflight/classification_reference.mmi` from the cleaned FASTA instead.

Validate and inspect the exact DAG:

```bash
conda activate nanopore_realdata_workflow

nanopore-realdata \
    --action validate \
    --config config/dundee_real_reads_20260806.yaml

bash scripts/run_dundee_real_reads_20260806.sh --dry-run
```

Start or resume the full run:

```bash
bash scripts/run_dundee_real_reads_20260806.sh
```

If Snakemake reports a stale lock after an interrupted controller process:

```bash
bash scripts/run_dundee_real_reads_20260806.sh --unlock
bash scripts/run_dundee_real_reads_20260806.sh
```

The Dundee profile uses the `barton` account and `barton` partition. Snakemake
provides restart behaviour: completed results are skipped and incomplete
stages are rerun.

KmerSutra starts only after Kraken2, Metabuli and minimap2 finish. Each sample
receives a normal three-day allocation and a 70-hour internal limit. The
internal limit leaves time to publish a compact failure record before Slurm
ends the job. Only if real evidence shows that longer execution is useful
should the available `4week` QoS be requested:

```bash
bash scripts/run_workflow.sh \
    --config config/dundee_real_reads_20260806.yaml \
    --profile workflow/profiles/dundee \
    --jobs 10 \
    --set-resource classify_kmersutra:runtime=40320 \
    --set-resource classify_kmersutra:slurm_qos=4week \
    --verbose
```

The KmerSutra `timeout_minutes` setting must also be increased below the chosen
Slurm runtime if this exceptional mode is used.

## Output layout

```text
00_preflight/
  preflight.json
  resolved_samples.tsv
  software_versions.tsv
  classification_reference.mmi
01_host_depletion/
  <sample_id>/
    non_host.fastq.gz
    host_removal_summary.tsv
    metadata.json
    complete.json
02_classification/
  kraken2/<sample_id>/
  metabuli/<sample_id>/
  minimap2/<sample_id>/
  kmersutra/<sample_id>/
03_final/
  sample_summary.tsv
  classifier_taxon_reports.tsv.gz
  minimap2_taxon_reports.tsv.gz
  kmersutra_species_calls.tsv.gz
  SHA256SUMS.tsv
  workflow.complete.json
```

For `host_removed` input, `01_host_depletion` is a compatibility stage name:
its metadata and summary explicitly state that depletion was not performed.
The prepared FASTQ disappears after successful downstream use when
`keep_non_host_fastq` is false. A failed KmerSutra sample retains only compact
failure evidence; its final table records the workflow status.

Controlled minimap2 mapping uses `map-ont`, retains secondary PAF alignments
for ambiguity assessment, and reports best alignments that meet MAPQ 15 and
alignment-block length 500 by default. Cross-taxon best-hit ties remain
ambiguous rather than being forced into a unique call.

## Quality checks

```bash
bash scripts/run_quality_checks.sh \
    --results-dir /absolute/path/to/test_results
```

The gate runs Ruff lint and formatting checks, unit and integration tests,
branch-aware coverage with a 90% minimum, Python compilation, shell syntax,
and a genuine Snakemake DAG dry-run using synthetic inputs and fake resources.

## Interpretation

These methods produce different evidence. Kraken2 and Metabuli reports are
harmonised into one taxon-report table. Controlled minimap2 results remain a
separate alignment-derived table, and KmerSutra calls remain a separate
species-evidence table. The workflow does not force them into a false
consensus. Any paper comparison must retain database and software versions and
apply pre-declared analysis thresholds consistently across samples.
