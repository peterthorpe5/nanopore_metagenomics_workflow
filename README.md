# Nanopore real-data classification workflow

This standalone Snakemake package processes real Oxford Nanopore FASTQ or
FASTQ.GZ samples. It is deliberately separate from the Plasmodium benchmarking
and KmerSutra paper workflows.

For each logical sample it:

1. validates the sample sheet, references, databases and software;
2. removes reads mapping to a configurable host reference with minimap2 and
   samtools;
3. classifies the non-host reads independently with Kraken2 and Metabuli;
4. optionally screens the same reads with KmerSutra; and
5. writes harmonised TSV summaries, compact provenance and checksums.

KmerSutra is optional and failure-tolerant by default. A KmerSutra error or
time limit is recorded per sample, but it does not invalidate successful
Kraken2 or Metabuli analysis and does not prevent final aggregation.

## Execution and storage model

The Dundee template uses `/tmp`, which has been verified on the Barton compute
nodes as a writable, job-specific XFS filesystem with approximately 7 TB free.
Each heavy rule:

- creates a bounded workspace beneath `/tmp`;
- stages its host index or classifier database once per job;
- stages the required FASTQ inputs;
- performs all tool work locally;
- validates outputs;
- atomically `rsync`s declared results and logs to the persistent run root; and
- removes its local workspace on exit.

Kraken2 and Metabuli process all samples sequentially within their respective
staged jobs. This avoids copying a roughly 55 GiB Kraken2 database once per
sample. A global Snakemake resource permits only one large database-staging job
at a time.

Non-host FASTQ files are Snakemake temporary outputs by default and are removed
after all three classifier branches finish. Kraken2 and Metabuli per-read
assignments are retained only as gzip-compressed declared evidence. This is
configurable.

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

Edit `config/real_data.yaml` to supply:

- the output directory;
- the sample-sheet path;
- a host FASTA and optionally an existing minimap2 MMI index;
- the Kraken2 database directory;
- the Metabuli database directory; and
- the KmerSutra panel path.

The host FASTA can represent any appropriate host. If several host genomes
must be removed, provide one combined reference FASTA with unique sequence
identifiers.

## Installation

From the repository root:

```bash
conda env create --file environment/environment.yml
conda activate nanopore_realdata_workflow
python -m pip install --editable '.[test]'
```

KmerSutra is not bundled into this independent package. Install the current
KmerSutra checkout into the same environment only when that branch is enabled:

```bash
python -m pip install --editable /absolute/path/to/kmersutra
kmersutra-screen --help >/dev/null
```

To run without KmerSutra, set `kmersutra.enabled: false` and optionally leave
`databases.kmersutra_panel` blank. Kraken2 and Metabuli still run normally.

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

## Run on Dundee Slurm

The normal profile uses the `barton` account and `barton` partition. It does
not request the `4week` QoS.

```bash
conda activate nanopore_realdata_workflow

bash scripts/run_workflow.sh \
    --config config/real_data.yaml \
    --profile workflow/profiles/dundee \
    --jobs 10 \
    --verbose
```

Snakemake provides restart behaviour. Rerun the same command after a failed or
interrupted execution. Valid completed results are skipped, incomplete outputs
are rerun and databases are staged again only for the stage that needs them.

KmerSutra runs as one serial, restartable job per sample, after Kraken2 and
Metabuli. Each sample currently receives a normal three-day allocation and a
70-hour internal time limit. The internal limit allows it to write a failure
record before Slurm ends the job. Only if real evidence shows that longer
execution is scientifically useful should the available `4week` QoS be
requested:

```bash
bash scripts/run_workflow.sh \
    --config config/real_data.yaml \
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
  host_reference.mmi
01_host_depletion/
  <sample_id>/
    non_host.fastq.gz
    host_removal_summary.tsv
    host_depletion.log
    metadata.json
    complete.json
02_classification/
  kraken2/<sample_id>/
  metabuli/<sample_id>/
  kmersutra/<sample_id>/
03_final/
  sample_summary.tsv
  classifier_taxon_reports.tsv.gz
  kmersutra_species_calls.tsv.gz
  SHA256SUMS.tsv
  workflow.complete.json
```

`non_host.fastq.gz` disappears after successful downstream use when
`keep_non_host_fastq` is false. A KmerSutra sample that fails has a compact
`failure.json` plus a retained failed-attempt log; the final KmerSutra table
records its workflow status.

## Quality checks

```bash
bash scripts/run_quality_checks.sh \
    --results-dir /absolute/path/to/test_results
```

The gate runs Ruff, the unit and integration tests, branch-aware coverage with
a 90% minimum, Python compilation, shell syntax validation and whitespace
checks. Tests use synthetic FASTQ and fake executables; they do not require the
real classifier databases.

## Important interpretation

The three tools produce different kinds of evidence. Kraken2 and Metabuli
reports are harmonised into a shared taxon-report table, while KmerSutra calls
remain a separate species-evidence table. The workflow does not force these
outputs into a false consensus call. Scientific interpretation and any
pathogen-specific reporting thresholds should be defined after the real sample
design and expanded database panel are known.
