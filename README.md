# Nanopore metagenomics workflow

Version 0.3.0 is a standalone, production-oriented Snakemake workflow for real
Oxford Nanopore FASTQ or FASTQ.GZ data. It prepares raw or already host-removed
reads, runs Kraken2, Metabuli, controlled-reference minimap2 and KmerSutra as
independent branches, and produces complete tabular plus offline HTML reports.

The workflow is suitable for the real-read classifier comparison supporting
the KmerSutra study, but it is not coupled to either manuscript repository.
It does not make a clinical diagnosis and does not force classifier outputs
into a single consensus call.

## What changed in v0.3.0

- Every classifier has an independent `continue` or `fail` policy.
- Kraken2, Metabuli, minimap2 and KmerSutra failures no longer have to block
  final reporting.
- Missing tools, missing resources, command errors, timeouts and invalid output
  are represented by explicit status records.
- Successful samples are retained when another sample fails inside a bundled
  classifier job.
- Final aggregation depends on compact terminal status records, not on every
  classifier report existing.
- A self-contained HTML report is written for each classifier.
- A cross-classifier comparison report exposes agreement and disagreement
  without treating one method as truth.
- A final run dashboard and one report per sample are generated.
- All HTML is offline: CSS, JavaScript and data are embedded, with no CDN or
  internet dependency.
- Tables can be filtered, sorted and exported as visible TSV from the browser.
- The final directory includes reusable normalised evidence and JSON report
  data alongside the native method tables.
- The prepared Dundee run gives minimap2 160 GB RAM.
- The prepared KmerSutra run uses the Barton `4week` QoS explicitly.
- `environment.yml` no longer tries to install `-e .[test]` from the wrong
  `environment/` directory.
- `scripts/start_new_dataset.sh` provides a safe starting point for a new run.

## Workflow outline

For each logical sample, the workflow:

1. validates the core configuration, sample sheet and FASTQ records;
2. records readiness separately for every classifier;
3. removes host reads when `inputs.read_state: raw`, or accepts previously
   depleted reads when `inputs.read_state: host_removed`;
4. stages each large database/reference and required FASTQ data into job-local
   `/tmp`;
5. runs Kraken2, Metabuli and controlled-reference minimap2 independently;
6. runs KmerSutra as one restartable job per sample;
7. validates and atomically publishes completed scientific outputs;
8. records failures and timeouts as compact JSON rather than publishing partial
   classifier tables; and
9. aggregates every available validated result into TSV, JSON and HTML reports.

Kraken2, Metabuli and minimap2 each process all samples sequentially inside one
resource-staged job. This avoids repeatedly copying a large database for every
sample. KmerSutra remains one serial job per sample because it is much slower
and more likely to require isolated restart or failure handling.

## Failure model

The default production configuration uses `failure_policy: continue` for all
four classifiers.

| Event | Other classifiers run? | Final HTML generated? | Recorded status |
|---|---:|---:|---|
| Classifier succeeds | Yes | Yes | `success` |
| Some samples succeed in a bundled classifier | Yes | Yes | stage `partial` |
| Command exits non-zero | Yes | Yes | `failed` |
| Internal time limit is reached | Yes | Yes | `timeout` |
| Executable or configured resource is absent | Yes | Yes | `unavailable` |
| KmerSutra is disabled | Yes | Yes | `disabled` |
| Terminal metadata cannot be read | Yes | Yes | `invalid` |
| No terminal record exists | Yes | Yes | `missing` |

The host-preparation stage is intentionally different. Raw reads cannot be
classified correctly if host depletion fails, so host preparation remains a
core workflow dependency. For `host_removed` input, the existing depleted
reads are validated and prepared without performing host mapping again.

If a classifier is configured with `failure_policy: fail`, its exception is
allowed to stop Snakemake. This stricter mode is available for projects that
require every configured method to finish.

Scheduler-level hard kills cannot execute cleanup code. Each classifier
therefore has an internal time limit below its requested Slurm runtime. This
leaves time to publish a compact status record and remove the `/tmp` workspace.
For an out-of-memory kill, Slurm may terminate the process immediately; rerun
with more memory if a status record cannot be written. Already completed
branches remain restartable and are not recalculated.

If the Snakemake controller nevertheless exits non-zero during a real run,
`scripts/run_workflow.sh` automatically attempts the equivalent of:

```bash
nanopore-realdata --action report --config config/my_run.yaml
```

This salvage action labels absent terminal records as `missing` and builds the
final reports from every result that did complete. The launcher then returns
the original non-zero exit code: reporting is rescued, but the workflow failure
is not hidden. Dry runs and unlock operations never trigger salvage reporting.

## Execution and `/tmp` contract

The Dundee configuration uses job-local `/tmp`. Each heavy rule:

- creates a unique workspace beneath `/tmp`;
- verifies scratch is local and has sufficient capacity;
- stages its database, panel or reference once per job;
- stages only the FASTQ inputs needed by that job;
- writes bulky working data inside `/tmp`;
- validates declared outputs before publication;
- atomically `rsync`s completed outputs to persistent storage; and
- removes the temporary workspace on success, failure or caught timeout.

Failed classifier attempts retain compact logs and `failure.json`; incomplete
large outputs are not published as valid results.

## Installation

From the repository root on Dundee:

```bash
cd ~/data/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow

conda env create --file environment/environment.yml
conda activate nanopore_realdata_workflow

python -m pip install --editable '.[test]'
```

The editable package installation must be run from the repository root, where
`pyproject.toml` exists. It is deliberately not embedded in
`environment/environment.yml` because pip resolves relative paths from the
environment file's directory.

To update an existing environment:

```bash
cd ~/data/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow

conda env update \
    --name nanopore_realdata_workflow \
    --file environment/environment.yml \
    --prune

conda activate nanopore_realdata_workflow
python -m pip install --editable '.[test]'
```

KmerSutra is an independent package. Install its current checkout into the same
environment when `kmersutra.enabled: true`:

```bash
KMERSUTRA_REPO="/absolute/path/to/kmersutra"

python -m pip install --editable "${KMERSUTRA_REPO}"
kmersutra-screen --help >/dev/null
```

Verify the complete environment:

```bash
kraken2 --version
metabuli --help >/dev/null
minimap2 --version
snakemake --version
nanopore-realdata --help >/dev/null
kmersutra-screen --help >/dev/null
python -m pytest -q
```

## Prepared Dundee real-read benchmark

The supplied run configuration covers the 11 MRC FASTQ.GZ files beneath:

```text
/home/pthorpe001/data/2026_plasmodium_kraken_sensitivity/real_reads
```

They are declared as already host-removed. The workflow does not perform host
depletion again.

The prepared configuration uses:

- the Plasmodium benchmarking Kraken2 database;
- the matching custom Metabuli database;
- `shared_bact_viral_plasmo_refs.cleaned.fa` for controlled minimap2 mapping;
- a new run-local minimap2 index built from that cleaned FASTA; and
- the KmerSutra v0.46 raw-ONT LOD-balanced `species_kmer_panel.tsv.gz`.

The existing `shared_bact_viral_plasmo_refs.mmi` is not used because its name
associates it with the uncleaned FASTA and no checksum sidecar proves that the
index matches the cleaned reference.

The prepared resources are:

| Stage | Memory | Runtime/QoS |
|---|---:|---|
| Kraken2 | 96 GB | normal Barton allocation |
| Metabuli | 160 GB | normal Barton allocation |
| minimap2 index + classification | 160 GB | normal Barton allocation |
| KmerSutra per sample | 96 GB | 40,320 minutes, `4week` QoS |

KmerSutra's internal limit is 40,100 minutes, leaving 220 minutes for status
publication and cleanup before the four-week scheduler allocation ends.

Validate the paths and inspect the DAG:

```bash
conda activate nanopore_realdata_workflow

nanopore-realdata \
    --action validate \
    --config config/dundee_real_reads_20260806.yaml

bash scripts/run_dundee_real_reads_20260806.sh --dry-run
```

Start or resume the full benchmark inside `tmux` or `screen`:

```bash
bash scripts/run_dundee_real_reads_20260806.sh
```

The launcher applies the `4week` QoS only to `classify_kmersutra`. The other
classifier jobs continue to use the default Barton partition and QoS.

If Snakemake reports a stale controller lock:

```bash
bash scripts/run_dundee_real_reads_20260806.sh --unlock
bash scripts/run_dundee_real_reads_20260806.sh
```

## Starting a new dataset

### 1. Create a protected config/sample-sheet pair

Use the new helper from the repository root:

```bash
bash scripts/start_new_dataset.sh \
    --action initialise \
    --config config/my_new_run.yaml
```

This creates, without overwriting existing files:

```text
config/my_new_run.yaml
config/my_new_run.samples.tsv
```

If either path already exists, the command stops with exit code 3 and leaves
both untouched.

You may choose the sample-sheet path explicitly:

```bash
bash scripts/start_new_dataset.sh \
    --action initialise \
    --config config/my_new_run.yaml \
    --samples config/my_inputs.tsv
```

### 2. Fill in the TSV sample sheet

The file is tab-separated. It must contain `sample_id` and `fastq`; optional
columns are `run_id`, `barcode` and `description`.

```text
sample_id	fastq	run_id	barcode	description
sample_A	/absolute/path/sample_A.fastq.gz	run_1	barcode01	first sample
sample_B	/absolute/path/sample_B.fastq.gz	run_1	barcode02	second sample
```

For several Nanopore chunks belonging to one logical sample, repeat the sample
ID on successive rows. Their row order is preserved during merging.

Every FASTQ path must be unique, non-empty and end with `.fastq`, `.fastq.gz`,
`.fq` or `.fq.gz`. Sample IDs may contain letters, numbers, dots, underscores
and hyphens.

### 3. Edit the YAML deliberately

At minimum, change:

- `run.id` and `run.output_directory`;
- `inputs.samples` to the new TSV;
- `inputs.read_state` to exactly `raw` or `host_removed`;
- host reference/index when using raw inputs;
- Kraken2 and Metabuli database paths;
- minimap2 classification reference and optional matching index;
- KmerSutra panel or `kmersutra.enabled: false`;
- memory, threads and runtime for the dataset size; and
- `reporting.focus_taxa` when the scientific focus is not Plasmodium.

For already depleted reads:

```yaml
inputs:
  samples: /absolute/path/to/my_new_run.samples.tsv
  read_state: host_removed

host:
  reference: ""
  index: ""
```

For raw reads, declare a real host reference. Leave `host.index` blank to build
a persistent run-local index:

```yaml
inputs:
  samples: /absolute/path/to/my_new_run.samples.tsv
  read_state: raw

host:
  reference: /absolute/path/to/host_reference.fasta
  index: ""
```

Classifier failure policies are independent:

```yaml
kraken2:
  confidence: 0.0
  failure_policy: continue

metabuli:
  min_score: 0.008
  max_ram_gb: 120
  failure_policy: continue

minimap2:
  reference: /absolute/path/to/classification_reference.cleaned.fa
  index: ""
  min_mapq: 15
  min_alignment: 500
  failure_policy: continue
```

### 4. Validate

```bash
bash scripts/start_new_dataset.sh \
    --action validate \
    --config config/my_new_run.yaml
```

Configuration loading checks core inputs. The workflow preflight records a
separate readiness assessment for each classifier so an unavailable method is
visible without making healthy classifiers unusable.

### 5. Dry-run the complete DAG

```bash
bash scripts/start_new_dataset.sh \
    --action dry-run \
    --config config/my_new_run.yaml \
    --profile workflow/profiles/dundee \
    --jobs 10
```

### 6. Start or resume

```bash
bash scripts/start_new_dataset.sh \
    --action run \
    --config config/my_new_run.yaml \
    --profile workflow/profiles/dundee \
    --jobs 10
```

The generic helper does not assume KmerSutra needs `4week`. If a new dataset
also requires it, use the explicit resource overrides with `run_workflow.sh`:

```bash
bash scripts/run_workflow.sh \
    --config config/my_new_run.yaml \
    --profile workflow/profiles/dundee \
    --jobs 10 \
    --set-resource classify_kmersutra:runtime=40320 \
    --set-resource classify_kmersutra:slurm_qos=4week \
    --verbose
```

The YAML `kmersutra.timeout_minutes` must remain below the scheduler runtime.

## Configuration notes

### Read state

`inputs.read_state` is provenance, not a convenience switch.

- Use `raw` only when this workflow must remove the configured host.
- Use `host_removed` only when the supplied FASTQs were already depleted.

The `01_host_depletion` directory is retained as a compatibility stage name for
both modes. In host-removed mode its summaries explicitly state
`host_depletion_performed=false`.

### Controlled minimap2 classification

The minimap2 classification reference is not the host reference. Mapping uses
the ONT preset, retains secondary PAF alignments for ambiguity assessment, and
applies configured MAPQ and alignment-length thresholds. Cross-taxon best-hit
ties remain ambiguous rather than being forced into a unique taxon.

If a supplied `.mmi` cannot be proven to match the configured FASTA, leave the
index blank and let the workflow build one from that exact reference.

### Reporting focus

The full TSV and JSON evidence remains complete. `reporting.max_table_rows`
limits only the number of rows embedded in each browser table. `focus_taxa`
contains case-insensitive taxon-name fragments highlighted in final,
classifier, comparison and sample reports.

```yaml
reporting:
  focus_taxa:
    - Plasmodium
  top_n: 20
  max_table_rows: 5000
```

## Output layout

```text
00_preflight/
  preflight.json
  resolved_samples.tsv
  software_versions.tsv
  classifier_readiness.tsv
  classification_reference.mmi
  classification_reference_index.complete.json

01_host_depletion/
  stage.complete.json
  <sample_id>/
    non_host.fastq.gz
    host_removal_summary.tsv
    metadata.json
    complete.json

02_classification/
  kraken2/
    stage.complete.json
    <sample_id>/
  metabuli/
    stage.complete.json
    <sample_id>/
  minimap2/
    stage.complete.json
    <sample_id>/
  kmersutra/
    <sample_id>/
      stage.complete.json

03_final/
  sample_summary.tsv
  classifier_status.tsv
  classifier_taxon_reports.tsv.gz
  minimap2_taxon_reports.tsv.gz
  kmersutra_species_calls.tsv.gz
  normalised_classifier_evidence.tsv.gz
  report_warnings.tsv
  report_data.json
  report_manifest.json
  README.txt
  SHA256SUMS.tsv
  workflow.complete.json
  reports/
    index.html
    comparison.html
    classifiers/
      kraken2.html
      metabuli.html
      minimap2.html
      kmersutra.html
    samples/
      <sample_id>.html
```

## Reading the HTML reports

Start with:

```text
03_final/reports/index.html
```

The final report shows:

- overall completeness;
- success/failure status by sample and classifier;
- links to each classifier and sample report;
- configured focus-taxon evidence;
- non-fatal parsing warnings; and
- core run provenance.

Each classifier report shows its own sample status, top reported taxa and
searchable evidence table. The comparison report uses positive,
species-comparable evidence for descriptive Jaccard overlap and an evidence
matrix. It does not represent overlap as accuracy and does not use one
classifier as a truth set.

All pages are portable and can be copied from the cluster with the rest of
`03_final`. No web server is required.

## Restart and repair

Snakemake skips valid completed work. To resume after a classifier failure,
correct the resource, executable, memory or runtime problem and run the same
launcher again. Successful method/sample outputs are not recalculated.

If a prior `failure.json` exists and the rerun succeeds, the atomic publication
replaces the failed terminal state with validated completion metadata.

To rebuild reports manually without rerunning classifiers:

```bash
nanopore-realdata --action report --config config/my_run.yaml
```

If a classifier is intentionally no longer required, disable it in the config
where supported or leave its failure policy as `continue`; the final report
will expose the incomplete evidence rather than hiding it.

## Quality checks

Run the complete gate from the repository root:

```bash
bash scripts/run_quality_checks.sh \
    --results-dir /absolute/path/to/test_results
```

The gate includes:

- Ruff lint and formatting checks;
- unit and integration tests;
- branch-aware coverage with a 90% minimum;
- Python byte-code compilation;
- Bash syntax checking;
- package builds; and
- a genuine Snakemake DAG dry run with synthetic inputs and resources.

## Scientific interpretation

Kraken2 and Metabuli taxonomic counts, controlled-reference minimap2
alignments and KmerSutra exact-k-mer evidence are different measurements.
Database and reference composition, software versions and thresholds can
materially alter apparent off-target calls. The workflow therefore retains
native tables, explicit provenance and a normalised descriptive layer while
avoiding a false consensus.
