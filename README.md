# Nanopore real-read classifier workflow

Version 0.4.1 is a restartable workflow for real Oxford Nanopore FASTQ data.
It runs Kraken2, Metabuli, masked-reference minimap2 and KmerSutra as
independent classifier branches. An independent PCR truth table can be added
for a benchmark, but is not required for routine real-sample classification.
The workflow does not force a consensus call and is not a clinical diagnostic
pipeline.

The repository is reusable across datasets. The prepared Dundee deployment is
one frozen 11-sample host-removed MRC benchmark preset in
`config/dundee_real_reads_v2_20260810.yaml`. Ten samples have species-level PCR
interpretations. `MRC1123_WuH001WB` remains in the classifier run but is marked
PCR `unknown` and excluded from the primary comparison.

## Why v0.4.1 exists

The previous partial run revealed two structural problems:

- several samples were bundled into one long job per classifier, and KmerSutra
  was downstream of the other classifier branches; and
- minimap2 used a very large mixed reference that split into about 30 index
  parts. Queries were mapped once per part, inflating mapped-read counts.

The v0.4 series removes both failure modes:

- every classifier is a separate per-sample Slurm array;
- KmerSutra depends only on prepared reads, never on another classifier;
- the shell that submits the DAG exits after `sbatch`, so there is no
  long-running Snakemake controller to lose;
- all dependencies use failure-aware completion and final aggregation waits for
  terminal classifier jobs even when some fail;
- outer Slurm failures are recorded per sample where the wrapper remains alive;
- missing hard-kill records remain `missing`, never biological negatives; and
- submission is journalled after every accepted Slurm job ID.

## Reduced masked-reference minimap2 contract

The prepared run uses the focused masked Plasmodium/outgroup FASTA already used
by the spike-in project:

```text
/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta
```

This is the exact default behind `MASKED_MINIMAP_DB_FASTA_DEFAULT`; it replaces
the impractically large shared bacterial/viral/Plasmodium FASTA. The workflow
does not reuse an old MMI. It fingerprints the FASTA and builds a checksum-bound
run-local index.

The run refuses to continue when any of these safeguards fails:

- the FASTA headers do not identify every species explicitly listed in
  `minimap2.required_species` (the prepared benchmark lists `Plasmodium inui`
  and `Plasmodium cynomolgi`);
- the reference exceeds 12,000,000,000 bases;
- minimap2 records anything other than one index part;
- the index exceeds 64,000,000,000 bytes;
- a PAF target is absent from the reference used to build the index;
- a query occurs in repeated non-consecutive PAF blocks; or
- mapped query groups exceed the validated input-read count.

The prepared run builds with `-I 16000000000`, which is larger than the hard
reference bound and therefore enforces a single index part.

## Detached Dundee DAG

The recommended backend is `nanopore_realdata.slurm`, not a persistent
Snakemake controller. The submitted stages are:

1. preflight, plus independent PCR validation only when a truth table is
   configured;
2. acceptance of already host-removed reads;
3. run-local single-part minimap2 index construction from the reduced masked
   reference;
4. four independent per-sample classifier arrays; and
5. failure-aware aggregation into TSV, JSON and offline HTML.

The prepared array limits are intentionally conservative:

| Method | Per-task resources | Maximum concurrent samples | QoS |
|---|---:|---:|---|
| Kraken2 | 12 CPUs, 96 GB | 1 | cluster default |
| Metabuli | 12 CPUs, 160 GB | 1 | cluster default |
| minimap2 | 12 CPUs, 160 GB | 2 | cluster default |
| KmerSutra | 24 CPUs, 96 GB | 1 | `4week` |

Only KmerSutra requests `4week`. Its internal timeout is 40,100 minutes, 220
minutes below the 40,320-minute Slurm allocation.

Snakemake remains available as a generic fallback through
`scripts/run_workflow.sh`, and its packaged DAG is also per sample. The prepared
Dundee benchmark should use the detached submitter.

## Installation

The authoritative source moves from the Mac checkout to GitHub and then to the
cluster checkout. Do not copy or extract release source directly onto the
cluster. The cluster repository must retain the canonical GitHub directory
name declared in the Dundee YAML:

```bash
PROJECT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity
REPOSITORY=${PROJECT}/nanopore_metagenomics_workflow

cd "${REPOSITORY}"
git status --short
git pull --ff-only
git remote get-url origin
git rev-parse HEAD

conda env update \
    --name nanopore_realdata_workflow \
    --file "${REPOSITORY}/environment/environment.yml"

conda run --name nanopore_realdata_workflow \
    python -m pip install --no-deps --editable "${REPOSITORY}"
```

KmerSutra 0.51.2 must already be installed in the same environment. The
preflight checks `kmersutra-screen`, Kraken2, Metabuli, minimap2, `pigz` and
`rsync` before scientific work begins.

The deployment identity guard refuses submission from a different repository
copy or package version. This prevents the mixed-checkout problem seen in the
earlier attempts.

## Prepared benchmark execution

First print the complete job plan without submitting anything. The plan now
validates the repository path and package version as well as the job graph:

```bash
cd /gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow
bash scripts/submit_dundee_real_reads_v2_20260810.sh --plan
```

Then submit once:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh
```

The command returns after the DAG is accepted by Slurm. It does not need
`tmux`, `screen` or an open SSH connection.

The durable submission journal is:

```text
runs/real_reads_classifier_benchmark_v2_20260810/workflow_control/slurm_submission.json
```

If the submitting shell is interrupted between `sbatch` calls, continue the
same journal with:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh --resume-submission
```

After every prior job has ended, a new attempt can safely reuse validated
completion records and retry unfinished work:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh --new-attempt
```

The submitter refuses an ambiguous repeat and checks `squeue` before allowing a
new attempt. It never deletes an earlier run.

## Final results

The final directory contains:

- `sample_summary.tsv` and `classifier_status.tsv`;
- native classifier tables and normalised method-specific evidence;
- `pcr_truth.tsv`, `pcr_concordance.tsv` and `pcr_method_summary.tsv` when PCR
  evaluation is configured; otherwise these retain schema headers only;
- checksums, provenance JSON and non-fatal reporting warnings; and
- a self-contained offline report at `reports/index.html`, including a dedicated
  PCR comparison page.

When PCR evaluation is configured, classifier failure or absence is
`unavailable`, never a negative call. All comparison counts expose their
available sample denominator.

## Starting another dataset

Create a protected YAML and sample manifest for routine real data:

```bash
bash scripts/start_new_dataset.sh \
    --action initialise \
    --config config/my_run.yaml
```

The helper refuses to overwrite either file. Edit every absolute path, set
`inputs.read_state` deliberately and validate before execution. Dataset size is
read from the manifest; no sample IDs or sample count are embedded in Python.

For a benchmark with independent PCR evidence, request the separate template:

```bash
bash scripts/start_new_dataset.sh \
    --action initialise \
    --config config/my_benchmark.yaml \
    --with-pcr-truth
```

Populate the truth TSV, set `inputs.pcr_truth` to its path, and explicitly list
all PCR-positive comparison species under `minimap2.required_species`.
PCR-unknown samples remain explicit with
`include_in_primary_comparison=false`. The truth table affects evaluation only;
it does not create classifier calls or change classifier thresholds.

Every future dataset whose manifest already contains host-removed reads can use
the same detached, failure-tolerant Slurm backend:

```bash
bash scripts/start_new_dataset.sh \
    --action plan-slurm \
    --config config/my_run.yaml

bash scripts/start_new_dataset.sh \
    --action submit-slurm \
    --config config/my_run.yaml
```

The number of array tasks comes from the manifest. Interrupted submissions can
use `--action resume-submission`; a later deliberate retry can use
`--action new-attempt` after the existing jobs have ended.

Raw-read configurations remain supported through the Snakemake path using
`--action dry-run` and `--action run` with a host reference. Version 0.4.1
deliberately refuses raw reads in the detached submitter because host depletion
is not yet implemented there as independent per-sample array tasks.

## Quality gate

Run the release checks from the repository root:

```bash
bash scripts/run_quality_checks.sh --results-dir /tmp/nanopore_realdata_quality
```

The gate runs Ruff, branch-aware coverage with a 90% minimum, all `unittest`
tests, Python compilation, package builds and Bash syntax validation. Synthetic
tests exercise failure continuation, atomic publication, Slurm dependency
planning, PCR exact matching, reference limits, multipart-index rejection and
impossible minimap2 read counts.

No software licence has yet been selected; see
`docs/LICENCE_SELECTION_REQUIRED.md` before public distribution.
