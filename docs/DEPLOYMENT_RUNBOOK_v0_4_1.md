# Deployment runbook: v0.4.1 real-read workflow and benchmark

Release date: 10 August 2026

This runbook follows the authoritative source route: Mac checkout to GitHub,
then GitHub to the cluster checkout. Source archives must not be copied or
extracted directly onto the cluster. The workflow submits a detached Slurm DAG
and does not alter the running ATCC job or the previous real-read result
directory.

## 1. Validate and push from the Mac

Apply the reviewed changes to the Mac checkout. Run the local quality gate and
inspect the exact changes before committing and pushing them yourself:

```bash
cd /path/to/nanopore_metagenomics_workflow

git status --short
git diff --check
bash scripts/run_quality_checks.sh \
    --results-dir /tmp/nanopore_realdata_quality \
    --run-label mac_acceptance_20260810

git diff
git add \
    README.md \
    config/dundee_real_reads_v2_20260810.yaml \
    config/real_data.template.yaml \
    docs/DEPLOYMENT_RUNBOOK_v0_4_1.md \
    docs/QUALITY_REPORT_v0_4_1.md \
    pyproject.toml \
    scripts/start_new_dataset.sh \
    src/nanopore_realdata/__init__.py \
    src/nanopore_realdata/config.py \
    src/nanopore_realdata/reference.py \
    src/nanopore_realdata/reporting.py \
    src/nanopore_realdata/slurm.py \
    src/nanopore_realdata/workflow.py \
    tests/helpers.py \
    tests/test_config_commands.py \
    tests/test_new_dataset_shell.py \
    tests/test_reference.py \
    tests/test_workflow_cli.py \
    tests/test_slurm.py
git commit -m "Make real-sample workflow reusable and fix deployment path"
git push origin main
```

Do not copy repository files directly to GPFS. GitHub is the transfer and audit
boundary.

## 2. Fast-forward the Dundee checkout

After the GitHub push has completed, run on a Dundee login node:

```bash
set -Eeuo pipefail

PROJECT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity
REPOSITORY=${PROJECT}/nanopore_metagenomics_workflow

cd "${REPOSITORY}"
ORIGIN_URL="$(git remote get-url origin)"
case "${ORIGIN_URL}" in
    "https://github.com/peterthorpe5/nanopore_metagenomics_workflow"|\
    "https://github.com/peterthorpe5/nanopore_metagenomics_workflow.git"|\
    "git@github.com:peterthorpe5/nanopore_metagenomics_workflow"|\
    "git@github.com:peterthorpe5/nanopore_metagenomics_workflow.git")
        ;;
    *)
        printf 'STOP: unexpected Git origin: %s\n' "${ORIGIN_URL}" >&2
        exit 3
        ;;
esac

if [[ -n "$(git status --porcelain)" ]]; then
    printf 'STOP: cluster checkout has local changes; inspect them before pulling.\n' >&2
    exit 3
fi

git fetch --prune origin
git pull --ff-only origin main
git rev-parse HEAD

test -f "${REPOSITORY}/pyproject.toml"
test -x "${REPOSITORY}/workflow/slurm/run_stage.sh"
```

The clean-tree and fast-forward-only checks prevent cluster-side edits or a
divergent checkout from being silently overwritten.

## 3. Update the environment

```bash
conda env update \
    --name nanopore_realdata_workflow \
    --file "${REPOSITORY}/environment/environment.yml"

conda run --name nanopore_realdata_workflow \
    python -m pip install --no-deps --editable "${REPOSITORY}"
```

Do not use `--prune` here: KmerSutra is installed separately in the same
environment.

Confirm the expected executables:

```bash
conda run --name nanopore_realdata_workflow kmersutra-screen --help >/dev/null
conda run --name nanopore_realdata_workflow kraken2 --version
conda run --name nanopore_realdata_workflow metabuli --help >/dev/null
conda run --name nanopore_realdata_workflow minimap2 --version
conda run --name nanopore_realdata_workflow snakemake --version
```

Confirm that the reduced masked reference selected for minimap2 is visible on
the login node. Do not substitute the much larger shared reference:

```bash
MASKED_REFERENCE=/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta
test -s "${MASKED_REFERENCE}"
ls -lh "${MASKED_REFERENCE}"
```

Preflight and index construction independently recheck the reference and refuse
the minimap2 branch unless its headers identify both `Plasmodium inui` and
`Plasmodium cynomolgi`.

## 4. Run the release gate

```bash
QUALITY_ROOT=${PROJECT}/quality/nanopore_realdata_v0_4_1

conda run --no-capture-output --name nanopore_realdata_workflow \
    bash "${REPOSITORY}/scripts/run_quality_checks.sh" \
    --results-dir "${QUALITY_ROOT}" \
    --run-label cluster_acceptance_20260810
```

Do not submit the benchmark if this command exits non-zero.

## 5. Inspect the detached plan

```bash
cd "${REPOSITORY}"
bash scripts/submit_dundee_real_reads_v2_20260810.sh --plan \
    | tee "${PROJECT}/real_reads_v2_submission_plan_20260810.json"
```

Check that the plan contains four classifier arrays over sample indices 0–10,
that KmerSutra alone uses `4week`, and that KmerSutra has no minimap2, Kraken2
or Metabuli dependency. With the prebuilt masked FASTA, the plan has eight jobs:
preflight, input acceptance, minimap2 index, four arrays and aggregation. It
must not contain a controlled-reference build job.

## 6. Submit once

```bash
cd "${REPOSITORY}"
bash scripts/submit_dundee_real_reads_v2_20260810.sh
```

The command returns after Slurm accepts the jobs. Closing SSH does not stop the
workflow.

The output and journal locations are:

```bash
OUTPUT=${PROJECT}/runs/real_reads_classifier_benchmark_v2_20260810
JOURNAL=${OUTPUT}/workflow_control/slurm_submission.json

python -m json.tool "${JOURNAL}"
squeue --user "${USER}" --name 'NRD_*'
```

Slurm job output is written beneath
`${OUTPUT}/workflow_control/slurm_logs/`. Resource-staging logs are kept beneath
`${OUTPUT}/workflow_control/resource_staging_logs/`.

## 7. Recovery rules

If the submitting shell was interrupted before every `sbatch` call was
journalled:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh --resume-submission
```

If the journal says `submitted`, do not use `--resume-submission`. Wait until
all recorded jobs have ended. Then a retry may be submitted with:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh --new-attempt
```

The new attempt reuses checksum-bound completion records and retries unfinished
samples. The submitter refuses a new attempt while any recorded job ID is still
active.

Do not delete `workflow_control`, completion JSON files or successful method
directories to make a retry. Those records are the restart and audit evidence.

## 8. Acceptance outputs

After the aggregate job ends, inspect:

```bash
FINAL=${OUTPUT}/03_final

python -m json.tool "${FINAL}/workflow.complete.json"
column -t -s $'\t' "${FINAL}/pcr_method_summary.tsv"
column -t -s $'\t' "${FINAL}/classifier_status.tsv" | less -S
```

Open `${FINAL}/reports/index.html` locally or through the usual cluster file
transfer route. The HTML is self-contained and needs no internet connection.

An overall `partial` result is valid when the status and PCR denominator tables
identify the unavailable branches. It must not be described as a complete
four-method comparison until all 44 sample–method runs are `success`.

## 9. Start later real-sample datasets

Create a new dataset configuration and manifest in the Mac checkout. PCR is
optional and is not created by default:

```bash
bash scripts/start_new_dataset.sh \
    --action initialise \
    --config config/new_real_samples.yaml
```

For a future benchmark with an independent PCR interpretation, add
`--with-pcr-truth`. Populate the manifest and YAML on the Mac, commit and push
them through GitHub, then fast-forward the clean Dundee checkout. Do not edit
the tracked configuration directly on the cluster.

On Dundee, validate and inspect the dataset-specific detached plan:

```bash
bash scripts/start_new_dataset.sh \
    --action validate \
    --config config/new_real_samples.yaml

bash scripts/start_new_dataset.sh \
    --action plan-slurm \
    --config config/new_real_samples.yaml
```

For host-removed inputs, submit with `--action submit-slurm`. The manifest
determines the number of samples and array tasks; no benchmark sample ID or
fixed sample count is coded into the workflow. Use
`--action resume-submission` only for an interrupted submission journal, and
`--action new-attempt` only after all earlier jobs have ended.

For raw FASTQs, configure `host.reference` and use the Snakemake `dry-run` and
`run` actions. Version 0.4.1 does not submit raw-read host depletion through the
detached backend; `plan-slurm` fails early instead of silently using an unsafe
bundled host-depletion job.
