# Deployment runbook: v0.4.0 real-read benchmark

Release date: 10 August 2026

This runbook installs the workflow from a versioned archive without a GitHub
push and submits a detached Slurm DAG. It does not alter the running ATCC job or
the previous real-read result directory.

## 1. Copy the release from the Mac

After downloading both release files, run on the Mac:

```bash
RELEASE=nanopore_metagenomics_workflow_v0_4_0_real_reads_v2_20260810.tar.gz
PROJECT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity

scp "${HOME}/Downloads/${RELEASE}" \
    "${HOME}/Downloads/${RELEASE}.sha256" \
    "pthorpe001@<dundee-login-host>:${PROJECT}/"
```

Replace only `<dundee-login-host>` with the normal SSH host name.

## 2. Verify and extract on Dundee

Run on a Dundee login node:

```bash
set -Eeuo pipefail

PROJECT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity
RELEASE=nanopore_metagenomics_workflow_v0_4_0_real_reads_v2_20260810.tar.gz
REPOSITORY=${PROJECT}/nanopore_metagenomics_workflow_v0_4_0

cd "${PROJECT}"
sha256sum --check "${RELEASE}.sha256"

if [[ -e "${REPOSITORY}" ]]; then
    printf 'STOP: release path already exists: %s\n' "${REPOSITORY}" >&2
    exit 3
fi

tar --extract --gzip --file "${RELEASE}"
test -f "${REPOSITORY}/pyproject.toml"
test -x "${REPOSITORY}/workflow/slurm/run_stage.sh"
```

The refusal to overwrite an existing release prevents a partial extraction from
silently replacing a working copy.

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
QUALITY_ROOT=${PROJECT}/quality/nanopore_realdata_v0_4_0

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
