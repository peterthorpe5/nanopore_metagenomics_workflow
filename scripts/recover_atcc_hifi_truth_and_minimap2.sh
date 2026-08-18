#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow"
CONFIG_PATH="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_srr9328980_20260818.yaml"
RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/benchmarks/comparators_atcc_msa1003_hifi/atcc_msa1003_hifi_comparators_srr9328980_20260817"
CONDA_ENVIRONMENT="nanopore_realdata_workflow"
STAGE_RUNNER="${REPOSITORY_ROOT}/workflow/slurm/run_stage.sh"
LOG_DIRECTORY="${RUN_ROOT}/workflow_control/slurm_logs"
BACKUP_DIRECTORY="${RUN_ROOT}/03_final.before_truth_fix_20260818"

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${path}" >&2
        exit 2
    }
}

require_file "${CONFIG_PATH}"
require_file "${STAGE_RUNNER}"
require_file "${RUN_ROOT}/02_classification/kraken2/SRR9328980/complete.json"
require_file "${RUN_ROOT}/02_classification/metabuli/SRR9328980/complete.json"

mkdir -p "${LOG_DIRECTORY}"
if [[ -d "${RUN_ROOT}/03_final" && ! -e "${BACKUP_DIRECTORY}" ]]; then
    cp -a "${RUN_ROOT}/03_final" "${BACKUP_DIRECTORY}"
fi

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    env "PYTHONPATH=${REPOSITORY_ROOT}/src" \
    python -m nanopore_realdata.cli \
    --action validate \
    --config "${CONFIG_PATH}" \
    --verbose

REFERENCE_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_ref_fix \
        --cpus-per-task=12 \
        --mem=160000M \
        --time=2880 \
        --signal=B:TERM@300 \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${STAGE_RUNNER}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --action build-minimap-reference \
        --config "${CONFIG_PATH}"
)"

INDEX_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_idx_fix \
        --cpus-per-task=12 \
        --mem=160000M \
        --time=2880 \
        --signal=B:TERM@300 \
        --dependency="afterok:${REFERENCE_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${STAGE_RUNNER}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --action build-minimap-index \
        --config "${CONFIG_PATH}"
)"

MINIMAP_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_mini_fix \
        --cpus-per-task=12 \
        --mem=160000M \
        --time=2880 \
        --signal=B:TERM@300 \
        --array=0-0%1 \
        --dependency="afterok:${INDEX_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${STAGE_RUNNER}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --action classify \
        --config "${CONFIG_PATH}" \
        --method minimap2 \
        --sample-index-from-slurm
)"

AGGREGATE_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_agg_fix \
        --cpus-per-task=2 \
        --mem=8192M \
        --time=180 \
        --signal=B:TERM@300 \
        --dependency="afterany:${MINIMAP_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${STAGE_RUNNER}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --action aggregate \
        --config "${CONFIG_PATH}"
)"

printf 'Submitted minimap2 recovery chain:\n'
printf '  reference: %s\n' "${REFERENCE_JOB_ID}"
printf '  index:     %s\n' "${INDEX_JOB_ID}"
printf '  minimap2:  %s\n' "${MINIMAP_JOB_ID}"
printf '  aggregate: %s\n' "${AGGREGATE_JOB_ID}"
printf 'Kraken2 and Metabuli were not resubmitted.\n'
