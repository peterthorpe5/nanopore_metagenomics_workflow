#!/usr/bin/env bash

set -Eeuo pipefail

MODE=""
KRAKEN2_JOB_ID=""

usage() {
    printf '%s\n' \
        'Usage: recover_atcc_hifi_metabuli_sweep.sh OPTIONS' \
        '' \
        'Required options:' \
        '  --mode plan|submit' \
        '  --kraken2-job-id JOB_ID' \
        '' \
        'This recovery submits only the five Metabuli operating points and' \
        'a summary job that also waits for the existing Kraken2 array.'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="${2:?Missing value for --mode}"
            shift 2
            ;;
        --kraken2-job-id)
            KRAKEN2_JOB_ID="${2:?Missing value for --kraken2-job-id}"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ "${MODE}" == plan || "${MODE}" == submit ]] || {
    printf 'ERROR: --mode must be plan or submit\n' >&2
    exit 2
}
[[ "${KRAKEN2_JOB_ID}" =~ ^[0-9]+$ ]] || {
    printf 'ERROR: --kraken2-job-id must be numeric\n' >&2
    exit 2
}

REPOSITORY_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow"
RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/benchmarks/comparators_atcc_msa1003_hifi/atcc_msa1003_hifi_comparators_srr9328980_20260817"
REFERENCE_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/databases/kmersutra_db/atcc_msa1003_heldout_v1"
KRAKEN2_DATABASE="${REFERENCE_ROOT}/kraken2_matched_475genomes_v1"
METABULI_DATABASE="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/databases/metabuli/custom_metabuli_db"
INPUT_FASTQ="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/benchmarks/kmersutra_atcc_msa1003/atcc_msa1003_hifi_srr9328980_heldout_v0512_20260731/stages/01_acquire_reads/fastq/SRR9328980.fastq.gz"
MANIFEST="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_classifier_operating_points_20260821.tsv"
TRUTH_MANIFEST="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_srr9328980_20260821.truth_species.tsv"
SWEEP_ROOT="${RUN_ROOT}/04_classifier_operating_point_sweep"
SUMMARY_DIRECTORY="${SWEEP_ROOT}/summary"
EXISTING_KRAKEN2_REPORT="${RUN_ROOT}/02_classification/kraken2/SRR9328980/report.tsv"
EXISTING_METABULI_REPORT="${RUN_ROOT}/02_classification/metabuli/SRR9328980/report.tsv"
CONDA_ENVIRONMENT="nanopore_realdata_workflow"
LOG_DIRECTORY="${RUN_ROOT}/workflow_control/slurm_logs"
SWEEP_SCRIPT="${REPOSITORY_ROOT}/scripts/atcc_hifi_classifier_sweep.py"
TASK_SCRIPT="${REPOSITORY_ROOT}/scripts/run_atcc_hifi_classifier_sweep_task.sh"
SUMMARY_SCRIPT="${REPOSITORY_ROOT}/scripts/run_atcc_hifi_classifier_sweep_summary.sh"

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${path}" >&2
        exit 2
    }
}

require_directory() {
    local path="$1"
    [[ -d "${path}" ]] || {
        printf 'ERROR: required directory is missing: %s\n' "${path}" >&2
        exit 2
    }
}

require_directory "${REPOSITORY_ROOT}"
require_directory "${RUN_ROOT}"
require_directory "${KRAKEN2_DATABASE}"
require_directory "${METABULI_DATABASE}"
require_file "${INPUT_FASTQ}"
require_file "${MANIFEST}"
require_file "${TRUTH_MANIFEST}"
require_file "${EXISTING_KRAKEN2_REPORT}"
require_file "${EXISTING_METABULI_REPORT}"
require_file "${SWEEP_SCRIPT}"
require_file "${TASK_SCRIPT}"
require_file "${SUMMARY_SCRIPT}"

METABULI_HELP="$(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
        metabuli classify -h 2>&1 || true
)"
for REQUIRED_OPTION in --seq-mode --precise --min-score --min-sp-score; do
    [[ "${METABULI_HELP}" == *"${REQUIRED_OPTION}"* ]] || {
        printf 'ERROR: the configured Metabuli lacks %s support\n' \
            "${REQUIRED_OPTION}" >&2
        exit 2
    }
done

python3 "${SWEEP_SCRIPT}" \
    --action validate-manifest \
    --manifest "${MANIFEST}" \
    --verbose

METABULI_TASK_COUNT="$(
    python3 "${SWEEP_SCRIPT}" \
        --action task-count \
        --manifest "${MANIFEST}" \
        --method metabuli
)"
[[ "${METABULI_TASK_COUNT}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: invalid Metabuli task count: %s\n' \
        "${METABULI_TASK_COUNT}" >&2
    exit 2
}

if [[ -e "${SUMMARY_DIRECTORY}" ]]; then
    printf 'ERROR: summary output already exists; refusing to overwrite: %s\n' \
        "${SUMMARY_DIRECTORY}" >&2
    exit 2
fi

mkdir -p "${LOG_DIRECTORY}" "${SWEEP_ROOT}"

printf 'ATCC HiFi Metabuli recovery plan:\n'
printf '  existing Kraken2 array: %s (not resubmitted)\n' "${KRAKEN2_JOB_ID}"
printf '  Metabuli sequence mode: 3 (long read)\n'
printf '  Metabuli operating points to run: %s\n' "${METABULI_TASK_COUNT}"
printf '  final summary waits for Kraken2 and Metabuli\n'
printf '  output root: %s\n' "${SWEEP_ROOT}"

if [[ "${MODE}" == plan ]]; then
    printf 'Plan checks passed. No jobs were submitted.\n'
    exit 0
fi

command -v sbatch >/dev/null 2>&1 || {
    printf 'ERROR: sbatch is unavailable\n' >&2
    exit 2
}

METABULI_LAST_INDEX="$((METABULI_TASK_COUNT - 1))"
METABULI_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_mbsweep_fix \
        --chdir="${REPOSITORY_ROOT}" \
        --array="0-${METABULI_LAST_INDEX}%2" \
        --cpus-per-task=12 \
        --mem=160000M \
        --time=2880 \
        --signal=B:TERM@300 \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${TASK_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --method metabuli \
        --manifest "${MANIFEST}" \
        --input-fastq "${INPUT_FASTQ}" \
        --database "${METABULI_DATABASE}" \
        --output-root "${SWEEP_ROOT}" \
        --threads 12 \
        --max-ram-gb 120
)"

SUMMARY_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_opsum_fix \
        --chdir="${REPOSITORY_ROOT}" \
        --cpus-per-task=2 \
        --mem=8192M \
        --time=180 \
        --signal=B:TERM@300 \
        --dependency="afterok:${KRAKEN2_JOB_ID}:${METABULI_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%j.out" \
        --error="${LOG_DIRECTORY}/%x.%j.err" \
        "${SUMMARY_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --manifest "${MANIFEST}" \
        --truth-manifest "${TRUTH_MANIFEST}" \
        --sweep-root "${SWEEP_ROOT}" \
        --existing-kraken2-report "${EXISTING_KRAKEN2_REPORT}" \
        --existing-metabuli-report "${EXISTING_METABULI_REPORT}" \
        --output-directory "${SUMMARY_DIRECTORY}"
)"

printf 'Submitted corrected Metabuli recovery:\n'
printf '  existing Kraken2 array: %s\n' "${KRAKEN2_JOB_ID}"
printf '  Metabuli array:         %s\n' "${METABULI_JOB_ID}"
printf '  replacement summary:    %s\n' "${SUMMARY_JOB_ID}"
printf 'Kraken2 and the existing baseline reports were not modified.\n'
