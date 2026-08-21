#!/usr/bin/env bash

set -Eeuo pipefail

MODE=""

usage() {
    printf '%s\n' \
        'Usage: submit_atcc_hifi_classifier_sweep.sh --mode plan|submit' \
        '' \
        'plan    validates every input and prints the intended Slurm jobs' \
        'submit  submits Kraken2, Metabuli and dependent summary jobs'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="${2:?Missing value for --mode}"
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
require_file "${KRAKEN2_DATABASE}/matched_database.complete.json"
require_file "${KRAKEN2_DATABASE}/hash.k2d"
require_file "${KRAKEN2_DATABASE}/opts.k2d"
require_file "${KRAKEN2_DATABASE}/taxo.k2d"
require_file "${INPUT_FASTQ}"
require_file "${MANIFEST}"
require_file "${TRUTH_MANIFEST}"
require_file "${EXISTING_KRAKEN2_REPORT}"
require_file "${EXISTING_METABULI_REPORT}"
require_file "${SWEEP_SCRIPT}"
require_file "${TASK_SCRIPT}"
require_file "${SUMMARY_SCRIPT}"

KRAKEN2_HELP="$(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
        kraken2 --help 2>&1 || true
)"
[[ "${KRAKEN2_HELP}" == *"--confidence"* ]] || {
    printf 'ERROR: the configured Kraken2 lacks --confidence support\n' >&2
    exit 2
}

METABULI_HELP="$(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
        metabuli classify -h 2>&1 || true
)"
[[ "${METABULI_HELP}" == *"--precise"* ]] || {
    printf 'ERROR: the configured Metabuli lacks the v1.2 --precise option\n' >&2
    exit 2
}
[[ "${METABULI_HELP}" == *"--min-sp-score"* ]] || {
    printf 'ERROR: the configured Metabuli lacks --min-sp-score support\n' >&2
    exit 2
}

python3 "${SWEEP_SCRIPT}" \
    --action validate-manifest \
    --manifest "${MANIFEST}" \
    --verbose

KRAKEN2_TASK_COUNT="$(
    python3 "${SWEEP_SCRIPT}" \
        --action task-count \
        --manifest "${MANIFEST}" \
        --method kraken2
)"
METABULI_TASK_COUNT="$(
    python3 "${SWEEP_SCRIPT}" \
        --action task-count \
        --manifest "${MANIFEST}" \
        --method metabuli
)"

[[ "${KRAKEN2_TASK_COUNT}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: invalid Kraken2 task count: %s\n' "${KRAKEN2_TASK_COUNT}" >&2
    exit 2
}
[[ "${METABULI_TASK_COUNT}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: invalid Metabuli task count: %s\n' "${METABULI_TASK_COUNT}" >&2
    exit 2
}

if [[ -e "${SUMMARY_DIRECTORY}" ]]; then
    printf 'ERROR: summary output already exists; refusing to overwrite: %s\n' \
        "${SUMMARY_DIRECTORY}" >&2
    exit 2
fi

mkdir -p "${LOG_DIRECTORY}" "${SWEEP_ROOT}"

printf 'ATCC HiFi classifier operating-point plan:\n'
printf '  matched Kraken2 database: %s\n' "${KRAKEN2_DATABASE}"
printf '  existing Kraken2 baseline: confidence 0.00\n'
printf '  new Kraken2 tasks: %s (0.05, 0.10, 0.20, 0.50)\n' "${KRAKEN2_TASK_COUNT}"
printf '  Metabuli database: %s\n' "${METABULI_DATABASE}"
printf '  existing Metabuli baseline: min-score 0.008, min-sp-score 0\n'
printf '  new Metabuli tasks: %s\n' "${METABULI_TASK_COUNT}"
printf '  explicit HiFi Metabuli thresholds: min-score 0.07, min-sp-score 0.3\n'
printf '  current Metabuli v1.2 HiFi mode: --precise 2\n'
printf '  output root: %s\n' "${SWEEP_ROOT}"

if [[ "${MODE}" == plan ]]; then
    printf 'Plan checks passed. No jobs were submitted.\n'
    exit 0
fi

command -v sbatch >/dev/null 2>&1 || {
    printf 'ERROR: sbatch is unavailable\n' >&2
    exit 2
}

KRAKEN2_LAST_INDEX="$((KRAKEN2_TASK_COUNT - 1))"
METABULI_LAST_INDEX="$((METABULI_TASK_COUNT - 1))"

KRAKEN2_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_k2sweep \
        --chdir="${REPOSITORY_ROOT}" \
        --array="0-${KRAKEN2_LAST_INDEX}%1" \
        --cpus-per-task=12 \
        --mem=409600M \
        --time=1440 \
        --signal=B:TERM@300 \
        --output="${LOG_DIRECTORY}/%x.%A_%a.out" \
        --error="${LOG_DIRECTORY}/%x.%A_%a.err" \
        "${TASK_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --method kraken2 \
        --manifest "${MANIFEST}" \
        --input-fastq "${INPUT_FASTQ}" \
        --database "${KRAKEN2_DATABASE}" \
        --output-root "${SWEEP_ROOT}" \
        --threads 12 \
        --max-ram-gb 120
)"

METABULI_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_mbsweep \
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
        --job-name=NRD_opsum \
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

printf 'Submitted ATCC HiFi classifier operating-point sweep:\n'
printf '  Kraken2 array: %s\n' "${KRAKEN2_JOB_ID}"
printf '  Metabuli array: %s\n' "${METABULI_JOB_ID}"
printf '  Summary job:    %s\n' "${SUMMARY_JOB_ID}"
printf 'Existing baseline results and both databases were preserved.\n'
