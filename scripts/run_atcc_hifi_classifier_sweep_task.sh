#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
METHOD=""
MANIFEST=""
INPUT_FASTQ=""
DATABASE=""
OUTPUT_ROOT=""
THREADS=""
MAX_RAM_GB="120"

usage() {
    printf '%s\n' \
        'Usage: run_atcc_hifi_classifier_sweep_task.sh [named options]' \
        '  --repository-root PATH' \
        '  --conda-environment NAME' \
        '  --method kraken2|metabuli' \
        '  --manifest PATH' \
        '  --input-fastq PATH' \
        '  --database PATH' \
        '  --output-root PATH' \
        '  --threads INTEGER' \
        '  --max-ram-gb INTEGER'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repository-root)
            REPOSITORY_ROOT="${2:?Missing value for --repository-root}"
            shift 2
            ;;
        --conda-environment)
            CONDA_ENVIRONMENT="${2:?Missing value for --conda-environment}"
            shift 2
            ;;
        --method)
            METHOD="${2:?Missing value for --method}"
            shift 2
            ;;
        --manifest)
            MANIFEST="${2:?Missing value for --manifest}"
            shift 2
            ;;
        --input-fastq)
            INPUT_FASTQ="${2:?Missing value for --input-fastq}"
            shift 2
            ;;
        --database)
            DATABASE="${2:?Missing value for --database}"
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT="${2:?Missing value for --output-root}"
            shift 2
            ;;
        --threads)
            THREADS="${2:?Missing value for --threads}"
            shift 2
            ;;
        --max-ram-gb)
            MAX_RAM_GB="${2:?Missing value for --max-ram-gb}"
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

[[ -d "${REPOSITORY_ROOT}" ]] || {
    printf 'ERROR: repository root is missing: %s\n' "${REPOSITORY_ROOT}" >&2
    exit 2
}
[[ "${METHOD}" == kraken2 || "${METHOD}" == metabuli ]] || {
    printf 'ERROR: --method must be kraken2 or metabuli\n' >&2
    exit 2
}
[[ -n "${CONDA_ENVIRONMENT}" ]] || {
    printf 'ERROR: --conda-environment is required\n' >&2
    exit 2
}
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: --threads must be a positive integer\n' >&2
    exit 2
}
[[ "${MAX_RAM_GB}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: --max-ram-gb must be a positive integer\n' >&2
    exit 2
}
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
    printf 'ERROR: SLURM_ARRAY_TASK_ID is not set\n' >&2
    exit 2
}

SWEEP_SCRIPT="${REPOSITORY_ROOT}/scripts/atcc_hifi_classifier_sweep.py"
for REQUIRED_FILE in "${MANIFEST}" "${INPUT_FASTQ}" "${SWEEP_SCRIPT}"; do
    [[ -s "${REQUIRED_FILE}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${REQUIRED_FILE}" >&2
        exit 2
    }
done
[[ -d "${DATABASE}" ]] || {
    printf 'ERROR: classifier database is missing: %s\n' "${DATABASE}" >&2
    exit 2
}

printf 'Started UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'Host: %s\n' "$(hostname --fqdn)"
printf 'Method: %s\n' "${METHOD}"
printf 'Array task: %s\n' "${SLURM_ARRAY_TASK_ID}"

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    python "${SWEEP_SCRIPT}" \
    --action run \
    --manifest "${MANIFEST}" \
    --method "${METHOD}" \
    --task-index "${SLURM_ARRAY_TASK_ID}" \
    --input-fastq "${INPUT_FASTQ}" \
    --database "${DATABASE}" \
    --output-root "${OUTPUT_ROOT}" \
    --threads "${THREADS}" \
    --max-ram-gb "${MAX_RAM_GB}" \
    --verbose

printf 'Completed UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
