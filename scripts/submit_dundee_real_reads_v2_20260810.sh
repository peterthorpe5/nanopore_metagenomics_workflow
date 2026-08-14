#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Submit the detached Dundee real-read classifier benchmark." \
        "" \
        "Options:" \
        "  --plan               Print the job arrays and dependencies only" \
        "  --resume-submission  Resume a submission interrupted between sbatch calls" \
        "  --new-attempt        Submit a new attempt after all previous jobs ended" \
        "  --retry-method NAME  Retry only one failed classifier plus aggregation" \
        "  --help"
}

ACTION="submit-slurm"
SUBMISSION_FLAG=""
RETRY_METHODS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --plan)
            ACTION="plan-slurm"
            shift
            ;;
        --resume-submission)
            [[ -z "${SUBMISSION_FLAG}" ]] || {
                printf 'ERROR: choose only one submission recovery flag\n' >&2
                exit 2
            }
            SUBMISSION_FLAG="--resume-submission"
            shift
            ;;
        --new-attempt)
            [[ -z "${SUBMISSION_FLAG}" ]] || {
                printf 'ERROR: choose only one submission recovery flag\n' >&2
                exit 2
            }
            SUBMISSION_FLAG="--new-attempt"
            shift
            ;;
        --retry-method)
            [[ $# -ge 2 ]] || {
                printf 'ERROR: --retry-method requires a classifier name\n' >&2
                exit 2
            }
            case "$2" in
                kraken2|metabuli|minimap2|kmersutra)
                    RETRY_METHODS+=("$2")
                    ;;
                *)
                    printf 'ERROR: unsupported retry method: %s\n' "$2" >&2
                    exit 2
                    ;;
            esac
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

if [[ "${ACTION}" == "plan-slurm" && -n "${SUBMISSION_FLAG}" ]]; then
    printf 'ERROR: --plan cannot be combined with a submission recovery flag\n' >&2
    exit 2
fi

if [[ -n "${SUBMISSION_FLAG}" && ${#RETRY_METHODS[@]} -gt 0 ]]; then
    printf 'ERROR: --retry-method cannot be combined with a submission recovery flag\n' >&2
    exit 2
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${REPOSITORY_ROOT}/config/dundee_real_reads_v2_20260810.yaml"
CONDA_ENVIRONMENT="nanopore_realdata_workflow"

[[ -f "${CONFIG_PATH}" ]] || {
    printf 'ERROR: workflow configuration is missing: %s\n' "${CONFIG_PATH}" >&2
    exit 2
}

COMMAND=(
    conda run
    --no-capture-output
    --name "${CONDA_ENVIRONMENT}"
    env "PYTHONPATH=${REPOSITORY_ROOT}/src"
    python -m nanopore_realdata.cli
    --action "${ACTION}"
    --config "${CONFIG_PATH}"
    --verbose
)
if [[ -n "${SUBMISSION_FLAG}" ]]; then
    COMMAND+=("${SUBMISSION_FLAG}")
fi
for RETRY_METHOD in "${RETRY_METHODS[@]}"; do
    COMMAND+=(--retry-method "${RETRY_METHOD}")
done

printf 'Repository: %s\n' "${REPOSITORY_ROOT}"
printf 'Configuration: %s\n' "${CONFIG_PATH}"
printf 'Mode: %s\n' "${ACTION}"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'
"${COMMAND[@]}"
