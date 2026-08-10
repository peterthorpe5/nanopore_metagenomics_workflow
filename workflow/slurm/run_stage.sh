#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
ACTION=""
CONFIG_PATH=""
METHOD=""
SAMPLE_INDEX_FROM_SLURM=false

usage() {
    printf '%s\n' \
        "Usage: run_stage.sh --repository-root PATH --conda-environment NAME" \
        "       --action ACTION --config PATH [--method METHOD --sample-index-from-slurm]"
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
        --action)
            ACTION="${2:?Missing value for --action}"
            shift 2
            ;;
        --config)
            CONFIG_PATH="${2:?Missing value for --config}"
            shift 2
            ;;
        --method)
            METHOD="${2:?Missing value for --method}"
            shift 2
            ;;
        --sample-index-from-slurm)
            SAMPLE_INDEX_FROM_SLURM=true
            shift
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
[[ -f "${CONFIG_PATH}" ]] || {
    printf 'ERROR: configuration is missing: %s\n' "${CONFIG_PATH}" >&2
    exit 2
}
[[ -n "${CONDA_ENVIRONMENT}" ]] || {
    printf 'ERROR: --conda-environment is required\n' >&2
    exit 2
}
[[ -n "${ACTION}" ]] || {
    printf 'ERROR: --action is required\n' >&2
    exit 2
}

SAMPLE_INDEX=""
if [[ "${SAMPLE_INDEX_FROM_SLURM}" == true ]]; then
    SAMPLE_INDEX="${SLURM_ARRAY_TASK_ID:-}"
    [[ "${SAMPLE_INDEX}" =~ ^[0-9]+$ ]] || {
        printf 'ERROR: SLURM_ARRAY_TASK_ID is unavailable or invalid\n' >&2
        exit 2
    }
    [[ -n "${METHOD}" ]] || {
        printf 'ERROR: --method is required for an array classifier task\n' >&2
        exit 2
    }
fi

BASE_COMMAND=(
    conda run
    --no-capture-output
    --name "${CONDA_ENVIRONMENT}"
    env "PYTHONPATH=${REPOSITORY_ROOT}/src"
    python -m nanopore_realdata.cli
)

record_outer_failure() {
    local exit_code="$1"
    local reason="$2"
    trap - ERR TERM INT
    set +e
    if [[ "${ACTION}" == "classify" \
        && -n "${METHOD}" \
        && -n "${SAMPLE_INDEX}" ]]; then
        "${BASE_COMMAND[@]}" \
            --action record-scheduler-failure \
            --config "${CONFIG_PATH}" \
            --method "${METHOD}" \
            --sample-index "${SAMPLE_INDEX}" \
            --message "${reason}; wrapper exit code ${exit_code}" \
            --status scheduler_failed \
            --verbose
    fi
    exit "${exit_code}"
}

trap 'record_outer_failure "$?" "stage command failed"' ERR
trap 'record_outer_failure 143 "Slurm sent TERM before the allocation ended"' TERM
trap 'record_outer_failure 130 "stage was interrupted"' INT

COMMAND=(
    "${BASE_COMMAND[@]}"
    --action "${ACTION}"
    --config "${CONFIG_PATH}"
    --verbose
)
if [[ -n "${METHOD}" ]]; then
    COMMAND+=(--method "${METHOD}")
fi
if [[ -n "${SAMPLE_INDEX}" ]]; then
    COMMAND+=(--sample-index "${SAMPLE_INDEX}")
fi

printf 'Started UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'Host: %s\n' "$(hostname)"
printf 'Slurm job: %s array task: %s\n' \
    "${SLURM_JOB_ID:-not_slurm}" \
    "${SLURM_ARRAY_TASK_ID:-not_array}"
printf 'Repository: %s\n' "${REPOSITORY_ROOT}"
printf 'Configuration: %s\n' "${CONFIG_PATH}"
printf 'Command:'
printf ' %q' "${COMMAND[@]}"
printf '\n'

"${COMMAND[@]}"

printf 'Completed UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
