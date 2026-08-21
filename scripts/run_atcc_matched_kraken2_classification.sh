#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
CONFIG_PATH=""
RUN_ROOT=""
ARCHIVE_SUFFIX=""

usage() {
    printf '%s\n' \
        'Usage: run_atcc_matched_kraken2_classification.sh [named options]' \
        '  --repository-root PATH' \
        '  --conda-environment NAME' \
        '  --config PATH' \
        '  --run-root PATH' \
        '  --archive-suffix TEXT'
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
        --config)
            CONFIG_PATH="${2:?Missing value for --config}"
            shift 2
            ;;
        --run-root)
            RUN_ROOT="${2:?Missing value for --run-root}"
            shift 2
            ;;
        --archive-suffix)
            ARCHIVE_SUFFIX="${2:?Missing value for --archive-suffix}"
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

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${path}" >&2
        exit 2
    }
}

[[ -d "${REPOSITORY_ROOT}" ]] || {
    printf 'ERROR: repository root is missing: %s\n' "${REPOSITORY_ROOT}" >&2
    exit 2
}
[[ -d "${RUN_ROOT}" ]] || {
    printf 'ERROR: run root is missing: %s\n' "${RUN_ROOT}" >&2
    exit 2
}
[[ -n "${CONDA_ENVIRONMENT}" ]] || {
    printf 'ERROR: --conda-environment is required\n' >&2
    exit 2
}
[[ "${ARCHIVE_SUFFIX}" =~ ^[A-Za-z0-9_.-]+$ ]] || {
    printf 'ERROR: --archive-suffix must be filesystem safe\n' >&2
    exit 2
}

STAGE_RUNNER="${REPOSITORY_ROOT}/workflow/slurm/run_stage.sh"
require_file "${CONFIG_PATH}"
require_file "${STAGE_RUNNER}"
require_file "${RUN_ROOT}/02_classification/metabuli/SRR9328980/complete.json"

CURRENT_RESULT="${RUN_ROOT}/02_classification/kraken2/SRR9328980"
ARCHIVED_RESULT="${RUN_ROOT}/02_classification/kraken2/SRR9328980.${ARCHIVE_SUFFIX}"
if [[ -e "${ARCHIVED_RESULT}" ]]; then
    printf 'ERROR: archive target already exists: %s\n' "${ARCHIVED_RESULT}" >&2
    exit 2
fi
if [[ -d "${CURRENT_RESULT}" ]]; then
    mv "${CURRENT_RESULT}" "${ARCHIVED_RESULT}"
    printf 'Archived superseded Kraken2 result: %s\n' "${ARCHIVED_RESULT}"
fi

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    env "PYTHONPATH=${REPOSITORY_ROOT}/src" \
    python -m nanopore_realdata.cli \
    --action validate \
    --config "${CONFIG_PATH}" \
    --verbose

export SLURM_ARRAY_TASK_ID=0
"${STAGE_RUNNER}" \
    --repository-root "${REPOSITORY_ROOT}" \
    --conda-environment "${CONDA_ENVIRONMENT}" \
    --action classify \
    --config "${CONFIG_PATH}" \
    --method kraken2 \
    --sample-index-from-slurm

require_file "${CURRENT_RESULT}/complete.json"
if [[ -s "${CURRENT_RESULT}/failure.json" ]]; then
    printf 'ERROR: Kraken2 wrote failure.json despite stage completion\n' >&2
    exit 2
fi
printf 'Matched Kraken2 classification completed: %s\n' "${CURRENT_RESULT}/complete.json"
