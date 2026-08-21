#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
CONFIG_PATH=""
RUN_ROOT=""
ARCHIVE_SUFFIX=""

usage() {
    printf '%s\n' \
        'Usage: run_atcc_matched_kraken2_aggregation.sh [named options]' \
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
require_file "${RUN_ROOT}/02_classification/kraken2/SRR9328980/complete.json"
require_file "${RUN_ROOT}/02_classification/metabuli/SRR9328980/complete.json"

CURRENT_FINAL="${RUN_ROOT}/03_final"
ARCHIVED_FINAL="${RUN_ROOT}/03_final.${ARCHIVE_SUFFIX}"
if [[ -e "${ARCHIVED_FINAL}" ]]; then
    printf 'ERROR: aggregate archive target already exists: %s\n' "${ARCHIVED_FINAL}" >&2
    exit 2
fi
if [[ -d "${CURRENT_FINAL}" ]]; then
    mv "${CURRENT_FINAL}" "${ARCHIVED_FINAL}"
    printf 'Archived superseded aggregate: %s\n' "${ARCHIVED_FINAL}"
fi

"${STAGE_RUNNER}" \
    --repository-root "${REPOSITORY_ROOT}" \
    --conda-environment "${CONDA_ENVIRONMENT}" \
    --action aggregate \
    --config "${CONFIG_PATH}"

require_file "${CURRENT_FINAL}/workflow.complete.json"
require_file "${CURRENT_FINAL}/pcr_method_summary.tsv"
printf 'Corrected aggregation completed: %s\n' "${CURRENT_FINAL}/workflow.complete.json"
