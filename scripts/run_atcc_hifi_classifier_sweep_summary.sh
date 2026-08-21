#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
MANIFEST=""
TRUTH_MANIFEST=""
SWEEP_ROOT=""
EXISTING_KRAKEN2_REPORT=""
EXISTING_METABULI_REPORT=""
OUTPUT_DIRECTORY=""

usage() {
    printf '%s\n' \
        'Usage: run_atcc_hifi_classifier_sweep_summary.sh [named options]' \
        '  --repository-root PATH' \
        '  --conda-environment NAME' \
        '  --manifest PATH' \
        '  --truth-manifest PATH' \
        '  --sweep-root PATH' \
        '  --existing-kraken2-report PATH' \
        '  --existing-metabuli-report PATH' \
        '  --output-directory PATH'
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
        --manifest)
            MANIFEST="${2:?Missing value for --manifest}"
            shift 2
            ;;
        --truth-manifest)
            TRUTH_MANIFEST="${2:?Missing value for --truth-manifest}"
            shift 2
            ;;
        --sweep-root)
            SWEEP_ROOT="${2:?Missing value for --sweep-root}"
            shift 2
            ;;
        --existing-kraken2-report)
            EXISTING_KRAKEN2_REPORT="${2:?Missing value for --existing-kraken2-report}"
            shift 2
            ;;
        --existing-metabuli-report)
            EXISTING_METABULI_REPORT="${2:?Missing value for --existing-metabuli-report}"
            shift 2
            ;;
        --output-directory)
            OUTPUT_DIRECTORY="${2:?Missing value for --output-directory}"
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
[[ -n "${CONDA_ENVIRONMENT}" ]] || {
    printf 'ERROR: --conda-environment is required\n' >&2
    exit 2
}

SWEEP_SCRIPT="${REPOSITORY_ROOT}/scripts/atcc_hifi_classifier_sweep.py"
for REQUIRED_FILE in \
    "${MANIFEST}" \
    "${TRUTH_MANIFEST}" \
    "${EXISTING_KRAKEN2_REPORT}" \
    "${EXISTING_METABULI_REPORT}" \
    "${SWEEP_SCRIPT}"
do
    [[ -s "${REQUIRED_FILE}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${REQUIRED_FILE}" >&2
        exit 2
    }
done

conda run --no-capture-output --name "${CONDA_ENVIRONMENT}" \
    python "${SWEEP_SCRIPT}" \
    --action summarise \
    --manifest "${MANIFEST}" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --sweep-root "${SWEEP_ROOT}" \
    --existing-kraken2-report "${EXISTING_KRAKEN2_REPORT}" \
    --existing-metabuli-report "${EXISTING_METABULI_REPORT}" \
    --output-directory "${OUTPUT_DIRECTORY}" \
    --focus-taxid 99158 \
    --verbose
