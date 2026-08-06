#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Initialise or run a new Nanopore dataset without changing package files." \
        "" \
        "Required:" \
        "  --action initialise|validate|dry-run|run" \
        "  --config PATH" \
        "" \
        "Options:" \
        "  --samples PATH   Sample-sheet path used by initialise" \
        "                   (default: CONFIG basename with .samples.tsv)" \
        "  --profile PATH   Snakemake profile for dry-run/run" \
        "  --jobs INT       Maximum jobs (default: 10)" \
        "  --help" \
        "" \
        "Examples:" \
        "  bash scripts/start_new_dataset.sh --action initialise --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action validate --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action dry-run --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action run --config config/my_run.yaml"
}

ACTION=""
CONFIG_PATH=""
SAMPLES_PATH=""
PROFILE_PATH=""
JOBS="10"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --action)
            [[ $# -ge 2 ]] || { printf 'ERROR: --action requires a value\n' >&2; exit 2; }
            ACTION="$2"
            shift 2
            ;;
        --config)
            [[ $# -ge 2 ]] || { printf 'ERROR: --config requires a value\n' >&2; exit 2; }
            CONFIG_PATH="$2"
            shift 2
            ;;
        --samples)
            [[ $# -ge 2 ]] || { printf 'ERROR: --samples requires a value\n' >&2; exit 2; }
            SAMPLES_PATH="$2"
            shift 2
            ;;
        --profile)
            [[ $# -ge 2 ]] || { printf 'ERROR: --profile requires a value\n' >&2; exit 2; }
            PROFILE_PATH="$2"
            shift 2
            ;;
        --jobs)
            [[ $# -ge 2 ]] || { printf 'ERROR: --jobs requires a value\n' >&2; exit 2; }
            JOBS="$2"
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

[[ -n "${ACTION}" ]] || { printf 'ERROR: --action is required\n' >&2; exit 2; }
[[ -n "${CONFIG_PATH}" ]] || { printf 'ERROR: --config is required\n' >&2; exit 2; }
[[ "${JOBS}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: --jobs must be a positive integer\n' >&2
    exit 2
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
if [[ -z "${SAMPLES_PATH}" ]]; then
    if [[ "${CONFIG_PATH}" == *.yaml ]]; then
        SAMPLES_PATH="${CONFIG_PATH%.yaml}.samples.tsv"
    elif [[ "${CONFIG_PATH}" == *.yml ]]; then
        SAMPLES_PATH="${CONFIG_PATH%.yml}.samples.tsv"
    else
        SAMPLES_PATH="${CONFIG_PATH}.samples.tsv"
    fi
fi

case "${ACTION}" in
    initialise)
        if [[ -e "${CONFIG_PATH}" || -e "${SAMPLES_PATH}" ]]; then
            printf 'ERROR: refusing to overwrite an existing config or sample sheet:\n' >&2
            printf '  %s\n  %s\n' "${CONFIG_PATH}" "${SAMPLES_PATH}" >&2
            exit 3
        fi
        mkdir -p -- "$(dirname -- "${CONFIG_PATH}")" "$(dirname -- "${SAMPLES_PATH}")"
        cp -- "${REPO_DIR}/config/real_data.template.yaml" "${CONFIG_PATH}"
        cp -- "${REPO_DIR}/config/samples.template.tsv" "${SAMPLES_PATH}"
        printf 'Created configuration: %s\n' "${CONFIG_PATH}"
        printf 'Created sample sheet: %s\n' "${SAMPLES_PATH}"
        printf '%s\n' \
            "Next: add one FASTQ path per row to the TSV, edit every absolute resource path" \
            "in the YAML, set inputs.read_state deliberately, then run --action validate."
        ;;
    validate)
        [[ -f "${CONFIG_PATH}" ]] || { printf 'ERROR: missing config: %s\n' "${CONFIG_PATH}" >&2; exit 2; }
        nanopore-realdata --action validate --config "${CONFIG_PATH}"
        ;;
    dry-run|run)
        [[ -f "${CONFIG_PATH}" ]] || { printf 'ERROR: missing config: %s\n' "${CONFIG_PATH}" >&2; exit 2; }
        COMMAND=(
            bash "${REPO_DIR}/scripts/run_workflow.sh"
            --config "${CONFIG_PATH}"
            --jobs "${JOBS}"
            --verbose
        )
        if [[ -n "${PROFILE_PATH}" ]]; then
            COMMAND+=(--profile "${PROFILE_PATH}")
        fi
        if [[ "${ACTION}" == "dry-run" ]]; then
            COMMAND+=(--dry-run)
        fi
        "${COMMAND[@]}"
        ;;
    *)
        printf 'ERROR: unsupported --action: %s\n' "${ACTION}" >&2
        usage >&2
        exit 2
        ;;
esac
