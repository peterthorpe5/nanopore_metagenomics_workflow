#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Validate and submit Kraken2, Metabuli and minimap2 for the ATCC MSA-1003 HiFi run." \
        "" \
        "Usage:" \
        "  bash scripts/submit_atcc_msa1003_hifi_comparators_20260817.sh [options]" \
        "" \
        "Options:" \
        "  --config PATH   Override the supplied benchmark configuration." \
        "  --plan-only     Validate and print the Slurm plan without submitting." \
        "  --help          Show this help." \
        "" \
        "KmerSutra is disabled in this comparator-only configuration. Raw per-read" \
        "outputs are retained for direct investigation of Hammondia hammondi evidence."
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_srr9328980_20260817.yaml"
ACTION="submit-slurm"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            [[ $# -ge 2 ]] || {
                printf 'ERROR: --config requires a path\n' >&2
                exit 2
            }
            CONFIG_PATH="$2"
            shift 2
            ;;
        --plan-only)
            ACTION="plan-slurm"
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

if [[ ! -s "${CONFIG_PATH}" ]]; then
    printf 'ERROR: comparator configuration is missing or empty: %s\n' "${CONFIG_PATH}" >&2
    exit 1
fi
if ! command -v nanopore-realdata >/dev/null 2>&1; then
    printf '%s\n' \
        'ERROR: nanopore-realdata is unavailable in the active environment.' \
        'Activate the nanopore_realdata_workflow conda environment and retry.' >&2
    exit 1
fi

printf 'Validating comparator configuration: %s\n' "${CONFIG_PATH}"
bash "${REPOSITORY_ROOT}/scripts/start_new_dataset.sh" \
    --action validate \
    --config "${CONFIG_PATH}"

printf 'Running action: %s\n' "${ACTION}"
bash "${REPOSITORY_ROOT}/scripts/start_new_dataset.sh" \
    --action "${ACTION}" \
    --config "${CONFIG_PATH}"
