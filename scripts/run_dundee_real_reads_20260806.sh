#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Run the Dundee host-removed real-read classifier benchmark." \
        "" \
        "Options:" \
        "  --dry-run   Validate and display the Snakemake DAG only" \
        "  --unlock    Remove a stale Snakemake lock before resuming" \
        "  --help"
}

DRY_RUN=false
UNLOCK=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --unlock)
            UNLOCK=true
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

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${REPO_DIR}/config/dundee_real_reads_20260806.yaml"
PROFILE_PATH="${REPO_DIR}/workflow/profiles/dundee"

[[ -f "${CONFIG_PATH}" ]] || {
    printf 'ERROR: workflow configuration is missing: %s\n' "${CONFIG_PATH}" >&2
    exit 2
}

COMMAND=(
    bash "${REPO_DIR}/scripts/run_workflow.sh"
    --config "${CONFIG_PATH}"
    --profile "${PROFILE_PATH}"
    --jobs 10
    --verbose
)

if [[ "${DRY_RUN}" == true ]]; then
    COMMAND+=(--dry-run)
fi
if [[ "${UNLOCK}" == true ]]; then
    COMMAND+=(--unlock)
fi

printf 'Repository: %s\n' "${REPO_DIR}"
printf 'Configuration: %s\n' "${CONFIG_PATH}"
printf 'Starting at: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
"${COMMAND[@]}"
