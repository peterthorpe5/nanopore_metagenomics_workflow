#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Usage:" \
        "  bash scripts/run_workflow.sh --config PATH [options]" \
        "" \
        "Options:" \
        "  --profile PATH   Snakemake profile; omit for local execution" \
        "  --cores INT|all  Available cores (default: all)" \
        "  --jobs INT       Maximum concurrent jobs (default: 10)" \
        "  --target NAME    Optional Snakemake target; repeatable" \
        "  --set-resource RULE:RESOURCE=VALUE  Repeatable Snakemake override" \
        "  --dry-run" \
        "  --unlock" \
        "  --verbose" \
        "  --help"
}

CONFIG_PATH=""
PROFILE_PATH=""
CORES="all"
JOBS="10"
DRY_RUN=false
UNLOCK=false
VERBOSE=false
TARGETS=()
RESOURCE_OVERRIDES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            [[ $# -ge 2 ]] || { printf 'ERROR: --config requires a value\n' >&2; exit 2; }
            CONFIG_PATH="$2"
            shift 2
            ;;
        --profile)
            [[ $# -ge 2 ]] || { printf 'ERROR: --profile requires a value\n' >&2; exit 2; }
            PROFILE_PATH="$2"
            shift 2
            ;;
        --cores)
            [[ $# -ge 2 ]] || { printf 'ERROR: --cores requires a value\n' >&2; exit 2; }
            CORES="$2"
            shift 2
            ;;
        --jobs)
            [[ $# -ge 2 ]] || { printf 'ERROR: --jobs requires a value\n' >&2; exit 2; }
            JOBS="$2"
            shift 2
            ;;
        --target)
            [[ $# -ge 2 ]] || { printf 'ERROR: --target requires a value\n' >&2; exit 2; }
            TARGETS+=(--target "$2")
            shift 2
            ;;
        --set-resource)
            [[ $# -ge 2 ]] || { printf 'ERROR: --set-resource requires a value\n' >&2; exit 2; }
            RESOURCE_OVERRIDES+=(--set-resource "$2")
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --unlock)
            UNLOCK=true
            shift
            ;;
        --verbose)
            VERBOSE=true
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

[[ -n "${CONFIG_PATH}" ]] || { printf 'ERROR: --config is required\n' >&2; exit 2; }
[[ -f "${CONFIG_PATH}" ]] || { printf 'ERROR: config does not exist: %s\n' "${CONFIG_PATH}" >&2; exit 2; }

COMMAND=(
    nanopore-realdata
    --action run
    --config "${CONFIG_PATH}"
    --cores "${CORES}"
    --jobs "${JOBS}"
)
if [[ -n "${PROFILE_PATH}" ]]; then
    COMMAND+=(--profile "${PROFILE_PATH}")
fi
if [[ "${DRY_RUN}" == true ]]; then
    COMMAND+=(--dry-run)
fi
if [[ "${UNLOCK}" == true ]]; then
    COMMAND+=(--unlock)
fi
if [[ "${VERBOSE}" == true ]]; then
    COMMAND+=(--verbose)
fi
COMMAND+=("${TARGETS[@]}")
COMMAND+=("${RESOURCE_OVERRIDES[@]}")

printf 'Starting real-data workflow at %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
set +e
"${COMMAND[@]}"
WORKFLOW_EXIT_CODE=$?
set -e

if [[ "${WORKFLOW_EXIT_CODE}" -ne 0 \
    && "${DRY_RUN}" == false \
    && "${UNLOCK}" == false ]]; then
    printf 'Workflow exited with code %s; attempting a truthful partial report.\n' \
        "${WORKFLOW_EXIT_CODE}" >&2
    if nanopore-realdata --action report --config "${CONFIG_PATH}" --verbose; then
        printf 'Partial report generated; the original workflow failure remains unresolved.\n' >&2
    else
        printf 'Partial reporting was not possible, usually because core input preparation failed.\n' >&2
    fi
fi

exit "${WORKFLOW_EXIT_CODE}"
