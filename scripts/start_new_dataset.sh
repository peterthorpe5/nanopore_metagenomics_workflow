#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

usage() {
    printf '%s\n' \
        "Initialise or run a new Nanopore dataset without changing package files." \
        "" \
        "Required:" \
        "  --action initialise|validate|plan-slurm|submit-slurm|" \
        "           resume-submission|new-attempt|dry-run|run" \
        "  --config PATH" \
        "" \
        "Options:" \
        "  --samples PATH   Sample-sheet path used by initialise" \
        "                   (default: CONFIG basename with .samples.tsv)" \
        "  --with-pcr-truth Create an optional PCR truth template" \
        "  --pcr-truth PATH PCR truth path used by initialise; implies" \
        "                   --with-pcr-truth" \
        "  --profile PATH   Snakemake profile for dry-run/run" \
        "  --jobs INT       Maximum jobs (default: 10)" \
        "  --retry-method NAME  With plan-slurm or submit-slurm, retry only" \
        "                       this classifier and then aggregate; repeatable" \
        "  --help" \
        "" \
        "Examples:" \
        "  bash scripts/start_new_dataset.sh --action initialise --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action validate --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action plan-slurm --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action submit-slurm --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action dry-run --config config/my_run.yaml" \
        "  bash scripts/start_new_dataset.sh --action run --config config/my_run.yaml"
}

ACTION=""
CONFIG_PATH=""
SAMPLES_PATH=""
PCR_TRUTH_PATH=""
WITH_PCR_TRUTH="false"
PROFILE_PATH=""
JOBS="10"
RETRY_METHODS=()

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
        --pcr-truth)
            [[ $# -ge 2 ]] || { printf 'ERROR: --pcr-truth requires a value\n' >&2; exit 2; }
            PCR_TRUTH_PATH="$2"
            WITH_PCR_TRUTH="true"
            shift 2
            ;;
        --with-pcr-truth)
            WITH_PCR_TRUTH="true"
            shift
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

[[ -n "${ACTION}" ]] || { printf 'ERROR: --action is required\n' >&2; exit 2; }
[[ -n "${CONFIG_PATH}" ]] || { printf 'ERROR: --config is required\n' >&2; exit 2; }
[[ "${JOBS}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: --jobs must be a positive integer\n' >&2
    exit 2
}
if [[ ${#RETRY_METHODS[@]} -gt 0 \
    && "${ACTION}" != "plan-slurm" \
    && "${ACTION}" != "submit-slurm" ]]; then
    printf 'ERROR: --retry-method requires --action plan-slurm or submit-slurm\n' >&2
    exit 2
fi

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
if [[ "${WITH_PCR_TRUTH}" == "true" && -z "${PCR_TRUTH_PATH}" ]]; then
    if [[ "${CONFIG_PATH}" == *.yaml ]]; then
        PCR_TRUTH_PATH="${CONFIG_PATH%.yaml}.pcr_truth.tsv"
    elif [[ "${CONFIG_PATH}" == *.yml ]]; then
        PCR_TRUTH_PATH="${CONFIG_PATH%.yml}.pcr_truth.tsv"
    else
        PCR_TRUTH_PATH="${CONFIG_PATH}.pcr_truth.tsv"
    fi
fi

case "${ACTION}" in
    initialise)
        if [[ -e "${CONFIG_PATH}" || -e "${SAMPLES_PATH}" \
            || ( "${WITH_PCR_TRUTH}" == "true" && -e "${PCR_TRUTH_PATH}" ) ]]; then
            printf 'ERROR: refusing to overwrite an existing run template:\n' >&2
            printf '  %s\n  %s\n' "${CONFIG_PATH}" "${SAMPLES_PATH}" >&2
            if [[ "${WITH_PCR_TRUTH}" == "true" ]]; then
                printf '  %s\n' "${PCR_TRUTH_PATH}" >&2
            fi
            exit 3
        fi
        mkdir -p -- \
            "$(dirname -- "${CONFIG_PATH}")" \
            "$(dirname -- "${SAMPLES_PATH}")"
        if [[ "${WITH_PCR_TRUTH}" == "true" ]]; then
            mkdir -p -- "$(dirname -- "${PCR_TRUTH_PATH}")"
        fi
        cp -- "${REPO_DIR}/config/real_data.template.yaml" "${CONFIG_PATH}"
        cp -- "${REPO_DIR}/config/samples.template.tsv" "${SAMPLES_PATH}"
        printf 'Created configuration: %s\n' "${CONFIG_PATH}"
        printf 'Created sample sheet: %s\n' "${SAMPLES_PATH}"
        if [[ "${WITH_PCR_TRUTH}" == "true" ]]; then
            cp -- "${REPO_DIR}/config/pcr_truth.template.tsv" "${PCR_TRUTH_PATH}"
            printf 'Created optional PCR truth table: %s\n' "${PCR_TRUTH_PATH}"
            printf '%s\n' \
                "Set inputs.pcr_truth in the YAML to this file after populating it."
        fi
        printf '%s\n' \
            "Next: populate the sample TSV, edit every absolute resource path in the YAML," \
            "set inputs.read_state deliberately, then run --action validate."
        ;;
    validate)
        [[ -f "${CONFIG_PATH}" ]] || { printf 'ERROR: missing config: %s\n' "${CONFIG_PATH}" >&2; exit 2; }
        nanopore-realdata --action validate --config "${CONFIG_PATH}"
        ;;
    plan-slurm|submit-slurm|resume-submission|new-attempt)
        [[ -f "${CONFIG_PATH}" ]] || { printf 'ERROR: missing config: %s\n' "${CONFIG_PATH}" >&2; exit 2; }
        CLI_ACTION="${ACTION}"
        SUBMISSION_FLAG=""
        if [[ "${ACTION}" == "resume-submission" ]]; then
            CLI_ACTION="submit-slurm"
            SUBMISSION_FLAG="--resume-submission"
        elif [[ "${ACTION}" == "new-attempt" ]]; then
            CLI_ACTION="submit-slurm"
            SUBMISSION_FLAG="--new-attempt"
        fi
        COMMAND=(
            nanopore-realdata
            --action "${CLI_ACTION}"
            --config "${CONFIG_PATH}"
            --verbose
        )
        if [[ -n "${SUBMISSION_FLAG}" ]]; then
            COMMAND+=("${SUBMISSION_FLAG}")
        fi
        for RETRY_METHOD in "${RETRY_METHODS[@]}"; do
            COMMAND+=(--retry-method "${RETRY_METHOD}")
        done
        "${COMMAND[@]}"
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
