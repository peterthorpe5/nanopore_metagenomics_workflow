#!/usr/bin/env bash

set -Eeuo pipefail

REPOSITORY_ROOT=""
CONDA_ENVIRONMENT=""
REFERENCE_ROOT=""
TAXONOMY_ROOT=""
TRUTH_MANIFEST=""
DATABASE_ROOT=""
THREADS=""

usage() {
    printf '%s\n' \
        'Usage: build_atcc_matched_kraken2_database.sh [named options]' \
        '  --repository-root PATH' \
        '  --conda-environment NAME' \
        '  --reference-root PATH' \
        '  --taxonomy-root PATH' \
        '  --truth-manifest PATH' \
        '  --database-root PATH' \
        '  --threads INTEGER'
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
        --reference-root)
            REFERENCE_ROOT="${2:?Missing value for --reference-root}"
            shift 2
            ;;
        --taxonomy-root)
            TAXONOMY_ROOT="${2:?Missing value for --taxonomy-root}"
            shift 2
            ;;
        --truth-manifest)
            TRUTH_MANIFEST="${2:?Missing value for --truth-manifest}"
            shift 2
            ;;
        --database-root)
            DATABASE_ROOT="${2:?Missing value for --database-root}"
            shift 2
            ;;
        --threads)
            THREADS="${2:?Missing value for --threads}"
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

require_directory() {
    local path="$1"
    [[ -d "${path}" ]] || {
        printf 'ERROR: required directory is missing: %s\n' "${path}" >&2
        exit 2
    }
}

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${path}" >&2
        exit 2
    }
}

require_directory "${REPOSITORY_ROOT}"
require_directory "${REFERENCE_ROOT}"
require_directory "${TAXONOMY_ROOT}"
[[ -n "${CONDA_ENVIRONMENT}" ]] || {
    printf 'ERROR: --conda-environment is required\n' >&2
    exit 2
}
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: --threads must be a positive integer\n' >&2
    exit 2
}
[[ "${DATABASE_ROOT}" = /* ]] || {
    printf 'ERROR: --database-root must be an absolute path\n' >&2
    exit 2
}

REFERENCE_TOOL="${REPOSITORY_ROOT}/scripts/matched_kraken2_reference.py"
GENOME_CONFIG="${REFERENCE_ROOT}/kmersutra_genome_config.tsv"
GATE_SUMMARY="${REFERENCE_ROOT}/reference_audit/atcc_reference_gate_summary.tsv"
require_file "${REFERENCE_TOOL}"
require_file "${GENOME_CONFIG}"
require_file "${GATE_SUMMARY}"
require_file "${TRUTH_MANIFEST}"
require_file "${TAXONOMY_ROOT}/names.dmp"
require_file "${TAXONOMY_ROOT}/nodes.dmp"

CONDA_COMMAND=(
    conda run
    --no-capture-output
    --name "${CONDA_ENVIRONMENT}"
)

if [[ -d "${DATABASE_ROOT}" ]]; then
    require_file "${DATABASE_ROOT}/matched_database.complete.json"
    require_file "${DATABASE_ROOT}/provenance/kraken2_database_inspect.tsv"
    require_file "${DATABASE_ROOT}/provenance/library_preparation.json"
    "${CONDA_COMMAND[@]}" python "${REFERENCE_TOOL}" \
        --action validate-database \
        --database "${DATABASE_ROOT}" \
        --inspect-report "${DATABASE_ROOT}/provenance/kraken2_database_inspect.tsv" \
        --truth-manifest "${TRUTH_MANIFEST}" \
        --preparation-summary "${DATABASE_ROOT}/provenance/library_preparation.json" \
        --verbose
    printf 'Matched Kraken2 database already exists and passed validation: %s\n' \
        "${DATABASE_ROOT}"
    exit 0
fi

DATABASE_PARENT="$(dirname "${DATABASE_ROOT}")"
require_directory "${DATABASE_PARENT}"
BUILD_ROOT="${DATABASE_ROOT}.building.${SLURM_JOB_ID:-$$}"
[[ ! -e "${BUILD_ROOT}" ]] || {
    printf 'ERROR: build workspace already exists: %s\n' "${BUILD_ROOT}" >&2
    exit 2
}

on_error() {
    local exit_code="$?"
    trap - ERR
    printf 'ERROR: matched Kraken2 database build failed with exit code %s\n' \
        "${exit_code}" >&2
    printf 'The incomplete build was retained for diagnosis: %s\n' "${BUILD_ROOT}" >&2
    exit "${exit_code}"
}
trap on_error ERR

mkdir -p "${BUILD_ROOT}/taxonomy" "${BUILD_ROOT}/provenance"
cp -p "${TAXONOMY_ROOT}/names.dmp" "${BUILD_ROOT}/taxonomy/names.dmp"
cp -p "${TAXONOMY_ROOT}/nodes.dmp" "${BUILD_ROOT}/taxonomy/nodes.dmp"
for optional_taxonomy_file in merged.dmp delnodes.dmp; do
    if [[ -s "${TAXONOMY_ROOT}/${optional_taxonomy_file}" ]]; then
        cp -p \
            "${TAXONOMY_ROOT}/${optional_taxonomy_file}" \
            "${BUILD_ROOT}/taxonomy/${optional_taxonomy_file}"
    fi
done

COMBINED_FASTA="${BUILD_ROOT}/provenance/matched_475_genomes.fna"
ASSEMBLY_MANIFEST="${BUILD_ROOT}/provenance/matched_assembly_manifest.tsv"
PREPARATION_SUMMARY="${BUILD_ROOT}/provenance/library_preparation.json"

printf 'Started matched Kraken2 build UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'Host: %s\n' "$(hostname)"
printf 'Source reference: %s\n' "${REFERENCE_ROOT}"
printf 'Build workspace: %s\n' "${BUILD_ROOT}"

"${CONDA_COMMAND[@]}" python "${REFERENCE_TOOL}" \
    --action prepare-library \
    --genome-config "${GENOME_CONFIG}" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --gate-summary "${GATE_SUMMARY}" \
    --output-fasta "${COMBINED_FASTA}" \
    --output-manifest "${ASSEMBLY_MANIFEST}" \
    --output-summary "${PREPARATION_SUMMARY}" \
    --expected-genome-count 475 \
    --expected-source-species-count 396 \
    --expected-truth-species-count 20 \
    --verbose

"${CONDA_COMMAND[@]}" kraken2-build \
    --add-to-library "${COMBINED_FASTA}" \
    --db "${BUILD_ROOT}"

"${CONDA_COMMAND[@]}" kraken2-build \
    --build \
    --threads "${THREADS}" \
    --db "${BUILD_ROOT}"

INSPECT_REPORT="${BUILD_ROOT}/provenance/kraken2_database_inspect.tsv"
"${CONDA_COMMAND[@]}" kraken2-inspect --db "${BUILD_ROOT}" > "${INSPECT_REPORT}"

"${CONDA_COMMAND[@]}" python "${REFERENCE_TOOL}" \
    --action validate-database \
    --database "${BUILD_ROOT}" \
    --inspect-report "${INSPECT_REPORT}" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --preparation-summary "${PREPARATION_SUMMARY}" \
    --verbose

mv "${BUILD_ROOT}" "${DATABASE_ROOT}"
BUILD_ROOT="${DATABASE_ROOT}"

"${CONDA_COMMAND[@]}" python "${REFERENCE_TOOL}" \
    --action validate-database \
    --database "${DATABASE_ROOT}" \
    --inspect-report "${DATABASE_ROOT}/provenance/kraken2_database_inspect.tsv" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --preparation-summary "${DATABASE_ROOT}/provenance/library_preparation.json" \
    --output-summary "${DATABASE_ROOT}/matched_database.complete.json" \
    --verbose

trap - ERR
printf 'Completed matched Kraken2 build UTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'Validated database: %s\n' "${DATABASE_ROOT}"
