#!/usr/bin/env bash

set -Eeuo pipefail

MODE=""

usage() {
    printf '%s\n' \
        'Usage: submit_atcc_hifi_matched_kraken2.sh --mode plan|submit' \
        '' \
        'plan   checks paths and prints the intended three-job chain' \
        'submit checks paths and submits database, classification and aggregation jobs'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="${2:?Missing value for --mode}"
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

[[ "${MODE}" == plan || "${MODE}" == submit ]] || {
    printf 'ERROR: --mode must be plan or submit\n' >&2
    exit 2
}

REPOSITORY_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow"
CONFIG_PATH="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_srr9328980_20260821_matched_kraken2.yaml"
TRUTH_MANIFEST="${REPOSITORY_ROOT}/config/atcc_msa1003_hifi_srr9328980_20260821.kraken_reference_truth.tsv"
RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/benchmarks/comparators_atcc_msa1003_hifi/atcc_msa1003_hifi_comparators_srr9328980_20260817"
REFERENCE_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/databases/kmersutra_db/atcc_msa1003_heldout_v1"
TAXONOMY_ROOT="$(dirname "${REFERENCE_ROOT}")/ncbi_taxonomy"
DATABASE_ROOT="${REFERENCE_ROOT}/kraken2_matched_475genomes_v1"
CONDA_ENVIRONMENT="nanopore_realdata_workflow"
LOG_DIRECTORY="${RUN_ROOT}/workflow_control/slurm_logs"
ARCHIVE_SUFFIX="superseded_unmatched_database_20260821"
BUILD_SCRIPT="${REPOSITORY_ROOT}/scripts/build_atcc_matched_kraken2_database.sh"
CLASSIFY_SCRIPT="${REPOSITORY_ROOT}/scripts/run_atcc_matched_kraken2_classification.sh"
AGGREGATE_SCRIPT="${REPOSITORY_ROOT}/scripts/run_atcc_matched_kraken2_aggregation.sh"

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || {
        printf 'ERROR: required file is missing or empty: %s\n' "${path}" >&2
        exit 2
    }
}

require_directory() {
    local path="$1"
    [[ -d "${path}" ]] || {
        printf 'ERROR: required directory is missing: %s\n' "${path}" >&2
        exit 2
    }
}

require_directory "${REPOSITORY_ROOT}"
require_directory "${RUN_ROOT}"
require_directory "${REFERENCE_ROOT}"
require_directory "${TAXONOMY_ROOT}"
require_file "${CONFIG_PATH}"
require_file "${TRUTH_MANIFEST}"
require_file "${REFERENCE_ROOT}/kmersutra_genome_config.tsv"
require_file "${REFERENCE_ROOT}/reference_audit/atcc_reference_gate_summary.tsv"
require_file "${TAXONOMY_ROOT}/names.dmp"
require_file "${TAXONOMY_ROOT}/nodes.dmp"
require_file "${BUILD_SCRIPT}"
require_file "${CLASSIFY_SCRIPT}"
require_file "${AGGREGATE_SCRIPT}"
require_file "${RUN_ROOT}/02_classification/metabuli/SRR9328980/complete.json"

if [[ -e "${RUN_ROOT}/02_classification/kraken2/SRR9328980.${ARCHIVE_SUFFIX}" ]]; then
    printf 'ERROR: the superseded-result archive already exists; refusing a duplicate rerun: %s\n' \
        "${RUN_ROOT}/02_classification/kraken2/SRR9328980.${ARCHIVE_SUFFIX}" >&2
    exit 2
fi
if [[ -e "${RUN_ROOT}/03_final.${ARCHIVE_SUFFIX}" ]]; then
    printf 'ERROR: the superseded aggregate archive already exists; refusing a duplicate rerun: %s\n' \
        "${RUN_ROOT}/03_final.${ARCHIVE_SUFFIX}" >&2
    exit 2
fi

mkdir -p "${LOG_DIRECTORY}"

printf 'Matched Kraken2 recovery plan:\n'
printf '  source genomes: %s\n' "${REFERENCE_ROOT}/kmersutra_genome_config.tsv"
printf '  matched database: %s\n' "${DATABASE_ROOT}"
printf '  truth manifest: %s\n' "${TRUTH_MANIFEST}"
printf '  reads: SRR9328980 (existing sample manifest)\n'
printf '  preserved comparator: Metabuli\n'
printf '  rerun comparator: Kraken2 only\n'
printf '  final stage: corrected aggregation\n'

if [[ "${MODE}" == plan ]]; then
    printf 'Plan checks passed. No jobs were submitted.\n'
    exit 0
fi

command -v sbatch >/dev/null 2>&1 || {
    printf 'ERROR: sbatch is unavailable\n' >&2
    exit 2
}

BUILD_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_k2db475 \
        --chdir="${REPOSITORY_ROOT}" \
        --cpus-per-task=24 \
        --mem=256000M \
        --time=2880 \
        --signal=B:TERM@300 \
        --output="${LOG_DIRECTORY}/%x.%j.out" \
        --error="${LOG_DIRECTORY}/%x.%j.err" \
        "${BUILD_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --reference-root "${REFERENCE_ROOT}" \
        --taxonomy-root "${TAXONOMY_ROOT}" \
        --truth-manifest "${TRUTH_MANIFEST}" \
        --database-root "${DATABASE_ROOT}" \
        --threads 24
)"

CLASSIFY_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_k2fair \
        --chdir="${REPOSITORY_ROOT}" \
        --cpus-per-task=12 \
        --mem=409600M \
        --time=2880 \
        --signal=B:TERM@300 \
        --dependency="afterok:${BUILD_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%j.out" \
        --error="${LOG_DIRECTORY}/%x.%j.err" \
        "${CLASSIFY_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --config "${CONFIG_PATH}" \
        --run-root "${RUN_ROOT}" \
        --archive-suffix "${ARCHIVE_SUFFIX}"
)"

AGGREGATE_JOB_ID="$(
    sbatch --parsable \
        --account=barton \
        --partition=barton \
        --job-name=NRD_k2agg \
        --chdir="${REPOSITORY_ROOT}" \
        --cpus-per-task=2 \
        --mem=8192M \
        --time=180 \
        --signal=B:TERM@300 \
        --dependency="afterok:${CLASSIFY_JOB_ID}" \
        --output="${LOG_DIRECTORY}/%x.%j.out" \
        --error="${LOG_DIRECTORY}/%x.%j.err" \
        "${AGGREGATE_SCRIPT}" \
        --repository-root "${REPOSITORY_ROOT}" \
        --conda-environment "${CONDA_ENVIRONMENT}" \
        --config "${CONFIG_PATH}" \
        --run-root "${RUN_ROOT}" \
        --archive-suffix "${ARCHIVE_SUFFIX}"
)"

printf 'Submitted matched Kraken2 correction chain:\n'
printf '  database:       %s\n' "${BUILD_JOB_ID}"
printf '  classification: %s\n' "${CLASSIFY_JOB_ID}"
printf '  aggregation:    %s\n' "${AGGREGATE_JOB_ID}"
printf 'Metabuli and KmerSutra were not resubmitted.\n'
