#!/usr/bin/env bash

set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "${LINENO}" "${BASH_COMMAND}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_ROOT="${REPO_ROOT}/test_results"
RUN_LABEL="$(date -u +%Y%m%dT%H%M%SZ)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir)
            RESULTS_ROOT="${2:?Missing value for --results-dir}"
            shift 2
            ;;
        --run-label)
            RUN_LABEL="${2:?Missing value for --run-label}"
            shift 2
            ;;
        --help|-h)
            printf 'Usage: bash scripts/run_quality_checks.sh [--results-dir PATH] [--run-label LABEL]\n'
            exit 0
            ;;
        *)
            printf 'ERROR: unknown option: %s\n' "$1" >&2
            exit 2
            ;;
    esac
done

RUN_DIR="${RESULTS_ROOT%/}/${RUN_LABEL}"
mkdir -p "${RUN_DIR}/coverage_html"
cd "${REPO_ROOT}"

python -m ruff check src tests 2>&1 | tee "${RUN_DIR}/ruff.log"
python -m ruff format --check src tests 2>&1 | tee "${RUN_DIR}/ruff_format.log"
python -m coverage erase
python -m coverage run -m unittest discover -s tests -p 'test_*.py' \
    2>&1 | tee "${RUN_DIR}/unittest.log"
python -m coverage report 2>&1 | tee "${RUN_DIR}/coverage.txt"
python -m coverage html -d "${RUN_DIR}/coverage_html"
python -m compileall -q src tests
python -m build --outdir "${RUN_DIR}/dist" . \
    2>&1 | tee "${RUN_DIR}/package_build.log"
while IFS= read -r script_path; do
    bash -n "${script_path}"
done < <(find scripts -maxdepth 1 -type f -name '*.sh' -print | sort)
GIT_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ "${GIT_ROOT}" == "${REPO_ROOT}" ]]; then
    git diff --check 2>&1 | tee "${RUN_DIR}/git_diff_check.log"
else
    printf 'Not a Git worktree; skipped git diff --check.\n' \
        | tee "${RUN_DIR}/git_diff_check.log"
fi

printf 'Quality results: %s\n' "${RUN_DIR}"
