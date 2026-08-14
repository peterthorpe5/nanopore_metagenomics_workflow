# v0.4.2 selective-retry recovery

Release date: 13 August 2026

## Scope

This recovery applies to the prepared run at:

```text
/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/runs/real_reads_classifier_benchmark_v2_20260810
```

The first attempt produced valid results for all 11 Metabuli, minimap2 and
KmerSutra tasks. All 11 Kraken2 tasks exceeded their 96,000 MB allocations.
Final aggregation then stopped on an integer-formatting error in the optional
PCR comparison layer.

Version 0.4.2 fixes the aggregation edge case and adds a generic selective
classifier retry. The Dundee benchmark preset assigns Kraken2 409,600 MB
(400 GiB). This does not change the reusable template's default allocation.

PCR remains optional and affects evaluation only. A routine dataset requires a
sample manifest and configured classifier resources, not a PCR truth table.

## Deployment route

Apply and validate the v0.4.2 patch in the Mac Git checkout. Commit and push it
to `peterthorpe5/nanopore_metagenomics_workflow`. Do not copy source files
directly to the cluster.

On Dundee, fast-forward the canonical checkout and reinstall it:

```bash
set -Eeuo pipefail

PROJECT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity"
REPOSITORY="${PROJECT}/nanopore_metagenomics_workflow"

cd "${REPOSITORY}"
git status --short
git pull --ff-only

conda run --name nanopore_realdata_workflow \
    python -m pip install --no-deps --editable "${REPOSITORY}"

conda run --name nanopore_realdata_workflow \
    python -c 'import nanopore_realdata; print(nanopore_realdata.__version__)'
```

The printed version must be `0.4.2` and `git status --short` must be empty.

## Safe plan

Do not use `--new-attempt`; it represents a complete DAG retry. Plan only the
failed method:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh \
    --plan \
    --retry-method kraken2
```

The plan must contain exactly:

1. `retry_classify_kraken2`: array `0-10%1`, 12 CPUs, 409,600 MB; and
2. `retry_aggregate`: dependent on the Kraken2 retry array with `afterany`.

It must not contain Metabuli, minimap2, KmerSutra, preflight, input acceptance
or index construction.

## Submit

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh \
    --retry-method kraken2
```

The submitter verifies that all 11 classification-ready FASTQs still exist,
checks that the previous recorded jobs are no longer active, archives the old
journal and writes a new selective-retry journal. Existing successful
classifier directories are not modified by submission.

If this submission is interrupted between its two `sbatch` calls, resume the
stored selective plan with:

```bash
bash scripts/submit_dundee_real_reads_v2_20260810.sh --resume-submission
```

## Final checks

After the retry aggregation job ends:

```bash
OUTPUT="${PROJECT}/runs/real_reads_classifier_benchmark_v2_20260810"
FINAL="${OUTPUT}/03_final"

python -m json.tool "${FINAL}/workflow.complete.json"

awk -F $'\t' '
    NR > 1 {counts[$2 "\t" $3]++}
    END {
        print "method\tstatus\tcount"
        for (key in counts) {
            print key "\t" counts[key]
        }
    }
' "${FINAL}/classifier_status.tsv" | sort | column -t -s $'\t'
```

The intended final state is `status: success` with 11 successful results for
each of Kraken2, Metabuli, minimap2 and KmerSutra. If Kraken2 fails again,
aggregation must still complete and report a truthful partial result.
