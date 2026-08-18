# ATCC MSA-1003 HiFi truth and minimap2 recovery

This patch corrects two issues in the first comparator run:

- truth concordance was restricted to Plasmodium taxa and could not evaluate
  the 20-organism bacterial mock community; and
- the minimap2 reference requirement used three legacy product names rather
  than the current taxonomy names in the controlled genome configuration.

Kraken2 and Metabuli completed successfully and are not resubmitted. The
recovery builds the controlled minimap2 reference and index, runs minimap2 for
SRR9328980, then regenerates the final reports using the reference-composition
truth table.

## Apply and run

1. Extract the supplied archive beside the HPC repository.
2. Copy the overlay into the repository:

   ```bash
   REPO=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow

   rsync -av nanopore_realdata_truth_fix_v0_4_3/ "${REPO}/"
   ```

3. Activate the existing environment and refresh the editable installation:

   ```bash
   conda activate nanopore_realdata_workflow
   cd "${REPO}"
   python -m pip install --no-deps -e .
   python -m pip show nanopore-realdata-workflow | awk '/^Version:/ {print $2}'
   ```

   The printed version must be `0.4.3`.

4. Run the focused tests:

   ```bash
   PYTHONPATH=src python -m unittest -v \
       tests.test_pcr \
       tests.test_atcc_hifi_truth
   ```

5. Submit only the minimap2 recovery chain and final aggregation:

   ```bash
   bash scripts/recover_atcc_hifi_truth_and_minimap2.sh
   ```

6. Monitor the four new jobs:

   ```bash
   squeue -u pthorpe001
   ```

7. When the aggregation job has finished, check:

   ```bash
   RUN=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/benchmarks/comparators_atcc_msa1003_hifi/atcc_msa1003_hifi_comparators_srr9328980_20260817

   column -t -s $'\t' "${RUN}/03_final/classifier_status.tsv"
   column -t -s $'\t' "${RUN}/03_final/pcr_method_summary.tsv"
   sed -n '1,12p' "${RUN}/03_final/pcr_concordance.tsv"
   ```

The recovery script makes one non-destructive copy of the original final
report at `03_final.before_truth_fix_20260818` before regenerating it.
