ATCC HiFi Kraken2 and Metabuli operating-point sweep

This overlay adds a non-destructive classifier sensitivity analysis. It does not
replace the matched databases or the existing confidence-0.00 Kraken2 and
min-score-0.008 Metabuli results.

Mac repository target:
/Users/PThorpe001/github_repos/nanopore_metagenomics_workflow

After applying the overlay, run the unit tests from the repository root:

PYTHONPATH=scripts python -m unittest -v tests.test_atcc_hifi_classifier_sweep

On the cluster, always run plan mode before submission:

bash scripts/submit_atcc_hifi_classifier_sweep.sh --mode plan
bash scripts/submit_atcc_hifi_classifier_sweep.sh --mode submit

See docs/ATCC_HIFI_CLASSIFIER_OPERATING_POINT_SWEEP_V0_4_5.md for the complete
locked grid, outputs and interpretation rules.
