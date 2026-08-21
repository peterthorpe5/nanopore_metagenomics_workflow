# ATCC MSA-1003 HiFi classifier operating-point sweep

## Purpose

This patch replaces single-setting interpretation of Kraken2 and Metabuli with a
predeclared operating-point analysis. It reuses the existing matched-database
Kraken2 confidence-0.00 result and existing Metabuli min-score-0.008 result.
It does not rebuild either database and does not overwrite the baseline outputs.

## Locked settings

Kraken2 is run at confidence 0.05, 0.10, 0.20 and 0.50. Together with the
existing documented default of 0.00, this gives five confidence settings.

Metabuli is evaluated at:

- its documented maximum-sensitivity default, min-score 0 and min-sp-score 0;
- the existing min-score 0.008 and min-sp-score 0 result, which corresponds to
  the published Oxford Nanopore score recommendation;
- the documented PacBio HiFi min-score of 0.07;
- the explicit PacBio HiFi thresholds, min-score 0.07 and min-sp-score 0.3;
- the current Metabuli v1.2 HiFi precise preset, `--precise 2`;
- a stricter sensitivity-analysis point, min-score 0.10 and min-sp-score 0.3.

The settings are locked in
`config/atcc_msa1003_hifi_classifier_operating_points_20260821.tsv` before the
new results are observed.

## Outputs

New classifier outputs are written below:

```text
04_classifier_operating_point_sweep/
├── kraken2/<setting_id>/
├── metabuli/<setting_id>/
└── summary/
```

Each newly run setting contains its report, gzip-compressed per-read
classifications, stdout and stderr logs, metadata and a completion record. The
summary contains:

- `classifier_operating_point_summary.tsv`;
- `classifier_species_direct_read_counts.tsv`;
- the locked operating-point manifest;
- the locked 20-species truth manifest;
- a README and completion record.

Expected and additional species are counted at direct species rank using support
thresholds of at least 1, 2, 10 and 100 reads. *Hammondia hammondi* taxid 99158
is reported separately. The counts describe a known-composition mock community;
they are not estimates of clinical specificity.

## Submission

Run the non-mutating preflight first:

```bash
bash scripts/submit_atcc_hifi_classifier_sweep.sh --mode plan
```

Submit only if the plan passes:

```bash
bash scripts/submit_atcc_hifi_classifier_sweep.sh --mode submit
```

The launcher submits independent Kraken2 and Metabuli arrays, followed by one
summary job with an `afterok` dependency on both arrays.

## Primary documentation

- Kraken2 confidence scoring:
  <https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring>
- Metabuli v1.2 classify options and HiFi precise preset:
  <https://jaebeom-kim.github.io/metabuli-doc/modules/classify/>
