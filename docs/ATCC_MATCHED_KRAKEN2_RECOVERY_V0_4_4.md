# ATCC MSA-1003 matched Kraken2 correction

## Purpose

The earlier Kraken2 result for SRR9328980 used a database that did not contain the 20 bacterial members of the ATCC MSA-1003 community. Its apparent 0/20 recovery is therefore superseded and must not be used in the manuscript.

This correction builds a Kraken2 database from the same held-out 475-genome source collection used to construct the KmerSutra reference. This makes the source sequence collection comparable while leaving the algorithms and their reporting behaviour unchanged.

## Locked reference design

The builder requires all of the following before it invokes `kraken2-build`:

- exactly 475 source genome assemblies;
- exactly 396 source species;
- all 20 expected ATCC species represented, allowing the locked accepted-name mappings;
- an explicit PASS in the existing KmerSutra reference gate;
- no truth assembly accession in the source manifest;
- no truth sequence accession in any source FASTA header;
- valid positive NCBI taxids for every source assembly; and
- local NCBI `names.dmp` and `nodes.dmp` taxonomy files.

Every source FASTA header is rewritten in Kraken2's explicit taxid form. The database is then built with Kraken2's standard k-mer and minimizer settings. No additional genomes are downloaded and the published truth sequences are not added.

The final database is accepted only if `hash.k2d`, `opts.k2d` and `taxo.k2d` are non-empty and `kraken2-inspect` contains the taxid of every expected species. The database retains the source assembly manifest, build summary, inspect report and completion record.

## Recovery behaviour

The submission script creates three dependent jobs:

1. build and validate `kraken2_matched_475genomes_v1`;
2. archive the old unmatched Kraken2 result and classify SRR9328980 with the matched database; and
3. archive the old aggregate and regenerate the comparison tables.

The dependencies use `afterok`. Classification cannot start after an invalid database build, and aggregation cannot start after a failed classification. Metabuli and KmerSutra are not resubmitted.

The superseded files are moved, not deleted:

- `02_classification/kraken2/SRR9328980.superseded_unmatched_database_20260821`
- `03_final.superseded_unmatched_database_20260821`

## Deployment

After the overlay has been copied into the Mac repository, commit and push it. On the cluster, pull the commit, reinstall the editable package and run the targeted tests:

```bash
cd /gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow
git pull
conda activate nanopore_realdata_workflow
python -m pip install --no-deps -e .
python -m unittest -v tests.test_matched_kraken2_reference
bash -n scripts/build_atcc_matched_kraken2_database.sh
bash -n scripts/run_atcc_matched_kraken2_classification.sh
bash -n scripts/run_atcc_matched_kraken2_aggregation.sh
bash -n scripts/submit_atcc_hifi_matched_kraken2.sh
```

Check the plan without submitting jobs:

```bash
bash scripts/submit_atcc_hifi_matched_kraken2.sh --mode plan
```

Submit exactly one correction chain:

```bash
bash scripts/submit_atcc_hifi_matched_kraken2.sh --mode submit
```

Do not interpret the new Kraken2 result until the aggregate job has completed and the following files exist:

```text
kraken2_matched_475genomes_v1/matched_database.complete.json
02_classification/kraken2/SRR9328980/complete.json
03_final/workflow.complete.json
03_final/pcr_method_summary.tsv
```
