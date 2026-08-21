# ATCC HiFi Metabuli long-read recovery v0.4.6

## Cause

The first Metabuli sweep supplied one PacBio HiFi FASTQ without selecting
Metabuli long-read sequence mode. Metabuli requires `--seq-mode 3` for this
input layout. The `--precise 2` option is a separate HiFi classification
preset and does not select the input layout.

## Changes

- Add `--seq-mode 3` to every generated Metabuli classification command.
- Test long-read mode for explicit score thresholds and the HiFi preset.
- Validate `--seq-mode` support during preflight.
- Add a recovery launcher that submits only the five failed Metabuli tasks.
- Make the replacement summary wait for the existing Kraken2 array and the
  corrected Metabuli array.

Failed diagnostic directories are retained. Existing Kraken2 results,
Metabuli baseline results, databases and reads are not overwritten.

## Mac installation

Extract this archive, preview the overlay with `rsync -avn`, then apply it to
the Mac repository with `rsync -av`. Commit and push from the Mac repository,
then pull on the cluster.

## Cluster recovery

Cancel only the unsatisfied original summary job:

```bash
scancel 257602
```

Validate without submitting:

```bash
bash scripts/recover_atcc_hifi_metabuli_sweep.sh \
    --mode plan \
    --kraken2-job-id 257600
```

Submit the five corrected Metabuli tasks and replacement summary:

```bash
bash scripts/recover_atcc_hifi_metabuli_sweep.sh \
    --mode submit \
    --kraken2-job-id 257600
```

## Verification

The overlay passes 17 unit tests, Python byte-code compilation and Bash syntax
validation.
