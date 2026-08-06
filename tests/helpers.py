"""Synthetic input helpers for real-data workflow tests."""

from __future__ import annotations

import gzip
from pathlib import Path

import yaml


def write_fastq(*, path: Path, read_id: str = "read_1") -> None:
    """Write one valid synthetic FASTQ record."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"@{read_id}\nACGT\n+\nIIII\n"
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def build_test_project(
    *,
    root: Path,
    kmersutra_enabled: bool = True,
    failure_policy: str = "continue",
    stage_resources: bool = False,
    input_read_state: str = "raw",
) -> Path:
    """Create a minimal valid project and return its YAML path."""
    root.mkdir(parents=True, exist_ok=True)
    fastq = root / "reads" / "sample.fastq.gz"
    write_fastq(path=fastq)
    samples = root / "samples.tsv"
    samples.write_text(
        "sample_id\tfastq\trun_id\tbarcode\tdescription\n"
        f"sample_1\t{fastq}\trun_1\tbarcode01\tSynthetic\n",
        encoding="utf-8",
    )
    host_reference = root / "host.fasta"
    host_reference.write_text(">host\nACGTACGT\n", encoding="utf-8")
    minimap_reference = root / "classification.fasta"
    minimap_reference.write_text(
        ">kraken:taxid|1|reference_a Species alpha\nACGTACGT\n",
        encoding="utf-8",
    )
    kraken = root / "kraken"
    kraken.mkdir()
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (kraken / name).write_text(name, encoding="utf-8")
    metabuli = root / "metabuli"
    metabuli.mkdir()
    (metabuli / "db").write_text("database", encoding="utf-8")
    panel = root / "species_kmer_panel.tsv.gz"
    with gzip.open(panel, "wt", encoding="utf-8") as handle:
        handle.write("kmer\tk\tspecies_name\nAAAA\t51\tSpecies alpha\n")
    config = {
        "schema_version": 2,
        "run": {
            "id": "test_run",
            "output_directory": str(root / "results"),
        },
        "inputs": {"samples": str(samples), "read_state": input_read_state},
        "host": {"reference": str(host_reference), "index": ""},
        "databases": {
            "kraken2": str(kraken),
            "metabuli": str(metabuli),
            "kmersutra_panel": str(panel),
        },
        "minimap2": {
            "reference": str(minimap_reference),
            "index": "",
            "min_mapq": 15,
            "min_alignment": 500,
        },
        "execution": {
            "scratch_root": str(root),
            "stage_resources": stage_resources,
            "minimum_scratch_gb": 1,
            "keep_non_host_fastq": True,
            "keep_per_read_classifications": True,
        },
        "resources": {
            "host_threads": 2,
            "kraken2_threads": 3,
            "metabuli_threads": 4,
            "minimap2_threads": 4,
            "kmersutra_threads": 5,
            "host_memory_mb": 1000,
            "kraken2_memory_mb": 2000,
            "metabuli_memory_mb": 3000,
            "minimap2_memory_mb": 3000,
            "kmersutra_memory_mb": 4000,
            "host_runtime_minutes": 10,
            "kraken2_runtime_minutes": 20,
            "metabuli_runtime_minutes": 30,
            "minimap2_runtime_minutes": 30,
            "kmersutra_runtime_minutes": 40,
        },
        "kraken2": {"confidence": 0.0},
        "metabuli": {"min_score": 0.008, "max_ram_gb": 120},
        "kmersutra": {
            "enabled": kmersutra_enabled,
            "failure_policy": failure_policy,
            "timeout_minutes": 35,
            "screen_preset": "exact",
            "call_preset": "lineage_aware",
            "same_genus_reportable_min_fraction": 0.05,
            "write_parquet_outputs": False,
        },
        "provenance": {"checksum_inputs": False},
    }
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path
