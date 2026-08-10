"""Configuration and sample-manifest validation for the real-data workflow."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml


SCHEMA_VERSION = 3
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


@dataclass(frozen=True)
class Sample:
    """One logical sample and its ordered FASTQ parts."""

    sample_id: str
    fastq_paths: tuple[Path, ...]
    run_id: str = ""
    barcode: str = ""
    description: str = ""


@dataclass(frozen=True)
class WorkflowConfig:
    """Validated configuration needed by Snakemake and execution adapters."""

    config_path: Path
    run_id: str
    output_directory: Path
    samples_path: Path
    pcr_truth_path: Path
    input_read_state: str
    expected_repository_root: Path
    expected_package_version: str
    conda_environment: str
    host_reference: Path | None
    host_index: Path | None
    kraken_database: Path
    metabuli_database: Path
    kmersutra_panel: Path | None
    minimap_reference: Path
    minimap_genome_config: Path | None
    minimap_maximum_reference_bases: int
    minimap_index_batch_size_bases: int
    minimap_maximum_index_bytes: int
    minimap_index: Path | None
    scratch_root: Path
    stage_resources: bool
    minimum_scratch_gb: int
    keep_non_host_fastq: bool
    keep_per_read_classifications: bool
    threads_host: int
    threads_kraken2: int
    threads_metabuli: int
    threads_minimap2: int
    threads_kmersutra: int
    memory_host_mb: int
    memory_kraken2_mb: int
    memory_metabuli_mb: int
    memory_minimap2_mb: int
    memory_kmersutra_mb: int
    runtime_host_minutes: int
    runtime_kraken2_minutes: int
    runtime_metabuli_minutes: int
    runtime_minimap2_minutes: int
    runtime_kmersutra_minutes: int
    kraken_confidence: float
    kraken_failure_policy: str
    metabuli_min_score: float
    metabuli_max_ram_gb: int
    metabuli_failure_policy: str
    minimap_min_mapq: int
    minimap_min_alignment: int
    minimap_failure_policy: str
    kmersutra_screen_preset: str
    kmersutra_call_preset: str
    kmersutra_same_genus_fraction: float
    kmersutra_write_parquet: bool
    kmersutra_enabled: bool
    kmersutra_failure_policy: str
    kmersutra_timeout_minutes: int
    report_focus_taxa: tuple[str, ...]
    report_top_n: int
    report_max_table_rows: int
    checksum_inputs: bool
    slurm_account: str
    slurm_partition: str
    slurm_default_qos: str
    slurm_kmersutra_qos: str
    slurm_kraken2_concurrency: int
    slurm_metabuli_concurrency: int
    slurm_minimap2_concurrency: int
    slurm_kmersutra_concurrency: int
    samples: tuple[Sample, ...]


def load_workflow_config(*, config_path: Path) -> WorkflowConfig:
    """Load and validate a real-data workflow YAML file.

    Args:
        config_path: YAML configuration file.

    Returns:
        A fully resolved and validated workflow configuration.

    Raises:
        ValueError: If configuration fields are missing or invalid.
        FileNotFoundError: If a required configured path does not exist.
    """
    resolved_config = config_path.expanduser().resolve()
    if not resolved_config.is_file():
        raise FileNotFoundError(f"Configuration file does not exist: {resolved_config}")
    raw = yaml.safe_load(resolved_config.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Workflow configuration must contain a YAML mapping")
    if raw.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"schema_version must be {SCHEMA_VERSION}")

    base = resolved_config.parent
    run = _mapping(raw, "run")
    inputs = _mapping(raw, "inputs")
    deployment = _mapping(raw, "deployment")
    host = _mapping(raw, "host")
    databases = _mapping(raw, "databases")
    execution = _mapping(raw, "execution")
    resources = _mapping(raw, "resources")
    kraken2 = _mapping(raw, "kraken2")
    metabuli = _mapping(raw, "metabuli")
    kmersutra = _mapping(raw, "kmersutra")
    provenance = _mapping(raw, "provenance")
    slurm = _mapping(raw, "slurm")
    slurm_concurrency = _mapping(slurm, "array_concurrency")
    reporting = _optional_mapping(raw, "reporting")
    kmersutra_enabled = _boolean(kmersutra, "enabled")
    minimap2 = _mapping(raw, "minimap2")

    run_id = _identifier(run, "id")
    output_directory = _path(run, "output_directory", base=base, must_exist=False)
    samples_path = _path(inputs, "samples", base=base, must_exist=True)
    pcr_truth_path = _path(inputs, "pcr_truth", base=base, must_exist=True)
    input_read_state = _choice(
        inputs,
        "read_state",
        choices={"raw", "host_removed"},
    )
    host_reference = _optional_path(host, "reference", base=base, must_exist=True)
    host_index = _optional_path(host, "index", base=base, must_exist=True)
    expected_repository_root = _path(
        deployment,
        "expected_repository_root",
        base=base,
        must_exist=False,
    )
    expected_package_version = _non_empty_string(
        deployment,
        "expected_package_version",
    )
    conda_environment = _identifier(deployment, "conda_environment")
    # Classifier resources are deliberately resolved without requiring them to
    # exist at configuration-load time. Their readiness is assessed per branch
    # during preflight so one unavailable classifier cannot suppress all other
    # results or the final partial-results report.
    kraken_database = _path(databases, "kraken2", base=base, must_exist=False)
    metabuli_database = _path(databases, "metabuli", base=base, must_exist=False)
    kmersutra_panel = _optional_path(
        databases,
        "kmersutra_panel",
        base=base,
        must_exist=False,
    )
    minimap_reference_source = _optional_path(
        minimap2,
        "reference",
        base=base,
        must_exist=False,
    )
    minimap_genome_config = _optional_path(
        minimap2,
        "genome_config",
        base=base,
        must_exist=True,
    )
    if (minimap_reference_source is None) == (minimap_genome_config is None):
        raise ValueError("Configure exactly one of minimap2.reference or minimap2.genome_config")
    minimap_reference = (
        minimap_reference_source
        if minimap_reference_source is not None
        else output_directory / "00_preflight" / "controlled_minimap_reference.fa"
    )
    minimap_index = _optional_path(
        minimap2,
        "index",
        base=base,
        must_exist=False,
    )
    if minimap_index is not None:
        raise ValueError(
            "Prebuilt minimap2.index values are not accepted in schema v3; "
            "the workflow must build and checksum a single-part index"
        )
    scratch_root = _path(execution, "scratch_root", base=base, must_exist=True)

    if kraken_database.exists() and not kraken_database.is_dir():
        raise ValueError(f"Kraken2 database must be a directory: {kraken_database}")
    if metabuli_database.exists() and not metabuli_database.is_dir():
        raise ValueError(f"Metabuli database must be a directory: {metabuli_database}")
    if kmersutra_enabled and kmersutra_panel is None:
        raise ValueError("kmersutra_panel is required when KmerSutra is enabled")
    if kmersutra_panel is not None and kmersutra_panel.exists() and not kmersutra_panel.is_file():
        raise ValueError(f"KmerSutra panel must be a file: {kmersutra_panel}")
    if input_read_state == "raw" and host_reference is None:
        raise ValueError("host.reference is required when inputs.read_state is raw")
    if host_index is not None and not host_index.is_file():
        raise ValueError(f"Host minimap2 index must be a file: {host_index}")
    if host_reference is not None and not host_reference.is_file():
        raise ValueError(f"Host reference must be a file: {host_reference}")
    if minimap_reference.exists() and not minimap_reference.is_file():
        raise ValueError(f"Classification minimap2 reference must be a file: {minimap_reference}")
    if not scratch_root.is_dir():
        raise ValueError(f"Scratch root must be a directory: {scratch_root}")
    maximum_reference_bases = _positive_int(minimap2, "maximum_reference_bases")
    index_batch_size_bases = _positive_int(minimap2, "index_batch_size_bases")
    if index_batch_size_bases < maximum_reference_bases:
        raise ValueError(
            "minimap2.index_batch_size_bases must be at least as large as "
            "minimap2.maximum_reference_bases to guarantee a single index part"
        )

    samples = load_samples(samples_path=samples_path)
    workflow = WorkflowConfig(
        config_path=resolved_config,
        run_id=run_id,
        output_directory=output_directory,
        samples_path=samples_path,
        pcr_truth_path=pcr_truth_path,
        input_read_state=input_read_state,
        expected_repository_root=expected_repository_root,
        expected_package_version=expected_package_version,
        conda_environment=conda_environment,
        host_reference=host_reference,
        host_index=host_index,
        kraken_database=kraken_database,
        metabuli_database=metabuli_database,
        kmersutra_panel=kmersutra_panel,
        minimap_reference=minimap_reference,
        minimap_genome_config=minimap_genome_config,
        minimap_maximum_reference_bases=maximum_reference_bases,
        minimap_index_batch_size_bases=index_batch_size_bases,
        minimap_maximum_index_bytes=_positive_int(
            minimap2,
            "maximum_index_bytes",
        ),
        minimap_index=minimap_index,
        scratch_root=scratch_root,
        stage_resources=_boolean(execution, "stage_resources"),
        minimum_scratch_gb=_positive_int(execution, "minimum_scratch_gb"),
        keep_non_host_fastq=_boolean(execution, "keep_non_host_fastq"),
        keep_per_read_classifications=_boolean(
            execution,
            "keep_per_read_classifications",
        ),
        threads_host=_positive_int(resources, "host_threads"),
        threads_kraken2=_positive_int(resources, "kraken2_threads"),
        threads_metabuli=_positive_int(resources, "metabuli_threads"),
        threads_minimap2=_positive_int(resources, "minimap2_threads"),
        threads_kmersutra=_positive_int(resources, "kmersutra_threads"),
        memory_host_mb=_positive_int(resources, "host_memory_mb"),
        memory_kraken2_mb=_positive_int(resources, "kraken2_memory_mb"),
        memory_metabuli_mb=_positive_int(resources, "metabuli_memory_mb"),
        memory_minimap2_mb=_positive_int(resources, "minimap2_memory_mb"),
        memory_kmersutra_mb=_positive_int(resources, "kmersutra_memory_mb"),
        runtime_host_minutes=_positive_int(resources, "host_runtime_minutes"),
        runtime_kraken2_minutes=_positive_int(resources, "kraken2_runtime_minutes"),
        runtime_metabuli_minutes=_positive_int(resources, "metabuli_runtime_minutes"),
        runtime_minimap2_minutes=_positive_int(resources, "minimap2_runtime_minutes"),
        runtime_kmersutra_minutes=_positive_int(resources, "kmersutra_runtime_minutes"),
        kraken_confidence=_bounded_float(kraken2, "confidence", minimum=0.0, maximum=1.0),
        kraken_failure_policy=_choice_default(
            kraken2,
            "failure_policy",
            choices={"continue", "fail"},
            default="continue",
        ),
        metabuli_min_score=_bounded_float(
            metabuli,
            "min_score",
            minimum=0.0,
            maximum=1.0,
        ),
        metabuli_max_ram_gb=_positive_int(metabuli, "max_ram_gb"),
        metabuli_failure_policy=_choice_default(
            metabuli,
            "failure_policy",
            choices={"continue", "fail"},
            default="continue",
        ),
        minimap_min_mapq=_non_negative_int(minimap2, "min_mapq"),
        minimap_min_alignment=_positive_int(minimap2, "min_alignment"),
        minimap_failure_policy=_choice_default(
            minimap2,
            "failure_policy",
            choices={"continue", "fail"},
            default="continue",
        ),
        kmersutra_screen_preset=_choice(
            kmersutra,
            "screen_preset",
            choices={"exact", "raw_ont_sensitive"},
        ),
        kmersutra_call_preset=_choice(
            kmersutra,
            "call_preset",
            choices={"conservative", "lineage_aware", "raw_ont_sensitive", "strict"},
        ),
        kmersutra_same_genus_fraction=_bounded_float(
            kmersutra,
            "same_genus_reportable_min_fraction",
            minimum=0.0,
            maximum=1.0,
        ),
        kmersutra_write_parquet=_boolean(kmersutra, "write_parquet_outputs"),
        kmersutra_enabled=kmersutra_enabled,
        kmersutra_failure_policy=_choice(
            kmersutra,
            "failure_policy",
            choices={"continue", "fail"},
        ),
        kmersutra_timeout_minutes=_positive_int(kmersutra, "timeout_minutes"),
        report_focus_taxa=_string_tuple_default(
            reporting,
            "focus_taxa",
            default=("Plasmodium",),
        ),
        report_top_n=_positive_int_default(reporting, "top_n", default=20),
        report_max_table_rows=_positive_int_default(
            reporting,
            "max_table_rows",
            default=5000,
        ),
        checksum_inputs=_boolean(provenance, "checksum_inputs"),
        slurm_account=_identifier(slurm, "account"),
        slurm_partition=_identifier(slurm, "partition"),
        slurm_default_qos=_optional_identifier(slurm, "default_qos"),
        slurm_kmersutra_qos=_identifier(slurm, "kmersutra_qos"),
        slurm_kraken2_concurrency=_positive_int(slurm_concurrency, "kraken2"),
        slurm_metabuli_concurrency=_positive_int(slurm_concurrency, "metabuli"),
        slurm_minimap2_concurrency=_positive_int(slurm_concurrency, "minimap2"),
        slurm_kmersutra_concurrency=_positive_int(slurm_concurrency, "kmersutra"),
        samples=samples,
    )
    _validate_workflow_paths(workflow=workflow)
    return workflow


def load_samples(*, samples_path: Path) -> tuple[Sample, ...]:
    """Load a TSV sample sheet and group repeated sample identifiers.

    Args:
        samples_path: TSV with required ``sample_id`` and ``fastq`` columns.

    Returns:
        Samples in first-observed order, with FASTQ parts in row order.

    Raises:
        ValueError: If rows, identifiers or paths are invalid.
    """
    resolved = samples_path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Sample sheet does not exist: {resolved}")
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "fastq"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("Sample sheet requires sample_id and fastq columns")
        grouped: dict[str, dict[str, Any]] = {}
        seen_paths: set[Path] = set()
        for row_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not SAMPLE_ID_PATTERN.fullmatch(sample_id):
                raise ValueError(
                    f"Invalid sample_id on row {row_number}: {sample_id!r}; "
                    "use letters, numbers, dots, underscores or hyphens"
                )
            fastq_text = (row.get("fastq") or "").strip()
            if not fastq_text:
                raise ValueError(f"Blank fastq path on row {row_number}")
            fastq_path = _resolve_path(fastq_text, base=resolved.parent)
            if not any(str(fastq_path).lower().endswith(suffix) for suffix in FASTQ_SUFFIXES):
                raise ValueError(f"Unsupported FASTQ suffix on row {row_number}: {fastq_path}")
            if not fastq_path.is_file() or fastq_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FASTQ is missing or empty on row {row_number}: {fastq_path}"
                )
            if fastq_path in seen_paths:
                raise ValueError(f"FASTQ is listed more than once: {fastq_path}")
            seen_paths.add(fastq_path)
            metadata = {
                "run_id": (row.get("run_id") or "").strip(),
                "barcode": (row.get("barcode") or "").strip(),
                "description": (row.get("description") or "").strip(),
            }
            if sample_id not in grouped:
                grouped[sample_id] = {"paths": [], **metadata}
            else:
                for field, value in metadata.items():
                    current = grouped[sample_id][field]
                    if value and current and value != current:
                        raise ValueError(f"Conflicting {field} values for sample {sample_id!r}")
                    if value and not current:
                        grouped[sample_id][field] = value
            grouped[sample_id]["paths"].append(fastq_path)
    if not grouped:
        raise ValueError("Sample sheet contains no data rows")
    return tuple(
        Sample(
            sample_id=sample_id,
            fastq_paths=tuple(values["paths"]),
            run_id=values["run_id"],
            barcode=values["barcode"],
            description=values["description"],
        )
        for sample_id, values in grouped.items()
    )


def _validate_workflow_paths(*, workflow: WorkflowConfig) -> None:
    """Reject configurations that could overwrite source data or databases."""
    protected = {
        workflow.samples_path,
        workflow.pcr_truth_path,
        workflow.kraken_database,
        workflow.metabuli_database,
        *(sample.fastq_paths for sample in workflow.samples),
    }
    if workflow.minimap_genome_config is None:
        protected.add(workflow.minimap_reference)
    else:
        protected.add(workflow.minimap_genome_config)
    if workflow.host_reference is not None:
        protected.add(workflow.host_reference)
    if workflow.host_index is not None:
        protected.add(workflow.host_index)
    if workflow.minimap_index is not None:
        protected.add(workflow.minimap_index)
    if workflow.kmersutra_panel is not None:
        protected.add(workflow.kmersutra_panel)
    flattened: set[Path] = set()
    for item in protected:
        if isinstance(item, tuple):
            flattened.update(item)
        else:
            flattened.add(item)
    output = workflow.output_directory
    for path in flattened:
        if output == path or output in path.parents:
            raise ValueError(
                f"Output directory must not contain a configured input or database: {path}"
            )


def _mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Configuration section {key!r} must be a mapping")
    return value


def _optional_mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    """Return an optional configuration mapping.

    Args:
        parent: Parent configuration mapping.
        key: Section name.

    Returns:
        The configured mapping or an empty mapping when omitted.

    Raises:
        ValueError: If the optional section exists but is not a mapping.
    """
    value = parent.get(key, {})
    if not isinstance(value, dict):
        raise ValueError(f"Configuration section {key!r} must be a mapping")
    return value


def _identifier(mapping: Mapping[str, Any], key: str) -> str:
    value = str(mapping.get(key, "")).strip()
    if not SAMPLE_ID_PATTERN.fullmatch(value):
        raise ValueError(f"{key} must be a filesystem-safe identifier")
    return value


def _non_empty_string(mapping: Mapping[str, Any], key: str) -> str:
    """Return a required non-empty configuration string."""
    value = str(mapping.get(key, "")).strip()
    if not value:
        raise ValueError(f"{key} must be a non-empty string")
    return value


def _optional_string(mapping: Mapping[str, Any], key: str) -> str:
    """Return a stripped optional configuration string."""
    return str(mapping.get(key, "")).strip()


def _optional_identifier(mapping: Mapping[str, Any], key: str) -> str:
    """Return an optional scheduler-safe identifier."""
    value = _optional_string(mapping, key)
    if value and not SAMPLE_ID_PATTERN.fullmatch(value):
        raise ValueError(f"{key} must be a scheduler-safe identifier")
    return value


def _path(
    mapping: Mapping[str, Any],
    key: str,
    *,
    base: Path,
    must_exist: bool,
) -> Path:
    value = str(mapping.get(key, "")).strip()
    if not value:
        raise ValueError(f"Configuration path {key!r} is required")
    path = _resolve_path(value, base=base)
    if must_exist and not path.exists():
        raise FileNotFoundError(f"Configured path does not exist for {key}: {path}")
    return path


def _optional_path(
    mapping: Mapping[str, Any],
    key: str,
    *,
    base: Path,
    must_exist: bool,
) -> Path | None:
    value = str(mapping.get(key, "")).strip()
    if not value:
        return None
    path = _resolve_path(value, base=base)
    if must_exist and not path.exists():
        raise FileNotFoundError(f"Configured path does not exist for {key}: {path}")
    return path


def _resolve_path(value: str, *, base: Path) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(value))
    unresolved = re.findall(r"\$\{[^}]+\}|\$[A-Za-z_][A-Za-z0-9_]*", expanded)
    if unresolved:
        raise ValueError(f"Unresolved environment variable in path: {value}")
    path = Path(expanded)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _boolean(mapping: Mapping[str, Any], key: str) -> bool:
    value = mapping.get(key)
    if not isinstance(value, bool):
        raise ValueError(f"{key} must be true or false")
    return value


def _positive_int(mapping: Mapping[str, Any], key: str) -> int:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{key} must be a positive integer")
    return value


def _positive_int_default(
    mapping: Mapping[str, Any],
    key: str,
    *,
    default: int,
) -> int:
    """Return an optional positive integer with a validated default."""
    if key not in mapping:
        return default
    return _positive_int(mapping, key)


def _non_negative_int(mapping: Mapping[str, Any], key: str) -> int:
    """Return a validated integer greater than or equal to zero."""
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{key} must be a non-negative integer")
    return value


def _bounded_float(
    mapping: Mapping[str, Any],
    key: str,
    *,
    minimum: float,
    maximum: float,
) -> float:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{key} must be numeric")
    numeric = float(value)
    if not minimum <= numeric <= maximum:
        raise ValueError(f"{key} must be between {minimum} and {maximum}")
    return numeric


def _choice(mapping: Mapping[str, Any], key: str, *, choices: set[str]) -> str:
    value = str(mapping.get(key, "")).strip()
    if value not in choices:
        raise ValueError(f"{key} must be one of: {', '.join(sorted(choices))}")
    return value


def _choice_default(
    mapping: Mapping[str, Any],
    key: str,
    *,
    choices: set[str],
    default: str,
) -> str:
    """Return an optional constrained string with a validated default."""
    if key not in mapping:
        return default
    return _choice(mapping, key, choices=choices)


def _string_tuple_default(
    mapping: Mapping[str, Any],
    key: str,
    *,
    default: tuple[str, ...],
) -> tuple[str, ...]:
    """Return a non-empty, de-duplicated tuple of report focus terms."""
    if key not in mapping:
        return default
    value = mapping.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a YAML list of non-empty strings")
    cleaned: list[str] = []
    seen: set[str] = set()
    for item in value:
        text = str(item).strip()
        if not text:
            raise ValueError(f"{key} must not contain blank values")
        folded = text.casefold()
        if folded not in seen:
            seen.add(folded)
            cleaned.append(text)
    if not cleaned:
        raise ValueError(f"{key} must contain at least one value")
    return tuple(cleaned)
