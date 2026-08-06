"""Execution adapters and reporting for the real Nanopore Snakemake workflow."""

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from nanopore_realdata import __version__
from nanopore_realdata.commands import (
    kmersutra_command,
    kraken2_command,
    metabuli_command,
    minimap2_classification_command,
    minimap2_host_command,
    minimap2_index_command,
    pigz_command,
    required_executables,
    samtools_bam_command,
    samtools_count_command,
    samtools_non_host_fastq_command,
)
from nanopore_realdata.config import Sample, WorkflowConfig, load_workflow_config
from nanopore_realdata.minimap import summarise_minimap_paf
from nanopore_realdata.runtime import (
    capture_output,
    completion_is_valid,
    metadata_fingerprint,
    publish_directory,
    run_command,
    run_pipeline,
    scratch_workspace,
    sha256_file,
    stage_resource,
    task_signature,
    utc_now,
    validate_scratch,
    write_completion,
    write_json_atomic,
)


LOGGER = logging.getLogger(__name__)


def preflight(*, config_path: Path, output_path: Path) -> None:
    """Validate configuration, input records, databases and executables.

    Args:
        config_path: Workflow YAML path.
        output_path: Declared preflight JSON output.
    """
    workflow = load_workflow_config(config_path=config_path)
    actions = ["classify-kraken2", "classify-metabuli", "classify-minimap2"]
    if workflow.input_read_state == "raw":
        actions.extend(("build-host-index", "host-deplete"))
    else:
        actions.append("accept-host-removed")
    if workflow.kmersutra_enabled:
        actions.append("classify-kmersutra")
    for action in actions:
        _require_tools(action=action)
    _validate_kraken_database(database=workflow.kraken_database)
    _validate_non_empty_directory(label="Metabuli database", path=workflow.metabuli_database)
    if workflow.kmersutra_enabled:
        panel = _required_kmersutra_panel(workflow=workflow)
        _validate_panel(panel=panel)
    sample_rows = []
    for sample in workflow.samples:
        for part_number, fastq_path in enumerate(sample.fastq_paths, start=1):
            _validate_fastq_prefix(path=fastq_path)
            sample_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "part_number": part_number,
                    "fastq": str(fastq_path),
                    "size_bytes": fastq_path.stat().st_size,
                    "run_id": sample.run_id,
                    "barcode": sample.barcode,
                    "description": sample.description,
                }
            )

    tools = {
        "minimap2": _version(command=["minimap2", "--version"]),
        "samtools": _version(command=["samtools", "--version"]),
        "kraken2": _version(command=["kraken2", "--version"]),
        "metabuli": _version(command=["metabuli", "version"]),
        "kmersutra": (_kmersutra_version() if workflow.kmersutra_enabled else "disabled"),
        "nanopore_realdata_workflow": __version__,
        "snakemake": _version(command=["snakemake", "--version"]),
        "rsync": _version(command=["rsync", "--version"]),
        "pigz": _version(command=["pigz", "--version"]),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        path=output_path.parent / "resolved_samples.tsv",
        rows=sample_rows,
        fieldnames=(
            "sample_id",
            "part_number",
            "fastq",
            "size_bytes",
            "run_id",
            "barcode",
            "description",
        ),
    )
    _write_tsv(
        path=output_path.parent / "software_versions.tsv",
        rows=[{"software": name, "version": value} for name, value in tools.items()],
        fieldnames=("software", "version"),
    )
    payload = {
        "status": "success",
        "checked_at_utc": utc_now(),
        "run_id": workflow.run_id,
        "sample_count": len(workflow.samples),
        "fastq_part_count": len(sample_rows),
        "config": str(workflow.config_path),
        "config_sha256": sha256_file(path=workflow.config_path),
        "resources": {
            "input_read_state": workflow.input_read_state,
            "host_reference": (
                metadata_fingerprint(
                    paths=[_required_host_reference(workflow=workflow)],
                    checksum_files=True,
                )
                if workflow.input_read_state == "raw"
                else "not_applicable_already_host_removed"
            ),
            "kraken2_database": metadata_fingerprint(
                paths=[workflow.kraken_database],
                checksum_files=False,
            ),
            "metabuli_database": metadata_fingerprint(
                paths=[workflow.metabuli_database],
                checksum_files=False,
            ),
            "minimap2_reference": metadata_fingerprint(
                paths=[workflow.minimap_reference],
                checksum_files=True,
            ),
            "minimap2_index": (
                metadata_fingerprint(
                    paths=[workflow.minimap_index],
                    checksum_files=True,
                )
                if workflow.minimap_index is not None
                else "built_by_workflow_from_configured_reference"
            ),
            "kmersutra_panel": (
                metadata_fingerprint(
                    paths=[_required_kmersutra_panel(workflow=workflow)],
                    checksum_files=True,
                )
                if workflow.kmersutra_enabled
                else "disabled"
            ),
        },
        "software": tools,
    }
    write_json_atomic(path=output_path, payload=payload)


def build_host_index(
    *,
    config_path: Path,
    output_index: Path,
    output_completion: Path,
) -> None:
    """Build the host minimap2 index once using node-local scratch."""
    workflow = load_workflow_config(config_path=config_path)
    host_reference = _required_host_reference(workflow=workflow)
    signature = task_signature(
        task={"action": "build-host-index", "run_id": workflow.run_id},
        inputs=[host_reference, workflow.config_path],
        checksum_files=True,
    )
    if completion_is_valid(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_index],
    ):
        LOGGER.info("Host index is already complete: %s", output_index)
        return
    _require_tools(action="build-host-index")
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    with scratch_workspace(scratch_root=workflow.scratch_root, label="host_index") as workspace:
        local_reference = host_reference
        stage_log = output_index.parent / "host_index.log"
        if workflow.stage_resources:
            local_reference = stage_resource(
                source=host_reference,
                destination_root=workspace / "reference",
                log_path=stage_log,
            )
        local_index = workspace / "host_reference.mmi"
        run_command(
            command=minimap2_index_command(
                reference=local_reference,
                output_index=local_index,
            ),
            log_path=stage_log,
        )
        _publish_file(source=local_index, destination=output_index, log_path=stage_log)
    write_completion(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_index],
        extra={"action": "build-host-index"},
    )


def build_minimap_index(
    *,
    config_path: Path,
    output_index: Path,
    output_completion: Path,
) -> None:
    """Build a checksum-bound classification index from the configured FASTA."""
    workflow = load_workflow_config(config_path=config_path)
    signature = task_signature(
        task={"action": "build-minimap2-index", "run_id": workflow.run_id},
        inputs=[workflow.minimap_reference, workflow.config_path],
        checksum_files=True,
    )
    if completion_is_valid(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_index],
    ):
        LOGGER.info("Classification minimap2 index is already complete: %s", output_index)
        return
    _require_tools(action="build-host-index")
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    stage_log = output_index.parent / "minimap2_index.log"
    with scratch_workspace(
        scratch_root=workflow.scratch_root,
        label="classification_minimap_index",
    ) as workspace:
        local_reference = workflow.minimap_reference
        if workflow.stage_resources:
            local_reference = stage_resource(
                source=workflow.minimap_reference,
                destination_root=workspace / "reference",
                log_path=stage_log,
            )
        local_index = workspace / "classification_reference.mmi"
        run_command(
            command=minimap2_index_command(
                reference=local_reference,
                output_index=local_index,
            ),
            log_path=stage_log,
        )
        _publish_file(source=local_index, destination=output_index, log_path=stage_log)
    write_completion(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_index],
        extra={
            "action": "build-minimap2-index",
            "reference": str(workflow.minimap_reference),
            "reference_sha256": sha256_file(path=workflow.minimap_reference),
        },
    )


def host_deplete_batch(
    *,
    config_path: Path,
    host_index: Path,
    stage_completion: Path,
) -> None:
    """Remove host reads from every sample after staging the index once."""
    workflow = load_workflow_config(config_path=config_path)
    _require_tools(action="host-deplete")
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    stage_root = stage_completion.parent
    stage_root.mkdir(parents=True, exist_ok=True)
    stage_log = stage_root / "resource_staging.log"
    with scratch_workspace(scratch_root=workflow.scratch_root, label="host_depletion") as workspace:
        local_index = host_index
        if workflow.stage_resources:
            local_index = stage_resource(
                source=host_index,
                destination_root=workspace / "host_index",
                log_path=stage_log,
            )
        for sample in workflow.samples:
            local_fastqs = tuple(
                stage_resource(
                    source=fastq_path,
                    destination_root=(
                        workspace / "inputs" / sample.sample_id / f"part_{part_number:04d}"
                    ),
                    log_path=stage_log,
                )
                for part_number, fastq_path in enumerate(sample.fastq_paths, start=1)
            )
            _host_deplete_sample(
                workflow=workflow,
                sample=sample,
                host_index=local_index,
                signature_index=host_index,
                runtime_fastqs=local_fastqs,
                workspace=workspace,
            )
    declared = [
        _host_result_directory(workflow=workflow, sample=sample) for sample in workflow.samples
    ]
    write_json_atomic(
        path=stage_completion,
        payload={
            "status": "success",
            "completed_at_utc": utc_now(),
            "stage": "host_depletion",
            "sample_count": len(declared),
            "sample_directories": [str(path) for path in declared],
        },
    )


def accept_host_removed_batch(*, config_path: Path, stage_completion: Path) -> None:
    """Prepare already host-removed reads without repeating host depletion.

    Args:
        config_path: Workflow YAML path.
        stage_completion: Declared batch completion record.

    Raises:
        ValueError: If the configuration does not declare host-removed inputs.
    """
    workflow = load_workflow_config(config_path=config_path)
    if workflow.input_read_state != "host_removed":
        raise ValueError("accept_host_removed_batch requires inputs.read_state=host_removed")
    _require_tools(action="accept-host-removed")
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    stage_root = stage_completion.parent
    stage_root.mkdir(parents=True, exist_ok=True)
    stage_log = stage_root / "resource_staging.log"
    with scratch_workspace(
        scratch_root=workflow.scratch_root,
        label="accept_host_removed",
    ) as workspace:
        for sample in workflow.samples:
            local_fastqs = tuple(
                stage_resource(
                    source=fastq_path,
                    destination_root=(
                        workspace / "inputs" / sample.sample_id / f"part_{part_number:04d}"
                    ),
                    log_path=stage_log,
                )
                for part_number, fastq_path in enumerate(sample.fastq_paths, start=1)
            )
            _accept_host_removed_sample(
                workflow=workflow,
                sample=sample,
                runtime_fastqs=local_fastqs,
                workspace=workspace,
            )
    write_json_atomic(
        path=stage_completion,
        payload={
            "status": "success",
            "completed_at_utc": utc_now(),
            "stage": "accept_host_removed",
            "sample_count": len(workflow.samples),
            "host_depletion_performed": False,
        },
    )


def classify_batch(
    *,
    config_path: Path,
    method: str,
    stage_completion: Path,
    sample_id: str | None = None,
) -> None:
    """Classify selected non-host samples after staging one method resource.

    KmerSutra is normally invoked with one sample per Snakemake job so its
    timeout and failure state are isolated. Kraken2 and Metabuli deliberately
    remain bundled to stage each large database only once.
    """
    if method not in {"kraken2", "metabuli", "minimap2", "kmersutra"}:
        raise ValueError(f"Unsupported classifier method: {method}")
    workflow = load_workflow_config(config_path=config_path)
    selected_samples = workflow.samples
    if sample_id is not None:
        selected_samples = tuple(
            sample for sample in workflow.samples if sample.sample_id == sample_id
        )
        if not selected_samples:
            raise ValueError(f"Unknown sample_id for classification: {sample_id}")
    if method == "kmersutra" and not workflow.kmersutra_enabled:
        write_json_atomic(
            path=stage_completion,
            payload={
                "status": "skipped",
                "completed_at_utc": utc_now(),
                "stage": method,
                "reason": "KmerSutra is disabled in the workflow configuration",
                "sample_count": len(selected_samples),
            },
        )
        return
    _require_tools(action=f"classify-{method}")
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    reference_fasta: Path | None = None
    if method == "kraken2":
        resource = workflow.kraken_database
    elif method == "metabuli":
        resource = workflow.metabuli_database
    elif method == "minimap2":
        resource = _resolved_minimap_index(workflow=workflow)
        reference_fasta = workflow.minimap_reference
    else:
        resource = _required_kmersutra_panel(workflow=workflow)
    stage_root = stage_completion.parent
    stage_root.mkdir(parents=True, exist_ok=True)
    stage_log = stage_root / "resource_staging.log"
    with scratch_workspace(scratch_root=workflow.scratch_root, label=method) as workspace:
        local_resource = resource
        if workflow.stage_resources:
            local_resource = stage_resource(
                source=resource,
                destination_root=workspace / "resource",
                log_path=stage_log,
            )
        local_reference = reference_fasta
        if workflow.stage_resources and reference_fasta is not None:
            local_reference = stage_resource(
                source=reference_fasta,
                destination_root=workspace / "reference_fasta",
                log_path=stage_log,
            )
        for sample in selected_samples:
            try:
                _classify_sample(
                    workflow=workflow,
                    sample=sample,
                    method=method,
                    resource=local_resource,
                    signature_resources=tuple(
                        path for path in (resource, reference_fasta) if path is not None
                    ),
                    reference_fasta=local_reference,
                    workspace=workspace,
                )
            except Exception as error:
                if method != "kmersutra" or workflow.kmersutra_failure_policy == "fail":
                    raise
                LOGGER.exception(
                    "KmerSutra failed for %s; continuing as configured",
                    sample.sample_id,
                )
                failure_path = (
                    _method_result_directory(
                        workflow=workflow,
                        sample=sample,
                        method=method,
                    )
                    / "failure.json"
                )
                write_json_atomic(
                    path=failure_path,
                    payload={
                        "status": "failed",
                        "sample_id": sample.sample_id,
                        "method": method,
                        "error": str(error),
                        "failed_at_utc": utc_now(),
                    },
                )
    failures = [
        sample.sample_id
        for sample in selected_samples
        if (
            _method_result_directory(
                workflow=workflow,
                sample=sample,
                method=method,
            )
            / "failure.json"
        ).is_file()
    ]
    write_json_atomic(
        path=stage_completion,
        payload={
            "status": "partial" if failures else "success",
            "completed_at_utc": utc_now(),
            "stage": method,
            "sample_count": len(selected_samples),
            "failed_samples": failures,
        },
    )


def aggregate_results(*, config_path: Path, completion_path: Path) -> None:
    """Create harmonised TSV summaries and a checksummed result inventory."""
    workflow = load_workflow_config(config_path=config_path)
    final_root = completion_path.parent
    final_root.mkdir(parents=True, exist_ok=True)
    sample_rows: list[dict[str, Any]] = []
    taxon_rows: list[dict[str, Any]] = []
    minimap_rows: list[dict[str, Any]] = []
    kmersutra_rows: list[dict[str, Any]] = []
    for sample in workflow.samples:
        host_summary = _read_single_tsv(
            path=_host_result_directory(workflow=workflow, sample=sample)
            / "host_removal_summary.tsv"
        )
        sample_rows.append(
            {
                **host_summary,
                "input_read_state": workflow.input_read_state,
                "host_depletion_performed": str(workflow.input_read_state == "raw").lower(),
                "kraken2_status": _method_status(
                    workflow=workflow,
                    sample=sample,
                    method="kraken2",
                ),
                "metabuli_status": _method_status(
                    workflow=workflow,
                    sample=sample,
                    method="metabuli",
                ),
                "minimap2_status": _method_status(
                    workflow=workflow,
                    sample=sample,
                    method="minimap2",
                ),
                "kmersutra_status": _method_status(
                    workflow=workflow,
                    sample=sample,
                    method="kmersutra",
                ),
            }
        )
        for method in ("kraken2", "metabuli"):
            report = (
                _method_result_directory(workflow=workflow, sample=sample, method=method)
                / "report.tsv"
            )
            taxon_rows.extend(
                _parse_classifier_report(
                    path=report,
                    sample_id=sample.sample_id,
                    method=method,
                )
            )
        minimap_rows.extend(
            _read_rows_with_method(
                path=(
                    _method_result_directory(
                        workflow=workflow,
                        sample=sample,
                        method="minimap2",
                    )
                    / "taxon_report.tsv"
                ),
                sample_id=sample.sample_id,
                method="minimap2",
            )
        )
        calls = (
            _method_result_directory(
                workflow=workflow,
                sample=sample,
                method="kmersutra",
            )
            / "species_detection_calls.tsv"
        )
        if calls.is_file():
            kmersutra_rows.extend(_read_rows_with_prefix(path=calls, sample_id=sample.sample_id))
        else:
            kmersutra_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "method": "kmersutra",
                    "workflow_status": _method_status(
                        workflow=workflow,
                        sample=sample,
                        method="kmersutra",
                    ),
                }
            )

    summary_path = final_root / "sample_summary.tsv"
    taxon_path = final_root / "classifier_taxon_reports.tsv.gz"
    calls_path = final_root / "kmersutra_species_calls.tsv.gz"
    minimap_path = final_root / "minimap2_taxon_reports.tsv.gz"
    _write_tsv(
        path=summary_path,
        rows=sample_rows,
        fieldnames=(
            "sample_id",
            "input_reads",
            "non_host_reads",
            "host_reads_removed",
            "host_fraction",
            "input_read_state",
            "host_depletion_performed",
            "kraken2_status",
            "metabuli_status",
            "minimap2_status",
            "kmersutra_status",
        ),
    )
    _write_tsv(
        path=taxon_path,
        rows=taxon_rows,
        fieldnames=(
            "sample_id",
            "method",
            "fraction_percent",
            "clade_reads",
            "direct_reads",
            "rank_code",
            "tax_id",
            "taxon_name",
        ),
    )
    calls_fields = tuple(dict.fromkeys(key for row in kmersutra_rows for key in row))
    _write_tsv(path=calls_path, rows=kmersutra_rows, fieldnames=calls_fields)
    minimap_fields = tuple(dict.fromkeys(key for row in minimap_rows for key in row)) or (
        "sample_id",
        "method",
        "tax_id",
        "taxon_name",
        "best_read_count",
        "ambiguous_best_read_count",
        "alignment_count",
        "reference_count",
        "min_mapq",
        "min_alignment",
    )
    _write_tsv(path=minimap_path, rows=minimap_rows, fieldnames=minimap_fields)
    readme_path = final_root / "README.txt"
    readme_path.write_text(
        "Nanopore real-data workflow final results\n\n"
        "sample_summary.tsv contains host-removal totals and stage status.\n"
        "classifier_taxon_reports.tsv.gz harmonises Kraken2 and Metabuli reports.\n"
        "minimap2_taxon_reports.tsv.gz summarises controlled-reference best alignments.\n"
        "kmersutra_species_calls.tsv.gz combines KmerSutra species calls.\n"
        "SHA256SUMS.tsv records final-file checksums.\n",
        encoding="utf-8",
    )
    checksum_path = final_root / "SHA256SUMS.tsv"
    checksum_rows = [
        {"sha256": sha256_file(path=path), "file": path.name}
        for path in (summary_path, taxon_path, minimap_path, calls_path, readme_path)
    ]
    _write_tsv(
        path=checksum_path,
        rows=checksum_rows,
        fieldnames=("sha256", "file"),
    )
    write_json_atomic(
        path=completion_path,
        payload={
            "status": "success",
            "completed_at_utc": utc_now(),
            "run_id": workflow.run_id,
            "sample_count": len(workflow.samples),
            "outputs": [
                str(path)
                for path in (
                    summary_path,
                    taxon_path,
                    minimap_path,
                    calls_path,
                    readme_path,
                    checksum_path,
                )
            ],
        },
    )


def run_snakemake(
    *,
    config_path: Path,
    snakefile: Path,
    profile: Path | None,
    cores: str,
    jobs: int,
    dry_run: bool,
    unlock: bool,
    extra_arguments: Sequence[str] = (),
) -> int:
    """Run Snakemake with safe defaults and one database-staging resource."""
    if cores != "all":
        try:
            numeric_cores = int(cores)
        except ValueError as error:
            raise ValueError("Snakemake cores must be a positive integer or 'all'") from error
        if numeric_cores <= 0:
            raise ValueError("Snakemake cores must be positive")
    if jobs <= 0:
        raise ValueError("Snakemake jobs must be positive")
    command = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(config_path),
        "--config",
        f"workflow_config_path={config_path.resolve()}",
        "--cores",
        cores,
        "--jobs",
        str(jobs),
        "--resources",
        "database_stage=1",
        "kmersutra_worker=1",
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if profile is not None:
        command.extend(["--profile", str(profile)])
    if dry_run:
        command.append("--dry-run")
    if unlock:
        command.append("--unlock")
    command.extend(extra_arguments)
    return subprocess.run(command, check=False).returncode


def _host_deplete_sample(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    host_index: Path,
    signature_index: Path,
    runtime_fastqs: Sequence[Path],
    workspace: Path,
) -> None:
    final = _host_result_directory(workflow=workflow, sample=sample)
    output_fastq = final / "non_host.fastq.gz"
    summary = final / "host_removal_summary.tsv"
    metadata = final / "metadata.json"
    completion = final / "complete.json"
    signature = task_signature(
        task={
            "action": "host-deplete",
            "sample_id": sample.sample_id,
            "threads": workflow.threads_host,
        },
        inputs=[*sample.fastq_paths, signature_index, workflow.config_path],
        checksum_files=workflow.checksum_inputs,
    )
    if completion_is_valid(
        completion_path=completion,
        signature=signature,
        outputs=[output_fastq, summary, metadata],
    ):
        LOGGER.info("Skipping valid host-depletion result: %s", sample.sample_id)
        return
    local = workspace / "samples" / sample.sample_id
    shutil.rmtree(local, ignore_errors=True)
    local.mkdir(parents=True)
    local_bam = local / "host_alignment.bam"
    local_fastq = local / "non_host.fastq.gz"
    log_path = local / "host_depletion.log"
    try:
        run_pipeline(
            commands=[
                minimap2_host_command(
                    host_index=host_index,
                    fastq_paths=runtime_fastqs,
                    threads=workflow.threads_host,
                ),
                samtools_bam_command(
                    output_bam=local_bam,
                    threads=workflow.threads_host,
                ),
            ],
            log_path=log_path,
        )
        total_reads = _capture_integer(
            command=samtools_count_command(
                bam_path=local_bam,
                threads=workflow.threads_host,
                unmapped_only=False,
            )
        )
        non_host_reads = _capture_integer(
            command=samtools_count_command(
                bam_path=local_bam,
                threads=workflow.threads_host,
                unmapped_only=True,
            )
        )
        run_pipeline(
            commands=[
                samtools_non_host_fastq_command(
                    bam_path=local_bam,
                    threads=workflow.threads_host,
                ),
                pigz_command(threads=workflow.threads_host),
            ],
            log_path=log_path,
            stdout_path=local_fastq,
        )
        removed = total_reads - non_host_reads
        _write_tsv(
            path=local / "host_removal_summary.tsv",
            rows=[
                {
                    "sample_id": sample.sample_id,
                    "input_reads": total_reads,
                    "non_host_reads": non_host_reads,
                    "host_reads_removed": removed,
                    "host_fraction": f"{removed / total_reads:.8f}"
                    if total_reads
                    else "0.00000000",
                }
            ],
            fieldnames=(
                "sample_id",
                "input_reads",
                "non_host_reads",
                "host_reads_removed",
                "host_fraction",
            ),
        )
        write_json_atomic(
            path=local / "metadata.json",
            payload={
                "sample_id": sample.sample_id,
                "status": "success",
                "input_fastqs": [str(path) for path in sample.fastq_paths],
                "host_index": str(signature_index),
                "minimap2_version": _version(command=["minimap2", "--version"]),
                "samtools_version": _version(command=["samtools", "--version"]),
                "completed_at_utc": utc_now(),
            },
        )
        local_bam.unlink(missing_ok=True)
        publish_directory(source=local, destination=final, log_path=log_path)
        write_completion(
            completion_path=completion,
            signature=signature,
            outputs=[output_fastq, summary, metadata],
            extra={"sample_id": sample.sample_id, "action": "host-deplete"},
        )
    except Exception:
        _publish_failed_attempt(source=local, destination=final)
        raise


def _accept_host_removed_sample(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    runtime_fastqs: Sequence[Path],
    workspace: Path,
) -> None:
    """Validate and publish classification-ready host-removed FASTQ parts."""
    final = _host_result_directory(workflow=workflow, sample=sample)
    output_fastq = final / "non_host.fastq.gz"
    summary = final / "host_removal_summary.tsv"
    metadata = final / "metadata.json"
    completion = final / "complete.json"
    signature = task_signature(
        task={
            "action": "accept-host-removed",
            "sample_id": sample.sample_id,
            "input_read_state": workflow.input_read_state,
        },
        inputs=[*sample.fastq_paths, workflow.config_path],
        checksum_files=workflow.checksum_inputs,
    )
    if completion_is_valid(
        completion_path=completion,
        signature=signature,
        outputs=[output_fastq, summary, metadata],
    ):
        LOGGER.info("Skipping valid host-removed input preparation: %s", sample.sample_id)
        return
    local = workspace / "samples" / sample.sample_id
    shutil.rmtree(local, ignore_errors=True)
    local.mkdir(parents=True)
    try:
        read_count = _merge_fastq_parts(
            input_paths=runtime_fastqs,
            output_path=local / "non_host.fastq.gz",
        )
        _write_tsv(
            path=local / "host_removal_summary.tsv",
            rows=[
                {
                    "sample_id": sample.sample_id,
                    "input_reads": read_count,
                    "non_host_reads": read_count,
                    "host_reads_removed": 0,
                    "host_fraction": "0.00000000",
                    "input_read_state": "host_removed",
                    "host_depletion_performed": "false",
                }
            ],
            fieldnames=(
                "sample_id",
                "input_reads",
                "non_host_reads",
                "host_reads_removed",
                "host_fraction",
                "input_read_state",
                "host_depletion_performed",
            ),
        )
        write_json_atomic(
            path=local / "metadata.json",
            payload={
                "sample_id": sample.sample_id,
                "status": "success",
                "input_fastqs": [str(path) for path in sample.fastq_paths],
                "input_read_state": "host_removed",
                "host_depletion_performed": False,
                "read_count": read_count,
                "completed_at_utc": utc_now(),
            },
        )
        publish_directory(
            source=local,
            destination=final,
            log_path=local / "input_preparation.log",
        )
        write_completion(
            completion_path=completion,
            signature=signature,
            outputs=[output_fastq, summary, metadata],
            extra={
                "sample_id": sample.sample_id,
                "action": "accept-host-removed",
                "host_depletion_performed": False,
            },
        )
    except Exception:
        _publish_failed_attempt(source=local, destination=final)
        raise


def _classify_sample(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    method: str,
    resource: Path,
    signature_resources: Sequence[Path],
    reference_fasta: Path | None,
    workspace: Path,
) -> None:
    input_fastq = _host_result_directory(workflow=workflow, sample=sample) / "non_host.fastq.gz"
    if not input_fastq.is_file() or input_fastq.stat().st_size == 0:
        raise FileNotFoundError(f"Non-host FASTQ is missing for {sample.sample_id}: {input_fastq}")
    final = _method_result_directory(workflow=workflow, sample=sample, method=method)
    completion = final / "complete.json"
    expected = _method_expected_outputs(
        directory=final,
        method=method,
        keep_per_read=workflow.keep_per_read_classifications,
    )
    task = _method_task_settings(workflow=workflow, sample=sample, method=method)
    signature = task_signature(
        task=task,
        inputs=[input_fastq, *signature_resources, workflow.config_path],
        checksum_files=False,
    )
    if completion_is_valid(
        completion_path=completion,
        signature=signature,
        outputs=expected,
    ):
        LOGGER.info("Skipping valid %s result: %s", method, sample.sample_id)
        return
    local = workspace / "samples" / sample.sample_id
    shutil.rmtree(local, ignore_errors=True)
    local.mkdir(parents=True)
    log_path = local / f"{method}.log"
    runtime_input = stage_resource(
        source=input_fastq,
        destination_root=workspace / "inputs" / sample.sample_id,
        log_path=log_path,
    )
    try:
        if method == "kraken2":
            _run_kraken2(
                workflow=workflow,
                sample=sample,
                input_fastq=runtime_input,
                database=resource,
                output_directory=local,
                log_path=log_path,
            )
        elif method == "metabuli":
            _run_metabuli(
                workflow=workflow,
                sample=sample,
                input_fastq=runtime_input,
                database=resource,
                output_directory=local,
                log_path=log_path,
            )
        elif method == "minimap2":
            if reference_fasta is None:
                raise ValueError("Classification minimap2 requires a reference FASTA")
            _run_minimap2(
                workflow=workflow,
                sample=sample,
                input_fastq=runtime_input,
                reference_index=resource,
                reference_fasta=reference_fasta,
                output_directory=local,
                log_path=log_path,
            )
        else:
            _run_kmersutra(
                workflow=workflow,
                sample=sample,
                input_fastq=runtime_input,
                panel=resource,
                output_directory=local,
                log_path=log_path,
            )
        write_json_atomic(
            path=local / "metadata.json",
            payload={
                "sample_id": sample.sample_id,
                "method": method,
                "status": "success",
                "input_fastq": str(input_fastq),
                "resources": [str(path) for path in signature_resources],
                "settings": task,
                "tool_version": _method_version(method=method),
                "completed_at_utc": utc_now(),
            },
        )
        publish_directory(source=local, destination=final, log_path=log_path)
        final_expected = _method_expected_outputs(
            directory=final,
            method=method,
            keep_per_read=workflow.keep_per_read_classifications,
        )
        write_completion(
            completion_path=completion,
            signature=signature,
            outputs=final_expected,
            extra={"sample_id": sample.sample_id, "method": method},
        )
    except Exception:
        _publish_failed_attempt(source=local, destination=final)
        raise


def _run_kraken2(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    input_fastq: Path,
    database: Path,
    output_directory: Path,
    log_path: Path,
) -> None:
    raw = output_directory / "classifications.tsv"
    report = output_directory / "report.tsv"
    run_command(
        command=kraken2_command(
            input_fastq=input_fastq,
            database=database,
            classifications=raw,
            report=report,
            threads=workflow.threads_kraken2,
            confidence=workflow.kraken_confidence,
        ),
        log_path=log_path,
    )
    _handle_per_read_output(
        raw_path=raw,
        keep=workflow.keep_per_read_classifications,
        threads=workflow.threads_kraken2,
        log_path=log_path,
    )
    if not report.is_file() or report.stat().st_size == 0:
        raise RuntimeError(f"Kraken2 produced no report for {sample.sample_id}")


def _run_metabuli(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    input_fastq: Path,
    database: Path,
    output_directory: Path,
    log_path: Path,
) -> None:
    run_command(
        command=metabuli_command(
            input_fastq=input_fastq,
            database=database,
            output_directory=output_directory,
            output_prefix="metabuli",
            threads=workflow.threads_metabuli,
            maximum_ram_gb=workflow.metabuli_max_ram_gb,
            minimum_score=workflow.metabuli_min_score,
        ),
        log_path=log_path,
    )
    raw = output_directory / "metabuli_classifications.tsv"
    generated_report = output_directory / "metabuli_report.tsv"
    report = output_directory / "report.tsv"
    if generated_report.is_file():
        generated_report.replace(report)
    _handle_per_read_output(
        raw_path=raw,
        keep=workflow.keep_per_read_classifications,
        threads=workflow.threads_metabuli,
        log_path=log_path,
    )
    if not report.is_file() or report.stat().st_size == 0:
        raise RuntimeError(f"Metabuli produced no report for {sample.sample_id}")


def _run_minimap2(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    input_fastq: Path,
    reference_index: Path,
    reference_fasta: Path,
    output_directory: Path,
    log_path: Path,
) -> None:
    """Map reads to the controlled reference and write filtered taxon summaries."""
    paf_path = output_directory / "alignments.paf.gz"
    run_pipeline(
        commands=[
            minimap2_classification_command(
                reference_index=reference_index,
                input_fastq=input_fastq,
                threads=workflow.threads_minimap2,
            ),
            pigz_command(threads=workflow.threads_minimap2),
        ],
        log_path=log_path,
        stdout_path=paf_path,
    )
    summarise_minimap_paf(
        paf_path=paf_path,
        reference_fasta=reference_fasta,
        taxon_report_path=output_directory / "taxon_report.tsv",
        mapping_summary_path=output_directory / "mapping_summary.tsv",
        sample_id=sample.sample_id,
        minimum_mapq=workflow.minimap_min_mapq,
        minimum_alignment=workflow.minimap_min_alignment,
    )
    if not workflow.keep_per_read_classifications:
        paf_path.unlink()


def _run_kmersutra(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    input_fastq: Path,
    panel: Path,
    output_directory: Path,
    log_path: Path,
) -> None:
    run_command(
        command=kmersutra_command(
            input_fastq=input_fastq,
            panel=panel,
            sample_id=sample.sample_id,
            output_directory=output_directory,
            threads=workflow.threads_kmersutra,
            screen_preset=workflow.kmersutra_screen_preset,
            call_preset=workflow.kmersutra_call_preset,
            same_genus_fraction=workflow.kmersutra_same_genus_fraction,
            write_parquet=workflow.kmersutra_write_parquet,
        ),
        log_path=log_path,
        timeout_seconds=workflow.kmersutra_timeout_minutes * 60,
    )
    calls = output_directory / "species_detection_calls.tsv"
    if not calls.is_file() or calls.stat().st_size == 0:
        raise RuntimeError(f"KmerSutra produced no species calls for {sample.sample_id}")


def _handle_per_read_output(
    *,
    raw_path: Path,
    keep: bool,
    threads: int,
    log_path: Path,
) -> None:
    if not raw_path.is_file():
        raise RuntimeError(f"Classifier did not create its per-read output: {raw_path}")
    if not keep:
        raw_path.unlink()
        return
    compressed = raw_path.with_suffix(raw_path.suffix + ".gz")
    run_command(
        command=pigz_command(threads=threads) + [str(raw_path)],
        log_path=log_path,
        stdout_path=compressed,
    )
    raw_path.unlink()


def _method_expected_outputs(
    *,
    directory: Path,
    method: str,
    keep_per_read: bool,
) -> list[Path]:
    outputs = [directory / "metadata.json"]
    if method in {"kraken2", "metabuli"}:
        outputs.append(directory / "report.tsv")
        if keep_per_read:
            name = (
                "classifications.tsv.gz"
                if method == "kraken2"
                else "metabuli_classifications.tsv.gz"
            )
            outputs.append(directory / name)
    elif method == "minimap2":
        outputs.extend(
            [
                directory / "taxon_report.tsv",
                directory / "mapping_summary.tsv",
            ]
        )
        if keep_per_read:
            outputs.append(directory / "alignments.paf.gz")
    else:
        outputs.extend(
            [
                directory / "species_detection_calls.tsv",
                directory / "sample_species_kmer_evidence.tsv",
                directory / "sample_lineage_interpretation.tsv",
            ]
        )
    return outputs


def _method_task_settings(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    method: str,
) -> dict[str, Any]:
    common: dict[str, Any] = {
        "action": f"classify-{method}",
        "sample_id": sample.sample_id,
    }
    if method == "kraken2":
        common.update(
            threads=workflow.threads_kraken2,
            confidence=workflow.kraken_confidence,
        )
    elif method == "metabuli":
        common.update(
            threads=workflow.threads_metabuli,
            min_score=workflow.metabuli_min_score,
            max_ram_gb=workflow.metabuli_max_ram_gb,
        )
    elif method == "minimap2":
        common.update(
            threads=workflow.threads_minimap2,
            min_mapq=workflow.minimap_min_mapq,
            min_alignment=workflow.minimap_min_alignment,
        )
    else:
        common.update(
            threads=workflow.threads_kmersutra,
            screen_preset=workflow.kmersutra_screen_preset,
            call_preset=workflow.kmersutra_call_preset,
            same_genus_fraction=workflow.kmersutra_same_genus_fraction,
            write_parquet=workflow.kmersutra_write_parquet,
        )
    return common


def _host_result_directory(*, workflow: WorkflowConfig, sample: Sample) -> Path:
    return workflow.output_directory / "01_host_depletion" / sample.sample_id


def _method_result_directory(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    method: str,
) -> Path:
    return workflow.output_directory / "02_classification" / method / sample.sample_id


def _method_status(*, workflow: WorkflowConfig, sample: Sample, method: str) -> str:
    if method == "kmersutra" and not workflow.kmersutra_enabled:
        return "disabled"
    path = (
        _method_result_directory(
            workflow=workflow,
            sample=sample,
            method=method,
        )
        / "complete.json"
    )
    if not path.is_file():
        failure = path.with_name("failure.json")
        if failure.is_file():
            return "failed"
        return "missing"
    try:
        return str(json.loads(path.read_text(encoding="utf-8")).get("status", "invalid"))
    except (json.JSONDecodeError, OSError):
        return "invalid"


def _parse_classifier_report(
    *,
    path: Path,
    sample_id: str,
    method: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{method} report line {line_number} has fewer than six fields")
            rows.append(
                {
                    "sample_id": sample_id,
                    "method": method,
                    "fraction_percent": fields[0].strip(),
                    "clade_reads": fields[1].strip(),
                    "direct_reads": fields[2].strip(),
                    "rank_code": fields[3].strip(),
                    "tax_id": fields[4].strip(),
                    "taxon_name": fields[5].strip(),
                }
            )
    return rows


def _read_rows_with_prefix(*, path: Path, sample_id: str) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [{"sample_id": sample_id, "method": "kmersutra", **row} for row in rows]


def _read_rows_with_method(
    *,
    path: Path,
    sample_id: str,
    method: str,
) -> list[dict[str, Any]]:
    """Read a TSV and enforce its workflow sample and method fields."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [{**row, "sample_id": sample_id, "method": method} for row in rows]


def _read_single_tsv(*, path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one data row in {path}")
    return rows[0]


def _merge_fastq_parts(*, input_paths: Sequence[Path], output_path: Path) -> int:
    """Validate, merge and gzip one or more FASTQ parts.

    Args:
        input_paths: Ordered FASTQ or FASTQ.GZ parts for one sample.
        output_path: Gzip-compressed merged FASTQ path.

    Returns:
        Number of validated FASTQ records written.

    Raises:
        ValueError: If inputs are absent or a FASTQ record is malformed.
    """
    if not input_paths:
        raise ValueError("At least one FASTQ part is required")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    read_count = 0
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as output_handle:
        for input_path in input_paths:
            opener = gzip.open if str(input_path).lower().endswith(".gz") else open
            with opener(input_path, "rt", encoding="utf-8", newline="") as input_handle:
                while True:
                    lines = [input_handle.readline() for _ in range(4)]
                    if lines[0] == "":
                        if any(lines[1:]):
                            raise ValueError(f"Malformed FASTQ end-of-file: {input_path}")
                        break
                    if any(line == "" for line in lines[1:]):
                        raise ValueError(f"Truncated FASTQ record: {input_path}")
                    if not lines[0].startswith("@") or not lines[2].startswith("+"):
                        raise ValueError(f"Malformed FASTQ record: {input_path}")
                    if len(lines[1].rstrip("\r\n")) != len(lines[3].rstrip("\r\n")):
                        raise ValueError(f"FASTQ sequence and quality lengths differ: {input_path}")
                    output_handle.writelines(lines)
                    read_count += 1
    return read_count


def _write_tsv(
    *,
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _publish_file(*, source: Path, destination: Path, log_path: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.partial.{os.getpid()}")
    try:
        run_command(
            command=["rsync", "--archive", str(source), str(temporary)],
            log_path=log_path,
        )
        if not temporary.is_file() or temporary.stat().st_size != source.stat().st_size:
            raise RuntimeError(f"Published file validation failed: {destination}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _publish_failed_attempt(*, source: Path, destination: Path) -> None:
    if not source.is_dir() or not any(source.iterdir()):
        return
    logs = sorted(path for path in source.glob("*.log") if path.is_file())
    if not logs:
        return
    failed_root = destination.parent / "failed_attempts"
    failed_root.mkdir(parents=True, exist_ok=True)
    failed = failed_root / f"{destination.name}.{os.getpid()}"
    if failed.exists():
        shutil.rmtree(failed)
    failed.mkdir()
    for log_path in logs:
        shutil.copy2(log_path, failed / log_path.name)
    write_json_atomic(
        path=failed / "failure.json",
        payload={"status": "failed", "failed_at_utc": utc_now()},
    )


def _validate_kraken_database(*, database: Path) -> None:
    missing = [
        name for name in ("hash.k2d", "opts.k2d", "taxo.k2d") if not (database / name).is_file()
    ]
    if missing:
        raise ValueError(f"Kraken2 database is missing required files: {', '.join(missing)}")


def _validate_non_empty_directory(*, label: str, path: Path) -> None:
    if not path.is_dir() or not any(path.iterdir()):
        raise ValueError(f"{label} is missing or empty: {path}")


def _validate_panel(*, panel: Path) -> None:
    opener = gzip.open if str(panel).lower().endswith(".gz") else open
    with opener(panel, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    required = {"kmer", "k", "species_name"}
    if not required.issubset(header):
        raise ValueError(
            f"KmerSutra panel lacks required columns: {sorted(required - set(header))}"
        )


def _required_kmersutra_panel(*, workflow: WorkflowConfig) -> Path:
    panel = workflow.kmersutra_panel
    if panel is None:
        raise ValueError("KmerSutra is enabled but no panel is configured")
    return panel


def _required_host_reference(*, workflow: WorkflowConfig) -> Path:
    """Return the host reference required for raw-read processing."""
    reference = workflow.host_reference
    if reference is None:
        raise ValueError("Raw-read processing requires host.reference")
    return reference


def _resolved_minimap_index(*, workflow: WorkflowConfig) -> Path:
    """Return a configured or workflow-built classification minimap2 index."""
    return (
        workflow.minimap_index
        if workflow.minimap_index is not None
        else workflow.output_directory / "00_preflight" / "classification_reference.mmi"
    )


def _validate_fastq_prefix(*, path: Path) -> None:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = [handle.readline() for _ in range(4)]
    if any(line == "" for line in lines):
        raise ValueError(f"FASTQ contains no complete first record: {path}")
    if not lines[0].startswith("@") or not lines[2].startswith("+"):
        raise ValueError(f"FASTQ first record is malformed: {path}")
    sequence = lines[1].rstrip("\r\n")
    quality = lines[3].rstrip("\r\n")
    if len(sequence) != len(quality):
        raise ValueError(f"FASTQ first-record sequence and quality lengths differ: {path}")


def _require_tools(*, action: str) -> None:
    missing = [name for name in required_executables(action=action) if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"Required executables are not on PATH: {', '.join(missing)}")


def _capture_integer(*, command: Sequence[str]) -> int:
    value = capture_output(command=command)
    try:
        return int(value)
    except ValueError as error:
        raise RuntimeError(f"Expected integer command output, received: {value!r}") from error


def _version(*, command: Sequence[str]) -> str:
    try:
        return capture_output(command=command).splitlines()[0]
    except (RuntimeError, IndexError):
        return "unavailable"


def _method_version(*, method: str) -> str:
    if method == "kraken2":
        return _version(command=["kraken2", "--version"])
    if method == "metabuli":
        return _version(command=["metabuli", "version"])
    if method == "minimap2":
        return _version(command=["minimap2", "--version"])
    return _kmersutra_version()


def _kmersutra_version() -> str:
    return _version(
        command=[
            sys.executable,
            "-c",
            "import kmersutra; print(kmersutra.__version__)",
        ]
    )
