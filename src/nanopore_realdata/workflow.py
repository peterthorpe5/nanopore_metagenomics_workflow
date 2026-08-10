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
from nanopore_realdata.pcr import (
    PcrTruth,
    build_pcr_concordance,
    canonical_species_name,
    expected_pcr_species,
    load_pcr_truth,
    pcr_truth_rows,
)
from nanopore_realdata.reference import (
    build_controlled_reference,
    load_genome_sources,
    validate_required_reference_species,
    validate_required_species,
    validate_single_part_index_log,
)
from nanopore_realdata.reporting import (
    METHODS,
    build_normalised_evidence,
    generate_html_reports,
    serialisable_report_data,
)
from nanopore_realdata.runtime import (
    CommandTimeoutError,
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


def _load_optional_pcr_truth(*, workflow: WorkflowConfig) -> tuple[PcrTruth, ...]:
    """Load independent truth when configured, otherwise return no records.

    Args:
        workflow: Validated workflow configuration.

    Returns:
        PCR truth records in sample-manifest order, or an empty tuple when the
        run performs classification without an independent PCR comparison.
    """
    if workflow.pcr_truth_path is None:
        return ()
    return load_pcr_truth(
        path=workflow.pcr_truth_path,
        samples=workflow.samples,
    )


def _validate_pcr_reference_contract(
    *,
    workflow: WorkflowConfig,
    truth_records: Sequence[PcrTruth],
) -> None:
    """Require benchmark truth species to be declared in minimap2 settings.

    Classification is configured independently of PCR. When PCR evaluation is
    enabled, this guard prevents an expected species from being added to the
    reference contract implicitly or omitted accidentally.

    Args:
        workflow: Validated workflow configuration.
        truth_records: Optional independent PCR records.

    Raises:
        ValueError: If a PCR-positive comparison species is not explicitly
            listed under ``minimap2.required_species``.
    """
    configured = {
        canonical_species_name(value=species).casefold()
        for species in workflow.minimap_required_species
    }
    missing = [
        species
        for species in expected_pcr_species(records=truth_records)
        if species.casefold() not in configured
    ]
    if missing:
        raise ValueError(
            "PCR comparison species must be declared explicitly in "
            "minimap2.required_species: " + "; ".join(missing)
        )


def preflight(*, config_path: Path, output_path: Path) -> None:
    """Validate configuration, input records, databases and executables.

    Args:
        config_path: Workflow YAML path.
        output_path: Declared preflight JSON output.
    """
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    truth_records = _load_optional_pcr_truth(workflow=workflow)
    _validate_pcr_reference_contract(
        workflow=workflow,
        truth_records=truth_records,
    )
    minimap_sources = None
    if workflow.minimap_genome_config is not None:
        minimap_sources = load_genome_sources(
            config_path=workflow.minimap_genome_config,
        )
        validate_required_species(
            sources=minimap_sources,
            required_species=workflow.minimap_required_species,
        )
    # Host preparation is a core dependency: without classification-ready
    # reads no classifier can run. Classifier readiness is assessed separately
    # and recorded rather than raised, allowing healthy branches to continue.
    actions: list[str] = []
    if workflow.input_read_state == "raw":
        actions.extend(("build-host-index", "host-deplete"))
    else:
        actions.append("accept-host-removed")
    for action in actions:
        _require_tools(action=action)
    classifier_readiness = [
        _classifier_readiness(
            workflow=workflow,
            method=method,
            required_species=workflow.minimap_required_species,
        )
        for method in METHODS
    ]
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
    _write_tsv(
        path=output_path.parent / "classifier_readiness.tsv",
        rows=classifier_readiness,
        fieldnames=("method", "status", "message", "resource"),
    )
    _write_tsv(
        path=output_path.parent / "resolved_pcr_truth.tsv",
        rows=pcr_truth_rows(records=truth_records),
        fieldnames=(
            "sample_id",
            "pcr_status",
            "pcr_species",
            "pcr_species_canonical",
            "pcr_assay_or_source",
            "pcr_notes",
            "include_in_primary_comparison",
        ),
    )
    payload = {
        "status": "success",
        "checked_at_utc": utc_now(),
        "run_id": workflow.run_id,
        "sample_count": len(workflow.samples),
        "pcr_evaluation_configured": workflow.pcr_truth_path is not None,
        "primary_pcr_sample_count": sum(
            record.include_in_primary_comparison for record in truth_records
        ),
        "fastq_part_count": len(sample_rows),
        "config": str(workflow.config_path),
        "config_sha256": sha256_file(path=workflow.config_path),
        "resources": {
            "input_read_state": workflow.input_read_state,
            "minimap2_required_species": list(workflow.minimap_required_species),
            "host_reference": (
                metadata_fingerprint(
                    paths=[_required_host_reference(workflow=workflow)],
                    checksum_files=True,
                )
                if workflow.input_read_state == "raw"
                else "not_applicable_already_host_removed"
            ),
            "kraken2_database": _safe_resource_fingerprint(
                path=workflow.kraken_database,
                checksum_file=False,
            ),
            "metabuli_database": _safe_resource_fingerprint(
                path=workflow.metabuli_database,
                checksum_file=False,
            ),
            "pcr_truth": (
                _safe_resource_fingerprint(
                    path=workflow.pcr_truth_path,
                    checksum_file=True,
                )
                if workflow.pcr_truth_path is not None
                else "not_configured"
            ),
            "minimap2_reference": (
                "built_by_workflow_from_genome_config"
                if workflow.minimap_genome_config is not None
                else _safe_resource_fingerprint(
                    path=workflow.minimap_reference,
                    checksum_file=True,
                )
            ),
            "minimap2_genome_config": (
                _safe_resource_fingerprint(
                    path=workflow.minimap_genome_config,
                    checksum_file=True,
                )
                if workflow.minimap_genome_config is not None
                else "not_applicable_prebuilt_reference"
            ),
            "minimap2_source_genome_count": (
                len(minimap_sources) if minimap_sources is not None else "not_applicable"
            ),
            "minimap2_index": "built_by_workflow_with_single_part_validation",
            "kmersutra_panel": (
                _safe_resource_fingerprint(
                    path=_required_kmersutra_panel(workflow=workflow),
                    checksum_file=True,
                )
                if workflow.kmersutra_enabled
                else "disabled"
            ),
        },
        "software": tools,
        "classifier_readiness": classifier_readiness,
    }
    write_json_atomic(path=output_path, payload=payload)


def build_minimap_reference(
    *,
    config_path: Path,
    output_reference: Path,
    output_manifest: Path,
    output_completion: Path,
) -> None:
    """Build a bounded minimap2 reference from configured genome sources.

    Args:
        config_path: Workflow YAML path.
        output_reference: Controlled FASTA destination.
        output_manifest: Sequence-level provenance table destination.
        output_completion: Durable success record.
    """
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    if workflow.minimap_genome_config is None:
        raise ValueError("build_minimap_reference requires minimap2.genome_config")
    sources = load_genome_sources(config_path=workflow.minimap_genome_config)
    validate_required_species(
        sources=sources,
        required_species=workflow.minimap_required_species,
    )
    signature = task_signature(
        task={
            "action": "build-minimap2-reference",
            "run_id": workflow.run_id,
            "maximum_reference_bases": workflow.minimap_maximum_reference_bases,
            "required_species": list(workflow.minimap_required_species),
        },
        inputs=[
            workflow.config_path,
            workflow.minimap_genome_config,
            *(source.genome_fasta for source in sources),
        ],
        checksum_files=False,
    )
    if completion_is_valid(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_reference, output_manifest],
    ):
        LOGGER.info("Controlled minimap2 reference is already complete: %s", output_reference)
        return
    validate_scratch(
        scratch_root=workflow.scratch_root,
        minimum_gb=workflow.minimum_scratch_gb,
    )
    with scratch_workspace(
        scratch_root=workflow.scratch_root,
        label="controlled_minimap_reference",
    ) as workspace:
        local_reference = workspace / "controlled_minimap_reference.fa"
        local_manifest = workspace / "controlled_minimap_reference.manifest.tsv"
        stats = build_controlled_reference(
            sources=sources,
            output_fasta=local_reference,
            output_manifest=local_manifest,
            maximum_reference_bases=workflow.minimap_maximum_reference_bases,
        )
        publish_log = output_reference.parent / "controlled_minimap_reference.publish.log"
        _publish_file(
            source=local_reference,
            destination=output_reference,
            log_path=publish_log,
        )
        _publish_file(
            source=local_manifest,
            destination=output_manifest,
            log_path=publish_log,
        )
    write_completion(
        completion_path=output_completion,
        signature=signature,
        outputs=[output_reference, output_manifest],
        extra={
            "action": "build-minimap2-reference",
            "genome_config": str(workflow.minimap_genome_config),
            "genome_count": stats.genome_count,
            "reference_record_count": stats.reference_record_count,
            "total_bases": stats.total_bases,
            "reference_sha256": stats.reference_sha256,
            "manifest_sha256": stats.manifest_sha256,
            "required_species": list(workflow.minimap_required_species),
        },
    )


def build_host_index(
    *,
    config_path: Path,
    output_index: Path,
    output_completion: Path,
) -> None:
    """Build the host minimap2 index once using node-local scratch."""
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
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
    """Build a checksum-bound classification index without blocking other tools.

    A failed index build writes a terminal minimap2 readiness record when that
    classifier uses ``failure_policy: continue``. The downstream minimap2 rule
    will record matching per-sample failures, while aggregation remains able to
    report Kraken2, Metabuli and KmerSutra results.
    """
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    stage_log = output_index.parent / "minimap2_index.log"
    try:
        if not workflow.minimap_reference.is_file():
            raise FileNotFoundError(
                f"Classification minimap2 reference is missing: {workflow.minimap_reference}"
            )
        validate_required_reference_species(
            reference_fasta=workflow.minimap_reference,
            required_species=workflow.minimap_required_species,
        )
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
        stage_log.unlink(missing_ok=True)
        _require_tools(action="build-host-index")
        validate_scratch(
            scratch_root=workflow.scratch_root,
            minimum_gb=workflow.minimum_scratch_gb,
        )
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
                    index_batch_size_bases=workflow.minimap_index_batch_size_bases,
                ),
                log_path=stage_log,
                timeout_seconds=_internal_timeout_seconds(
                    runtime_minutes=workflow.runtime_minimap2_minutes
                ),
            )
            index_details = validate_single_part_index_log(
                log_path=stage_log,
                maximum_reference_bases=workflow.minimap_maximum_reference_bases,
            )
            index_size = local_index.stat().st_size
            if index_size > workflow.minimap_maximum_index_bytes:
                raise ValueError(
                    "minimap2 index exceeds the configured hard size limit: "
                    f"{index_size} > {workflow.minimap_maximum_index_bytes} bytes"
                )
            _publish_file(
                source=local_index,
                destination=output_index,
                log_path=stage_log,
            )
        write_completion(
            completion_path=output_completion,
            signature=signature,
            outputs=[output_index],
            extra={
                "action": "build-minimap2-index",
                "reference": str(workflow.minimap_reference),
                "reference_sha256": sha256_file(path=workflow.minimap_reference),
                "index_sha256": sha256_file(path=output_index),
                "index_size_bytes": output_index.stat().st_size,
                **index_details,
                "allocated_memory_mb": workflow.memory_minimap2_mb,
            },
        )
    except Exception as error:
        if workflow.minimap_failure_policy == "fail":
            raise
        output_index.unlink(missing_ok=True)
        LOGGER.exception("minimap2 index preparation failed; continuing with partial reporting")
        _write_stage_failure(
            path=output_completion,
            stage="minimap2_index",
            error=error,
            status=_failure_status(error=error),
            extra={
                "reference": str(workflow.minimap_reference),
                "allocated_memory_mb": workflow.memory_minimap2_mb,
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
    _validate_deployment_identity(workflow=workflow)
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
    _validate_deployment_identity(workflow=workflow)
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
    """Classify selected samples and always publish a terminal stage record.

    The production DAG invokes every method with one sample per job so resource,
    timeout and failure state are isolated. A multi-sample call remains available
    for local fallback execution and preserves successful samples when a later
    sample fails. With the default ``continue`` policy, tool, resource, timeout
    and output-validation failures become explicit status records rather than
    failed reporting dependencies.
    """
    if method not in {"kraken2", "metabuli", "minimap2", "kmersutra"}:
        raise ValueError(f"Unsupported classifier method: {method}")
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    selected_samples = workflow.samples
    if sample_id is not None:
        selected_samples = tuple(
            sample for sample in workflow.samples if sample.sample_id == sample_id
        )
        if not selected_samples:
            raise ValueError(f"Unknown sample_id for classification: {sample_id}")
    failure_policy = _method_failure_policy(workflow=workflow, method=method)
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
    stage_root = stage_completion.parent
    stage_root.mkdir(parents=True, exist_ok=True)
    stage_log = (
        workflow.output_directory
        / "workflow_control"
        / "resource_staging_logs"
        / f"{method}.{sample_id or 'batch'}.log"
    )
    successes: list[str] = []
    failures: list[dict[str, str]] = []
    try:
        _require_tools(action=f"classify-{method}")
        resource, reference_fasta = _classifier_resources(workflow=workflow, method=method)
        _validate_classifier_resource(
            workflow=workflow,
            method=method,
            resource=resource,
            reference_fasta=reference_fasta,
        )
        validate_scratch(
            scratch_root=workflow.scratch_root,
            minimum_gb=workflow.minimum_scratch_gb,
        )
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
                    successes.append(sample.sample_id)
                except Exception as error:
                    if failure_policy == "fail":
                        raise
                    LOGGER.exception(
                        "%s failed for %s; continuing as configured",
                        method,
                        sample.sample_id,
                    )
                    status = _failure_status(error=error)
                    _write_method_failure(
                        workflow=workflow,
                        sample=sample,
                        method=method,
                        error=error,
                        status=status,
                    )
                    failures.append(
                        {
                            "sample_id": sample.sample_id,
                            "status": status,
                            "message": str(error),
                        }
                    )
    except Exception as error:
        if failure_policy == "fail":
            raise
        LOGGER.exception("%s stage failed; continuing to partial reporting", method)
        status = _failure_status(error=error)
        already_terminal = set(successes) | {row["sample_id"] for row in failures}
        for sample in selected_samples:
            if sample.sample_id in already_terminal:
                continue
            _write_method_failure(
                workflow=workflow,
                sample=sample,
                method=method,
                error=error,
                status=status,
            )
            failures.append(
                {
                    "sample_id": sample.sample_id,
                    "status": status,
                    "message": str(error),
                }
            )

    stage_status = _batch_status(successes=successes, failures=failures)
    write_json_atomic(
        path=stage_completion,
        payload={
            "status": stage_status,
            "completed_at_utc": utc_now(),
            "stage": method,
            "sample_count": len(selected_samples),
            "successful_samples": successes,
            "failed_samples": [row["sample_id"] for row in failures],
            "sample_failures": failures,
            "failure_policy": failure_policy,
            "allocated_memory_mb": _method_memory_mb(workflow=workflow, method=method),
            "allocated_runtime_minutes": _method_runtime_minutes(
                workflow=workflow,
                method=method,
            ),
        },
    )


def aggregate_results(*, config_path: Path, completion_path: Path) -> None:
    """Create tabular and HTML reports from every validated available result.

    Classifier outputs are discovered only after their per-sample terminal
    status is read. Missing or malformed method outputs become non-fatal report
    warnings. Consequently, final reporting is still produced when any subset
    of classifiers succeeds.
    """
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    truth_records = _load_optional_pcr_truth(workflow=workflow)
    truth_by_sample = {record.sample_id: record for record in truth_records}
    final_root = completion_path.parent
    final_root.mkdir(parents=True, exist_ok=True)
    sample_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, str]] = []
    taxon_rows: list[dict[str, Any]] = []
    minimap_rows: list[dict[str, Any]] = []
    kmersutra_rows: list[dict[str, Any]] = []
    warning_rows: list[dict[str, str]] = []
    for sample in workflow.samples:
        truth = truth_by_sample.get(sample.sample_id)
        try:
            host_summary = _read_single_tsv(
                path=_host_result_directory(workflow=workflow, sample=sample)
                / "host_removal_summary.tsv"
            )
            if host_summary.get("sample_id") != sample.sample_id:
                raise ValueError(
                    "Host summary sample_id does not match the manifest: "
                    f"{host_summary.get('sample_id')!r} != {sample.sample_id!r}"
                )
        except (OSError, ValueError) as error:
            host_summary = {
                "sample_id": sample.sample_id,
                "input_reads": "",
                "non_host_reads": "",
                "host_reads_removed": "",
                "host_fraction": "",
                "input_read_state": workflow.input_read_state,
                "host_depletion_performed": "",
            }
            warning_rows.append(
                _report_warning(
                    sample_id=sample.sample_id,
                    source="host_preparation",
                    message=f"Host/input summary is unavailable: {error}",
                )
            )
        sample_statuses: dict[str, str] = {}
        for method in METHODS:
            record = _method_status_record(
                workflow=workflow,
                sample=sample,
                method=method,
            )
            status_rows.append(record)
            sample_statuses[method] = record["status"]
        for method in ("kraken2", "metabuli"):
            if sample_statuses[method] != "success":
                continue
            report = (
                _method_result_directory(workflow=workflow, sample=sample, method=method)
                / "report.tsv"
            )
            try:
                taxon_rows.extend(
                    _parse_classifier_report(
                        path=report,
                        sample_id=sample.sample_id,
                        method=method,
                    )
                )
            except (OSError, ValueError) as error:
                _invalidate_report_status(
                    status_rows=status_rows,
                    sample_statuses=sample_statuses,
                    sample_id=sample.sample_id,
                    method=method,
                    error=error,
                )
                warning_rows.append(
                    _report_warning(
                        sample_id=sample.sample_id,
                        source=method,
                        message=f"Validated completion had an unreadable report: {error}",
                    )
                )
        if sample_statuses["minimap2"] == "success":
            try:
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
            except (OSError, ValueError) as error:
                _invalidate_report_status(
                    status_rows=status_rows,
                    sample_statuses=sample_statuses,
                    sample_id=sample.sample_id,
                    method="minimap2",
                    error=error,
                )
                warning_rows.append(
                    _report_warning(
                        sample_id=sample.sample_id,
                        source="minimap2",
                        message=f"Validated completion had an unreadable report: {error}",
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
        if sample_statuses["kmersutra"] == "success":
            try:
                kmersutra_rows.extend(
                    _read_rows_with_prefix(path=calls, sample_id=sample.sample_id)
                )
            except (OSError, ValueError) as error:
                _invalidate_report_status(
                    status_rows=status_rows,
                    sample_statuses=sample_statuses,
                    sample_id=sample.sample_id,
                    method="kmersutra",
                    error=error,
                )
                warning_rows.append(
                    _report_warning(
                        sample_id=sample.sample_id,
                        source="kmersutra",
                        message=f"Validated completion had unreadable calls: {error}",
                    )
                )
        sample_rows.append(
            {
                **host_summary,
                "input_read_state": workflow.input_read_state,
                "pcr_status": truth.pcr_status if truth is not None else "not_configured",
                "pcr_species": truth.pcr_species_source_text if truth is not None else "",
                "include_in_primary_comparison": (
                    str(truth.include_in_primary_comparison).lower()
                    if truth is not None
                    else "false"
                ),
                **{f"{method}_status": sample_statuses[method] for method in METHODS},
            }
        )

    summary_path = final_root / "sample_summary.tsv"
    status_path = final_root / "classifier_status.tsv"
    taxon_path = final_root / "classifier_taxon_reports.tsv.gz"
    calls_path = final_root / "kmersutra_species_calls.tsv.gz"
    minimap_path = final_root / "minimap2_taxon_reports.tsv.gz"
    evidence_path = final_root / "normalised_classifier_evidence.tsv.gz"
    warning_path = final_root / "report_warnings.tsv"
    report_data_path = final_root / "report_data.json"
    pcr_truth_path = final_root / "pcr_truth.tsv"
    pcr_concordance_path = final_root / "pcr_concordance.tsv"
    pcr_method_summary_path = final_root / "pcr_method_summary.tsv"
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
            "pcr_status",
            "pcr_species",
            "include_in_primary_comparison",
            "kraken2_status",
            "metabuli_status",
            "minimap2_status",
            "kmersutra_status",
        ),
    )
    _write_tsv(
        path=status_path,
        rows=status_rows,
        fieldnames=(
            "sample_id",
            "method",
            "status",
            "message",
            "completed_at_utc",
            "allocated_memory_mb",
            "allocated_runtime_minutes",
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
    calls_fields = tuple(dict.fromkeys(key for row in kmersutra_rows for key in row)) or (
        "sample_id",
        "method",
        "species_name",
        "detection_call",
    )
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
    evidence_rows = build_normalised_evidence(
        classifier_rows=taxon_rows,
        minimap_rows=minimap_rows,
        kmersutra_rows=kmersutra_rows,
        focus_taxa=workflow.report_focus_taxa,
    )
    concordance_rows, pcr_method_summary_rows = build_pcr_concordance(
        truth_records=truth_records,
        status_rows=status_rows,
        evidence_rows=evidence_rows,
    )
    _write_tsv(
        path=evidence_path,
        rows=evidence_rows,
        fieldnames=(
            "sample_id",
            "method",
            "taxon_name",
            "tax_id",
            "rank",
            "evidence_count",
            "supporting_count",
            "fraction",
            "metric",
            "detected",
            "comparable",
            "is_focus",
        ),
    )
    _write_tsv(
        path=warning_path,
        rows=warning_rows,
        fieldnames=("sample_id", "source", "message"),
    )
    _write_tsv(
        path=pcr_truth_path,
        rows=pcr_truth_rows(records=truth_records),
        fieldnames=(
            "sample_id",
            "pcr_status",
            "pcr_species",
            "pcr_species_canonical",
            "pcr_assay_or_source",
            "pcr_notes",
            "include_in_primary_comparison",
        ),
    )
    _write_tsv(
        path=pcr_concordance_path,
        rows=concordance_rows,
        fieldnames=(
            "sample_id",
            "method",
            "pcr_status",
            "pcr_species",
            "pcr_species_canonical",
            "include_in_primary_comparison",
            "classifier_status",
            "detected_plasmodium_species",
            "detected_expected_species",
            "missed_expected_species",
            "additional_plasmodium_species",
            "expected_species_count",
            "detected_expected_species_count",
            "expected_species_evidence_count",
            "all_expected_species_detected",
            "species_exact_match",
            "comparison_status",
        ),
    )
    _write_tsv(
        path=pcr_method_summary_path,
        rows=pcr_method_summary_rows,
        fieldnames=(
            "method",
            "primary_sample_count",
            "available_sample_count",
            "unavailable_sample_count",
            "pcr_positive_available_count",
            "all_expected_species_detected_count",
            "exact_species_match_count",
            "pcr_negative_available_count",
            "concordant_negative_count",
        ),
    )
    write_json_atomic(
        path=report_data_path,
        payload=serialisable_report_data(
            sample_rows=sample_rows,
            status_rows=status_rows,
            evidence_rows=evidence_rows,
            warning_rows=warning_rows,
            pcr_truth_rows=pcr_truth_rows(records=truth_records),
            pcr_concordance_rows=concordance_rows,
            pcr_method_summary_rows=pcr_method_summary_rows,
        ),
    )
    html_paths = generate_html_reports(
        workflow=workflow,
        sample_rows=sample_rows,
        status_rows=status_rows,
        evidence_rows=evidence_rows,
        warning_rows=warning_rows,
        pcr_concordance_rows=concordance_rows,
        pcr_method_summary_rows=pcr_method_summary_rows,
        final_root=final_root,
    )
    report_manifest_path = final_root / "report_manifest.json"
    readme_path = final_root / "README.txt"
    pcr_readme = (
        "pcr_truth.tsv preserves the independent PCR interpretation used for evaluation.\n"
        "pcr_concordance.tsv compares each method with PCR without converting failures "
        "into biological non-detections.\n"
        "pcr_method_summary.tsv reports exact counts and denominators by method.\n"
        if truth_records
        else "PCR evaluation was not configured; PCR TSV files contain schema headers only.\n"
    )
    readme_path.write_text(
        "Nanopore real-data workflow final results\n\n"
        "Open reports/index.html for the offline final report.\n"
        "reports/classifiers/ contains one detailed HTML report per classifier.\n"
        "reports/comparison.html compares available method-specific evidence.\n"
        "reports/samples/ contains one navigable report per sample.\n"
        "sample_summary.tsv contains host-removal totals and classifier status.\n"
        "classifier_status.tsv records every successful, failed, timed-out, disabled, "
        "unavailable or missing sample-method result.\n"
        "classifier_taxon_reports.tsv.gz harmonises Kraken2 and Metabuli reports.\n"
        "minimap2_taxon_reports.tsv.gz summarises controlled-reference best alignments.\n"
        "kmersutra_species_calls.tsv.gz combines KmerSutra species calls.\n"
        "normalised_classifier_evidence.tsv.gz supports descriptive comparison without "
        "forcing consensus.\n"
        f"{pcr_readme}"
        "report_data.json is a reusable machine-readable reporting payload.\n"
        "report_warnings.tsv records non-fatal parsing or completeness problems.\n"
        "SHA256SUMS.tsv records final-file checksums.\n",
        encoding="utf-8",
    )
    checksum_path = final_root / "SHA256SUMS.tsv"
    checksum_targets = [
        summary_path,
        status_path,
        taxon_path,
        minimap_path,
        calls_path,
        evidence_path,
        warning_path,
        pcr_truth_path,
        pcr_concordance_path,
        pcr_method_summary_path,
        report_data_path,
        report_manifest_path,
        readme_path,
        *html_paths,
    ]
    checksum_rows = [
        {
            "sha256": sha256_file(path=path),
            "file": str(path.relative_to(final_root)),
        }
        for path in checksum_targets
    ]
    _write_tsv(
        path=checksum_path,
        rows=checksum_rows,
        fieldnames=("sha256", "file"),
    )
    enabled_statuses = [
        row["status"]
        for row in status_rows
        if not (row["method"] == "kmersutra" and not workflow.kmersutra_enabled)
    ]
    successful = sum(status == "success" for status in enabled_statuses)
    if successful == len(enabled_statuses):
        final_status = "success"
    elif successful:
        final_status = "partial"
    else:
        final_status = "failed"
    write_json_atomic(
        path=completion_path,
        payload={
            "status": final_status,
            "reporting_status": "success",
            "completed_at_utc": utc_now(),
            "run_id": workflow.run_id,
            "sample_count": len(workflow.samples),
            "pcr_evaluation_configured": workflow.pcr_truth_path is not None,
            "successful_classifier_runs": successful,
            "expected_classifier_runs": len(enabled_statuses),
            "reports_generated": True,
            "final_report": str(final_root / "reports" / "index.html"),
            "outputs": [
                str(path)
                for path in (
                    summary_path,
                    status_path,
                    taxon_path,
                    minimap_path,
                    calls_path,
                    evidence_path,
                    warning_path,
                    pcr_truth_path,
                    pcr_concordance_path,
                    pcr_method_summary_path,
                    report_data_path,
                    report_manifest_path,
                    readme_path,
                    checksum_path,
                    *html_paths,
                )
            ],
        },
    )


def _report_warning(*, sample_id: str, source: str, message: str) -> dict[str, str]:
    """Create a consistent non-fatal reporting warning row."""
    return {"sample_id": sample_id, "source": source, "message": message}


def _invalidate_report_status(
    *,
    status_rows: list[dict[str, str]],
    sample_statuses: dict[str, str],
    sample_id: str,
    method: str,
    error: Exception,
) -> None:
    """Prevent an unreadable successful table becoming a non-detection."""
    sample_statuses[method] = "invalid"
    matching = [
        row for row in status_rows if row["sample_id"] == sample_id and row["method"] == method
    ]
    if len(matching) != 1:
        raise RuntimeError(
            f"Expected one status row for {sample_id}/{method}; found {len(matching)}"
        )
    matching[0]["status"] = "invalid"
    matching[0]["message"] = f"Completed output could not be parsed: {error}"


def record_scheduler_failure(
    *,
    config_path: Path,
    method: str,
    sample_id: str,
    message: str,
    status: str = "scheduler_failed",
) -> None:
    """Write a terminal classifier record after an outer Slurm failure.

    Args:
        config_path: Workflow YAML path.
        method: Classifier branch.
        sample_id: Logical sample identifier.
        message: Scheduler or wrapper failure detail.
        status: Machine-readable terminal status.
    """
    if method not in METHODS:
        raise ValueError(f"Unsupported classifier method: {method}")
    if not message.strip():
        raise ValueError("Scheduler failure message must not be blank")
    workflow = load_workflow_config(config_path=config_path)
    _validate_deployment_identity(workflow=workflow)
    matching = [sample for sample in workflow.samples if sample.sample_id == sample_id]
    if len(matching) != 1:
        raise ValueError(f"Unknown sample_id for scheduler failure: {sample_id}")
    error = RuntimeError(message.strip())
    _write_method_failure(
        workflow=workflow,
        sample=matching[0],
        method=method,
        error=error,
        status=status,
    )
    stage = (
        _method_result_directory(
            workflow=workflow,
            sample=matching[0],
            method=method,
        )
        / "stage.complete.json"
    )
    write_json_atomic(
        path=stage,
        payload={
            "status": status,
            "completed_at_utc": utc_now(),
            "stage": method,
            "sample_count": 1,
            "successful_samples": [],
            "failed_samples": [sample_id],
            "sample_failures": [
                {"sample_id": sample_id, "status": status, "message": message.strip()}
            ],
            "failure_policy": _method_failure_policy(
                workflow=workflow,
                method=method,
            ),
            "allocated_memory_mb": _method_memory_mb(
                workflow=workflow,
                method=method,
            ),
            "allocated_runtime_minutes": _method_runtime_minutes(
                workflow=workflow,
                method=method,
            ),
        },
    )


def _validate_deployment_identity(*, workflow: WorkflowConfig) -> None:
    """Reject mixed repository copies or an unexpected package version."""
    actual_root = Path(__file__).resolve().parents[2]
    expected_root = workflow.expected_repository_root.resolve()
    if actual_root != expected_root:
        raise RuntimeError(
            "Workflow source root does not match deployment.expected_repository_root; "
            f"actual={actual_root}, expected={expected_root}"
        )
    if __version__ != workflow.expected_package_version:
        raise RuntimeError(
            "Workflow package version does not match deployment.expected_package_version; "
            f"actual={__version__}, expected={workflow.expected_package_version}"
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
        timeout_seconds=_internal_timeout_seconds(runtime_minutes=workflow.runtime_kraken2_minutes),
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
        timeout_seconds=_internal_timeout_seconds(
            runtime_minutes=workflow.runtime_metabuli_minutes
        ),
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
    input_read_count = _validated_non_host_read_count(
        workflow=workflow,
        sample=sample,
    )
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
        timeout_seconds=_internal_timeout_seconds(
            runtime_minutes=workflow.runtime_minimap2_minutes
        ),
    )
    summarise_minimap_paf(
        paf_path=paf_path,
        reference_fasta=reference_fasta,
        taxon_report_path=output_directory / "taxon_report.tsv",
        mapping_summary_path=output_directory / "mapping_summary.tsv",
        sample_id=sample.sample_id,
        minimum_mapq=workflow.minimap_min_mapq,
        minimum_alignment=workflow.minimap_min_alignment,
        input_read_count=input_read_count,
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
    _read_rows_with_prefix(path=calls, sample_id=sample.sample_id)


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


def _classifier_resources(
    *,
    workflow: WorkflowConfig,
    method: str,
) -> tuple[Path, Path | None]:
    """Return the primary resource and optional reference for a classifier."""
    if method == "kraken2":
        return workflow.kraken_database, None
    if method == "metabuli":
        return workflow.metabuli_database, None
    if method == "minimap2":
        return _resolved_minimap_index(workflow=workflow), workflow.minimap_reference
    if method == "kmersutra":
        return _required_kmersutra_panel(workflow=workflow), None
    raise ValueError(f"Unsupported classifier method: {method}")


def _validate_classifier_resource(
    *,
    workflow: WorkflowConfig,
    method: str,
    resource: Path,
    reference_fasta: Path | None,
) -> None:
    """Validate one classifier's resource without inspecting unrelated tools."""
    del workflow
    if method == "kraken2":
        _validate_kraken_database(database=resource)
    elif method == "metabuli":
        _validate_non_empty_directory(label="Metabuli database", path=resource)
    elif method == "minimap2":
        if not resource.is_file() or resource.stat().st_size == 0:
            raise FileNotFoundError(f"minimap2 index is missing or empty: {resource}")
        if reference_fasta is None or not reference_fasta.is_file():
            raise FileNotFoundError(
                f"minimap2 classification reference is missing: {reference_fasta}"
            )
    elif method == "kmersutra":
        _validate_panel(panel=resource)
    else:
        raise ValueError(f"Unsupported classifier method: {method}")


def _classifier_readiness(
    *,
    workflow: WorkflowConfig,
    method: str,
    required_species: Sequence[str] = (),
) -> dict[str, str]:
    """Return a non-raising preflight status for one independent classifier."""
    if method == "kmersutra" and not workflow.kmersutra_enabled:
        return {
            "method": method,
            "status": "disabled",
            "message": "Disabled by configuration",
            "resource": "",
        }
    try:
        _require_tools(action=f"classify-{method}")
        resource, reference = _classifier_resources(workflow=workflow, method=method)
        if method == "minimap2" and workflow.minimap_index is None:
            if workflow.minimap_genome_config is not None:
                resource_text = (
                    "controlled reference and index will be built from "
                    f"{workflow.minimap_genome_config}"
                )
            elif not workflow.minimap_reference.is_file():
                raise FileNotFoundError(
                    f"minimap2 classification reference is missing: {workflow.minimap_reference}"
                )
            else:
                validate_required_reference_species(
                    reference_fasta=workflow.minimap_reference,
                    required_species=required_species,
                )
                resource_text = "index will be built from configured reference"
        else:
            _validate_classifier_resource(
                workflow=workflow,
                method=method,
                resource=resource,
                reference_fasta=reference,
            )
            resource_text = str(resource)
        return {
            "method": method,
            "status": "ready",
            "message": "Executable and configured resource are available",
            "resource": resource_text,
        }
    except Exception as error:
        return {
            "method": method,
            "status": "unavailable",
            "message": str(error),
            "resource": str(_classifier_resource_hint(workflow=workflow, method=method)),
        }


def _classifier_resource_hint(*, workflow: WorkflowConfig, method: str) -> Path | str:
    """Return a resource label without requiring the resource to be valid."""
    if method == "kraken2":
        return workflow.kraken_database
    if method == "metabuli":
        return workflow.metabuli_database
    if method == "minimap2":
        return workflow.minimap_index or workflow.minimap_reference
    return workflow.kmersutra_panel or "not configured"


def _method_failure_policy(*, workflow: WorkflowConfig, method: str) -> str:
    policies = {
        "kraken2": workflow.kraken_failure_policy,
        "metabuli": workflow.metabuli_failure_policy,
        "minimap2": workflow.minimap_failure_policy,
        "kmersutra": workflow.kmersutra_failure_policy,
    }
    try:
        return policies[method]
    except KeyError as error:
        raise ValueError(f"Unsupported classifier method: {method}") from error


def _method_memory_mb(*, workflow: WorkflowConfig, method: str) -> int:
    values = {
        "kraken2": workflow.memory_kraken2_mb,
        "metabuli": workflow.memory_metabuli_mb,
        "minimap2": workflow.memory_minimap2_mb,
        "kmersutra": workflow.memory_kmersutra_mb,
    }
    return values[method]


def _method_runtime_minutes(*, workflow: WorkflowConfig, method: str) -> int:
    values = {
        "kraken2": workflow.runtime_kraken2_minutes,
        "metabuli": workflow.runtime_metabuli_minutes,
        "minimap2": workflow.runtime_minimap2_minutes,
        "kmersutra": workflow.runtime_kmersutra_minutes,
    }
    return values[method]


def _internal_timeout_seconds(*, runtime_minutes: int) -> int:
    """Leave ten minutes for failure publication before scheduler termination."""
    if runtime_minutes <= 0:
        raise ValueError("Runtime minutes must be positive")
    return max(60, (runtime_minutes - 10) * 60)


def _failure_status(*, error: Exception) -> str:
    """Classify a caught exception for machine-readable reporting."""
    if isinstance(error, CommandTimeoutError) or "time limit" in str(error).casefold():
        return "timeout"
    if isinstance(error, (FileNotFoundError, PermissionError)):
        return "unavailable"
    return "failed"


def _batch_status(
    *,
    successes: Sequence[str],
    failures: Sequence[Mapping[str, str]],
) -> str:
    """Summarise a classifier batch without discarding failure specificity."""
    if successes and failures:
        return "partial"
    if successes:
        return "success"
    statuses = {str(row.get("status", "failed")) for row in failures}
    if len(statuses) == 1:
        return next(iter(statuses))
    return "failed"


def _write_method_failure(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    method: str,
    error: Exception,
    status: str,
) -> None:
    """Write a compact per-sample failure record and invalidate stale success."""
    result = _method_result_directory(workflow=workflow, sample=sample, method=method)
    result.mkdir(parents=True, exist_ok=True)
    (result / "complete.json").unlink(missing_ok=True)
    write_json_atomic(
        path=result / "failure.json",
        payload={
            "status": status,
            "sample_id": sample.sample_id,
            "method": method,
            "error": str(error),
            "error_type": type(error).__name__,
            "failed_at_utc": utc_now(),
            "allocated_memory_mb": _method_memory_mb(workflow=workflow, method=method),
            "allocated_runtime_minutes": _method_runtime_minutes(
                workflow=workflow,
                method=method,
            ),
        },
    )


def _write_stage_failure(
    *,
    path: Path,
    stage: str,
    error: Exception,
    status: str,
    extra: Mapping[str, Any] | None = None,
) -> None:
    """Write a terminal stage token after a recoverable infrastructure failure."""
    payload: dict[str, Any] = {
        "status": status,
        "stage": stage,
        "error": str(error),
        "error_type": type(error).__name__,
        "completed_at_utc": utc_now(),
    }
    if extra:
        payload.update(extra)
    write_json_atomic(path=path, payload=payload)


def _safe_resource_fingerprint(*, path: Path, checksum_file: bool) -> str:
    """Fingerprint an available resource or return an explicit missing label."""
    if not path.exists():
        return f"missing:{path}"
    return metadata_fingerprint(paths=[path], checksum_files=checksum_file)


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
    result = _method_result_directory(
        workflow=workflow,
        sample=sample,
        method=method,
    )
    failure = result / "failure.json"
    if failure.is_file():
        try:
            return str(json.loads(failure.read_text(encoding="utf-8")).get("status", "failed"))
        except (json.JSONDecodeError, OSError):
            return "invalid"
    path = result / "complete.json"
    if not path.is_file():
        return "missing"
    try:
        return str(json.loads(path.read_text(encoding="utf-8")).get("status", "invalid"))
    except (json.JSONDecodeError, OSError):
        return "invalid"


def _method_status_record(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
    method: str,
) -> dict[str, str]:
    """Return a human- and machine-readable classifier terminal status."""
    status = _method_status(workflow=workflow, sample=sample, method=method)
    result = _method_result_directory(workflow=workflow, sample=sample, method=method)
    metadata_path = result / ("complete.json" if status == "success" else "failure.json")
    payload: dict[str, Any] = {}
    if metadata_path.is_file():
        try:
            value = json.loads(metadata_path.read_text(encoding="utf-8"))
            if isinstance(value, dict):
                payload = value
        except (json.JSONDecodeError, OSError):
            payload = {}
    message = ""
    if status == "disabled":
        message = "Disabled by configuration"
    elif status == "success":
        message = "Validated outputs available"
    elif status == "missing":
        message = "No terminal result record was found"
    else:
        message = str(payload.get("error", payload.get("reason", status)))
    return {
        "sample_id": sample.sample_id,
        "method": method,
        "status": status,
        "message": message,
        "completed_at_utc": str(payload.get("completed_at_utc", payload.get("failed_at_utc", ""))),
        "allocated_memory_mb": str(
            payload.get("allocated_memory_mb", _method_memory_mb(workflow=workflow, method=method))
        ),
        "allocated_runtime_minutes": str(
            payload.get(
                "allocated_runtime_minutes",
                _method_runtime_minutes(workflow=workflow, method=method),
            )
        ),
    }


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
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_name", "call"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValueError(f"KmerSutra calls table is missing columns: {missing}")
        rows = list(reader)
    for row in rows:
        observed_sample = (row.get("sample_id") or "").strip()
        if observed_sample and observed_sample != sample_id:
            raise ValueError(
                "KmerSutra calls sample_id does not match the manifest: "
                f"{observed_sample!r} != {sample_id!r}"
            )
    return [{**row, "sample_id": sample_id, "method": "kmersutra"} for row in rows]


def _read_rows_with_method(
    *,
    path: Path,
    sample_id: str,
    method: str,
) -> list[dict[str, Any]]:
    """Read a TSV and enforce its workflow sample and method fields."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "tax_id",
            "taxon_name",
            "best_read_count",
            "ambiguous_best_read_count",
            "alignment_count",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValueError(f"{method} taxon table is missing columns: {missing}")
        rows = list(reader)
    return [{**row, "sample_id": sample_id, "method": method} for row in rows]


def _read_single_tsv(*, path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one data row in {path}")
    return rows[0]


def _validated_non_host_read_count(
    *,
    workflow: WorkflowConfig,
    sample: Sample,
) -> int:
    """Return the validated read count used to bound minimap2 output.

    Args:
        workflow: Validated workflow configuration.
        sample: Logical sample.

    Returns:
        Non-negative number of reads supplied to every classifier.

    Raises:
        ValueError: If the host-removal summary is missing or invalid.
    """
    summary = _read_single_tsv(
        path=_host_result_directory(workflow=workflow, sample=sample) / "host_removal_summary.tsv"
    )
    if summary.get("sample_id") != sample.sample_id:
        raise ValueError(
            "Host summary sample_id does not match the requested sample: "
            f"{summary.get('sample_id')!r} != {sample.sample_id!r}"
        )
    value = summary.get("non_host_reads", "")
    try:
        count = int(value)
    except ValueError as error:
        raise ValueError(f"Invalid non_host_reads for {sample.sample_id}: {value!r}") from error
    if count < 0:
        raise ValueError(f"Negative non_host_reads for {sample.sample_id}: {count}")
    return count


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
