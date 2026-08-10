"""Detached Slurm submission for independent, restartable classifier arrays."""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from nanopore_realdata import __version__
from nanopore_realdata.config import WorkflowConfig, load_workflow_config
from nanopore_realdata.runtime import sha256_file, utc_now, write_json_atomic


JOB_ID_PATTERN = re.compile(r"^[0-9]+$")
METHODS = ("kraken2", "metabuli", "minimap2", "kmersutra")


@dataclass(frozen=True)
class JobResources:
    """Slurm resources for one submitted stage or array task."""

    cpus: int
    memory_mb: int
    runtime_minutes: int
    qos: str = ""


@dataclass(frozen=True)
class JobPlan:
    """One dependency-aware Slurm job specification."""

    key: str
    name: str
    action: str
    resources: JobResources
    dependency_mode: str = ""
    dependency_keys: tuple[str, ...] = ()
    method: str = ""
    array: str = ""


def build_submission_plan(*, workflow: WorkflowConfig) -> tuple[JobPlan, ...]:
    """Build the detached DAG used for the PCR-validated real-read run.

    Args:
        workflow: Validated workflow configuration.

    Returns:
        Ordered job specifications with symbolic dependencies.

    Raises:
        ValueError: If the current detached backend cannot prepare the inputs.
    """
    if workflow.input_read_state != "host_removed":
        raise ValueError(
            "Detached Slurm submission currently requires inputs.read_state=host_removed"
        )
    sample_max = len(workflow.samples) - 1
    if sample_max < 0:
        raise ValueError("At least one sample is required")
    default_qos = workflow.slurm_default_qos
    plans: list[JobPlan] = [
        JobPlan(
            key="preflight",
            name="NRD_preflight",
            action="preflight",
            resources=JobResources(2, 4096, 60, default_qos),
        ),
        JobPlan(
            key="accept_host_removed",
            name="NRD_accept",
            action="accept-host-removed",
            resources=JobResources(
                workflow.threads_host,
                workflow.memory_host_mb,
                workflow.runtime_host_minutes,
                default_qos,
            ),
            dependency_mode="afterany",
            dependency_keys=("preflight",),
        ),
    ]
    index_dependency = "preflight"
    if workflow.minimap_genome_config is not None:
        plans.append(
            JobPlan(
                key="build_minimap_reference",
                name="NRD_reference",
                action="build-minimap-reference",
                resources=JobResources(
                    workflow.threads_minimap2,
                    workflow.memory_minimap2_mb,
                    workflow.runtime_minimap2_minutes,
                    default_qos,
                ),
                dependency_mode="afterany",
                dependency_keys=("preflight",),
            )
        )
        index_dependency = "build_minimap_reference"
    plans.append(
        JobPlan(
            key="build_minimap_index",
            name="NRD_index",
            action="build-minimap-index",
            resources=JobResources(
                workflow.threads_minimap2,
                workflow.memory_minimap2_mb,
                workflow.runtime_minimap2_minutes,
                default_qos,
            ),
            dependency_mode="afterany",
            dependency_keys=(index_dependency,),
        )
    )
    method_resources = {
        "kraken2": JobResources(
            workflow.threads_kraken2,
            workflow.memory_kraken2_mb,
            workflow.runtime_kraken2_minutes,
            default_qos,
        ),
        "metabuli": JobResources(
            workflow.threads_metabuli,
            workflow.memory_metabuli_mb,
            workflow.runtime_metabuli_minutes,
            default_qos,
        ),
        "minimap2": JobResources(
            workflow.threads_minimap2,
            workflow.memory_minimap2_mb,
            workflow.runtime_minimap2_minutes,
            default_qos,
        ),
        "kmersutra": JobResources(
            workflow.threads_kmersutra,
            workflow.memory_kmersutra_mb,
            workflow.runtime_kmersutra_minutes,
            workflow.slurm_kmersutra_qos,
        ),
    }
    concurrency = {
        "kraken2": workflow.slurm_kraken2_concurrency,
        "metabuli": workflow.slurm_metabuli_concurrency,
        "minimap2": workflow.slurm_minimap2_concurrency,
        "kmersutra": workflow.slurm_kmersutra_concurrency,
    }
    classifier_keys: list[str] = []
    for method in METHODS:
        if method == "kmersutra" and not workflow.kmersutra_enabled:
            continue
        dependencies = ["accept_host_removed"]
        if method == "minimap2":
            dependencies.append("build_minimap_index")
        key = f"classify_{method}"
        classifier_keys.append(key)
        plans.append(
            JobPlan(
                key=key,
                name=f"NRD_{method[:8]}",
                action="classify",
                method=method,
                resources=method_resources[method],
                dependency_mode="afterany",
                dependency_keys=tuple(dependencies),
                array=f"0-{sample_max}%{concurrency[method]}",
            )
        )
    plans.append(
        JobPlan(
            key="aggregate",
            name="NRD_aggregate",
            action="aggregate",
            resources=JobResources(2, 8192, 180, default_qos),
            dependency_mode="afterany",
            dependency_keys=tuple(classifier_keys),
        )
    )
    return tuple(plans)


def _validated_repository_root(*, workflow: WorkflowConfig) -> Path:
    """Validate the executing checkout and configured package version.

    Args:
        workflow: Validated workflow configuration.

    Returns:
        Resolved root of the checkout providing the executing Python module.

    Raises:
        RuntimeError: If the executing checkout or package version differs from
            the frozen deployment configuration.
    """
    repository_root = Path(__file__).resolve().parents[2]
    expected_root = workflow.expected_repository_root.resolve()
    if repository_root != expected_root:
        raise RuntimeError(
            "Refusing Slurm execution from an unexpected repository copy: "
            f"actual={repository_root}; expected={expected_root}; "
            f"config={workflow.config_path}"
        )
    if __version__ != workflow.expected_package_version:
        raise RuntimeError(
            f"Refusing Slurm execution from version {__version__}; "
            f"expected {workflow.expected_package_version}"
        )
    return repository_root


def submit_workflow(
    *,
    config_path: Path,
    resume_submission: bool,
    new_attempt: bool,
) -> Path:
    """Submit the detached Slurm DAG and journal every accepted job ID.

    Args:
        config_path: Workflow configuration.
        resume_submission: Continue an interrupted submission journal.
        new_attempt: Archive a completed journal and submit a fresh attempt.

    Returns:
        Durable submission journal path.

    Raises:
        RuntimeError: If submission is unsafe or ``sbatch`` rejects a job.
        ValueError: If mutually exclusive controls are requested.
    """
    if resume_submission and new_attempt:
        raise ValueError("resume_submission and new_attempt are mutually exclusive")
    workflow = load_workflow_config(config_path=config_path)
    repository_root = _validated_repository_root(workflow=workflow)
    stage_script = repository_root / "workflow" / "slurm" / "run_stage.sh"
    if not stage_script.is_file():
        raise FileNotFoundError(f"Slurm stage wrapper is missing: {stage_script}")
    control_root = workflow.output_directory / "workflow_control"
    log_root = control_root / "slurm_logs"
    control_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)
    journal_path = control_root / "slurm_submission.json"
    config_digest = sha256_file(path=workflow.config_path)
    journal = _initial_journal(
        workflow=workflow,
        repository_root=repository_root,
        config_digest=config_digest,
    )
    if journal_path.is_file():
        current = _load_journal(path=journal_path)
        _validate_journal_identity(
            journal=current,
            workflow=workflow,
            repository_root=repository_root,
            config_digest=config_digest,
        )
        if resume_submission:
            if current.get("status") != "submitting":
                raise RuntimeError(
                    "--resume-submission is only valid for an interrupted submission"
                )
            journal = current
        elif new_attempt:
            _refuse_active_jobs(journal=current)
            history = control_root / "submission_history"
            history.mkdir(parents=True, exist_ok=True)
            archived = history / f"slurm_submission.{utc_now().replace(':', '')}.json"
            shutil.copy2(journal_path, archived)
        else:
            raise RuntimeError(
                "A submission journal already exists. Use --resume-submission only for "
                "an interrupted submission, or --new-attempt after confirming prior jobs ended."
            )
    write_json_atomic(path=journal_path, payload=journal)
    plans = build_submission_plan(workflow=workflow)
    jobs = journal.setdefault("jobs", {})
    assert isinstance(jobs, dict)
    for plan in plans:
        if plan.key in jobs:
            continue
        dependencies = _dependency_ids(plan=plan, jobs=jobs)
        command = _sbatch_command(
            workflow=workflow,
            plan=plan,
            dependency_ids=dependencies,
            stage_script=stage_script,
            repository_root=repository_root,
            log_root=log_root,
        )
        job_id = _submit(command=command)
        jobs[plan.key] = {
            "job_id": job_id,
            "job_name": plan.name,
            "action": plan.action,
            "method": plan.method,
            "array": plan.array,
            "dependency_mode": plan.dependency_mode,
            "dependency_job_ids": dependencies,
            "submitted_at_utc": utc_now(),
            "command": shlex.join(command),
        }
        write_json_atomic(path=journal_path, payload=journal)
    journal["status"] = "submitted"
    journal["completed_submission_at_utc"] = utc_now()
    write_json_atomic(path=journal_path, payload=journal)
    return journal_path


def planned_commands(*, config_path: Path) -> list[dict[str, Any]]:
    """Return a non-mutating human-readable detached submission plan.

    Args:
        config_path: Workflow configuration.

    Returns:
        Ordered symbolic job descriptions.
    """
    workflow = load_workflow_config(config_path=config_path)
    _validated_repository_root(workflow=workflow)
    return [
        {
            "key": plan.key,
            "action": plan.action,
            "method": plan.method,
            "array": plan.array,
            "dependency_mode": plan.dependency_mode,
            "dependency_keys": list(plan.dependency_keys),
            "cpus": plan.resources.cpus,
            "memory_mb": plan.resources.memory_mb,
            "runtime_minutes": plan.resources.runtime_minutes,
            "qos": plan.resources.qos or "default",
        }
        for plan in build_submission_plan(workflow=workflow)
    ]


def _initial_journal(
    *,
    workflow: WorkflowConfig,
    repository_root: Path,
    config_digest: str,
) -> dict[str, Any]:
    return {
        "status": "submitting",
        "created_at_utc": utc_now(),
        "run_id": workflow.run_id,
        "config": str(workflow.config_path),
        "config_sha256": config_digest,
        "repository_root": str(repository_root),
        "package_version": __version__,
        "account": workflow.slurm_account,
        "partition": workflow.slurm_partition,
        "sample_count": len(workflow.samples),
        "jobs": {},
    }


def _load_journal(*, path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Invalid Slurm submission journal: {path}") from error
    if not isinstance(payload, dict) or not isinstance(payload.get("jobs"), dict):
        raise RuntimeError(f"Malformed Slurm submission journal: {path}")
    return payload


def _validate_journal_identity(
    *,
    journal: Mapping[str, Any],
    workflow: WorkflowConfig,
    repository_root: Path,
    config_digest: str,
) -> None:
    expected = {
        "run_id": workflow.run_id,
        "config": str(workflow.config_path),
        "config_sha256": config_digest,
        "repository_root": str(repository_root),
        "package_version": __version__,
    }
    differences = {
        key: {"journal": journal.get(key), "current": value}
        for key, value in expected.items()
        if journal.get(key) != value
    }
    if differences:
        raise RuntimeError(f"Slurm journal identity differs from this run: {differences}")


def _dependency_ids(*, plan: JobPlan, jobs: Mapping[str, Any]) -> tuple[str, ...]:
    identifiers: list[str] = []
    for key in plan.dependency_keys:
        record = jobs.get(key)
        if not isinstance(record, Mapping):
            raise RuntimeError(f"Dependency has not been submitted: {key}")
        job_id = str(record.get("job_id", ""))
        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise RuntimeError(f"Dependency has an invalid job ID: {key}={job_id!r}")
        identifiers.append(job_id)
    return tuple(identifiers)


def _sbatch_command(
    *,
    workflow: WorkflowConfig,
    plan: JobPlan,
    dependency_ids: Sequence[str],
    stage_script: Path,
    repository_root: Path,
    log_root: Path,
) -> list[str]:
    command = [
        "sbatch",
        "--parsable",
        f"--account={workflow.slurm_account}",
        f"--partition={workflow.slurm_partition}",
        f"--job-name={plan.name}",
        f"--cpus-per-task={plan.resources.cpus}",
        f"--mem={plan.resources.memory_mb}M",
        f"--time={plan.resources.runtime_minutes}",
        "--signal=B:TERM@300",
        f"--output={log_root}/%x.%A_%a.out",
        f"--error={log_root}/%x.%A_%a.err",
    ]
    if plan.resources.qos:
        command.append(f"--qos={plan.resources.qos}")
    if plan.array:
        command.append(f"--array={plan.array}")
    if dependency_ids:
        command.append(f"--dependency={plan.dependency_mode}:{':'.join(dependency_ids)}")
    command.extend(
        [
            str(stage_script),
            "--repository-root",
            str(repository_root),
            "--conda-environment",
            workflow.conda_environment,
            "--action",
            plan.action,
            "--config",
            str(workflow.config_path),
        ]
    )
    if plan.method:
        command.extend(["--method", plan.method, "--sample-index-from-slurm"])
    return command


def _submit(*, command: Sequence[str]) -> str:
    completed = subprocess.run(
        list(command),
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    output = completed.stdout.strip()
    if completed.returncode != 0:
        raise RuntimeError(
            f"sbatch rejected a workflow job with exit code {completed.returncode}: {output}"
        )
    job_id = output.split(";", maxsplit=1)[0].strip()
    if not JOB_ID_PATTERN.fullmatch(job_id):
        raise RuntimeError(f"Unexpected sbatch response: {output!r}")
    return job_id


def _refuse_active_jobs(*, journal: Mapping[str, Any]) -> None:
    jobs = journal.get("jobs", {})
    if not isinstance(jobs, Mapping):
        raise RuntimeError("Existing Slurm journal has no valid jobs mapping")
    known = {
        str(record.get("job_id", ""))
        for record in jobs.values()
        if isinstance(record, Mapping) and JOB_ID_PATTERN.fullmatch(str(record.get("job_id", "")))
    }
    if not known:
        return
    user = os.environ.get("USER", "").strip()
    if not user:
        raise RuntimeError("USER is unavailable; cannot check for active Slurm jobs")
    completed = subprocess.run(
        ["squeue", "--user", user, "--noheader", "--format=%A"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Could not verify that the previous Slurm attempt ended: " + completed.stdout.strip()
        )
    active = sorted(known.intersection(completed.stdout.split()))
    if active:
        raise RuntimeError(
            "Previous workflow jobs are still queued or running: " + ", ".join(active)
        )
