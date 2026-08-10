"""Command-line interface for the real Nanopore classifier workflow."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Sequence

from nanopore_realdata.config import WorkflowConfig, load_workflow_config
from nanopore_realdata.workflow import (
    accept_host_removed_batch,
    aggregate_results,
    build_host_index,
    build_minimap_index,
    build_minimap_reference,
    classify_batch,
    host_deplete_batch,
    preflight,
    record_scheduler_failure,
    run_snakemake,
)
from nanopore_realdata.slurm import planned_commands, submit_workflow


ACTIONS = (
    "validate",
    "preflight",
    "build-host-index",
    "build-minimap-reference",
    "build-minimap-index",
    "accept-host-removed",
    "host-deplete",
    "classify",
    "aggregate",
    "report",
    "record-scheduler-failure",
    "plan-slurm",
    "submit-slurm",
    "run",
)


def build_parser() -> argparse.ArgumentParser:
    """Build the named-argument-only command parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Prepare real Nanopore FASTQs, classify them independently with "
            "Kraken2, Metabuli, minimap2 and KmerSutra, and build failure-aware "
            "offline reports."
        )
    )
    parser.add_argument("--action", required=True, choices=ACTIONS)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--completion", type=Path)
    parser.add_argument("--host-index", type=Path)
    parser.add_argument(
        "--method",
        choices=("kraken2", "metabuli", "minimap2", "kmersutra"),
    )
    parser.add_argument("--sample-id")
    parser.add_argument("--sample-index", type=int)
    parser.add_argument("--reference-manifest", type=Path)
    parser.add_argument("--message")
    parser.add_argument("--status", default="scheduler_failed")
    parser.add_argument("--snakefile", type=Path)
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--cores", default="all")
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--target", action="append", default=[])
    parser.add_argument(
        "--set-resource",
        action="append",
        default=[],
        help="Snakemake RULE:RESOURCE=VALUE override; repeatable",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--unlock", action="store_true")
    parser.add_argument("--resume-submission", action="store_true")
    parser.add_argument("--new-attempt", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(arguments: Sequence[str] | None = None) -> int:
    """Execute one validated workflow action.

    Args:
        arguments: Optional argument vector for tests.

    Returns:
        Process exit code.
    """
    parser = build_parser()
    args = parser.parse_args(arguments)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s\t%(levelname)s\t%(name)s\t%(message)s",
    )
    if args.action == "validate":
        workflow = load_workflow_config(config_path=args.config)
        logging.getLogger(__name__).warning(
            "Configuration is valid: run=%s samples=%s output=%s",
            workflow.run_id,
            len(workflow.samples),
            workflow.output_directory,
        )
        return 0
    if args.action == "preflight":
        output = args.output
        if output is None:
            workflow = load_workflow_config(config_path=args.config)
            output = workflow.output_directory / "00_preflight" / "preflight.json"
        preflight(
            config_path=args.config,
            output_path=output,
        )
        return 0
    if args.action == "build-host-index":
        build_host_index(
            config_path=args.config,
            output_index=_required_path(parser=parser, value=args.output, option="--output"),
            output_completion=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
        )
        return 0
    if args.action == "build-minimap-reference":
        workflow = load_workflow_config(config_path=args.config)
        if workflow.minimap_genome_config is None:
            parser.error("--action build-minimap-reference requires minimap2.genome_config")
        build_minimap_reference(
            config_path=args.config,
            output_reference=args.output or workflow.minimap_reference,
            output_manifest=(
                args.reference_manifest
                or workflow.output_directory
                / "00_preflight"
                / "controlled_minimap_reference.manifest.tsv"
            ),
            output_completion=(
                args.completion
                or workflow.output_directory
                / "00_preflight"
                / "controlled_minimap_reference.complete.json"
            ),
        )
        return 0
    if args.action == "build-minimap-index":
        output = args.output
        completion = args.completion
        if output is None or completion is None:
            workflow = load_workflow_config(config_path=args.config)
            output = output or (
                workflow.output_directory / "00_preflight" / "classification_reference.mmi"
            )
            completion = completion or (
                workflow.output_directory
                / "00_preflight"
                / "classification_reference_index.complete.json"
            )
        build_minimap_index(
            config_path=args.config,
            output_index=output,
            output_completion=completion,
        )
        return 0
    if args.action == "accept-host-removed":
        completion = args.completion
        if completion is None:
            workflow = load_workflow_config(config_path=args.config)
            completion = workflow.output_directory / "01_host_depletion" / "stage.complete.json"
        accept_host_removed_batch(
            config_path=args.config,
            stage_completion=completion,
        )
        return 0
    if args.action == "host-deplete":
        host_deplete_batch(
            config_path=args.config,
            host_index=_required_path(
                parser=parser,
                value=args.host_index,
                option="--host-index",
            ),
            stage_completion=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
        )
        return 0
    if args.action == "classify":
        if args.method is None:
            parser.error("--method is required for --action classify")
        workflow = None
        sample_id = None
        if args.sample_id is not None or args.sample_index is not None:
            workflow = load_workflow_config(config_path=args.config)
            sample_id = _resolve_sample_id(
                parser=parser,
                workflow=workflow,
                sample_id=args.sample_id,
                sample_index=args.sample_index,
            )
        completion = args.completion
        if completion is None:
            if sample_id is None:
                parser.error("--completion is required when --action classify runs a batch")
            assert workflow is not None
            completion = (
                workflow.output_directory
                / "02_classification"
                / args.method
                / sample_id
                / "stage.complete.json"
            )
        classify_batch(
            config_path=args.config,
            method=args.method,
            stage_completion=completion,
            sample_id=sample_id,
        )
        return 0
    if args.action == "aggregate":
        completion = args.completion
        if completion is None:
            workflow = load_workflow_config(config_path=args.config)
            completion = workflow.output_directory / "03_final" / "workflow.complete.json"
        aggregate_results(
            config_path=args.config,
            completion_path=completion,
        )
        return 0
    if args.action == "record-scheduler-failure":
        if args.method is None:
            parser.error("--method is required for scheduler failure recording")
        if args.message is None or not args.message.strip():
            parser.error("--message is required for scheduler failure recording")
        workflow = load_workflow_config(config_path=args.config)
        sample_id = _resolve_sample_id(
            parser=parser,
            workflow=workflow,
            sample_id=args.sample_id,
            sample_index=args.sample_index,
        )
        if sample_id is None:
            parser.error("--sample-id or --sample-index is required")
        record_scheduler_failure(
            config_path=args.config,
            method=args.method,
            sample_id=sample_id,
            message=args.message,
            status=args.status,
        )
        return 0
    if args.action == "plan-slurm":
        print(json.dumps(planned_commands(config_path=args.config), indent=2))
        return 0
    if args.action == "submit-slurm":
        journal = submit_workflow(
            config_path=args.config,
            resume_submission=args.resume_submission,
            new_attempt=args.new_attempt,
        )
        logging.getLogger(__name__).warning("Slurm submission journal: %s", journal)
        return 0
    if args.action == "report":
        workflow = load_workflow_config(config_path=args.config)
        aggregate_results(
            config_path=args.config,
            completion_path=(workflow.output_directory / "03_final" / "workflow.complete.json"),
        )
        return 0
    snakefile = args.snakefile or Path(__file__).with_name("Snakefile")
    extra_arguments: tuple[str, ...] = tuple(args.target)
    if args.set_resource:
        extra_arguments += ("--set-resources", *args.set_resource)
    return run_snakemake(
        config_path=args.config,
        snakefile=snakefile,
        profile=args.profile,
        cores=args.cores,
        jobs=args.jobs,
        dry_run=args.dry_run,
        unlock=args.unlock,
        extra_arguments=extra_arguments,
    )


def _required_path(
    *,
    parser: argparse.ArgumentParser,
    value: Path | None,
    option: str,
) -> Path:
    if value is None:
        parser.error(f"{option} is required for the selected action")
    return value


def _resolve_sample_id(
    *,
    parser: argparse.ArgumentParser,
    workflow: WorkflowConfig,
    sample_id: str | None,
    sample_index: int | None,
) -> str | None:
    """Resolve mutually exclusive sample selectors against the manifest."""
    if sample_id is not None and sample_index is not None:
        parser.error("--sample-id and --sample-index are mutually exclusive")
    samples = getattr(workflow, "samples")
    if sample_index is not None:
        if sample_index < 0 or sample_index >= len(samples):
            parser.error(f"--sample-index must be between 0 and {len(samples) - 1}")
        return str(samples[sample_index].sample_id)
    if sample_id is not None:
        matches = [sample.sample_id for sample in samples if sample.sample_id == sample_id]
        if len(matches) != 1:
            parser.error(f"Unknown --sample-id: {sample_id}")
        return sample_id
    return None


if __name__ == "__main__":
    raise SystemExit(main())
