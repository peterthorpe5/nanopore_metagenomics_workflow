"""Command-line interface for the real Nanopore Snakemake workflow."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.workflow import (
    accept_host_removed_batch,
    aggregate_results,
    build_host_index,
    build_minimap_index,
    classify_batch,
    host_deplete_batch,
    preflight,
    run_snakemake,
)


ACTIONS = (
    "validate",
    "preflight",
    "build-host-index",
    "build-minimap-index",
    "accept-host-removed",
    "host-deplete",
    "classify",
    "aggregate",
    "report",
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
        preflight(
            config_path=args.config,
            output_path=_required_path(parser=parser, value=args.output, option="--output"),
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
    if args.action == "build-minimap-index":
        build_minimap_index(
            config_path=args.config,
            output_index=_required_path(parser=parser, value=args.output, option="--output"),
            output_completion=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
        )
        return 0
    if args.action == "accept-host-removed":
        accept_host_removed_batch(
            config_path=args.config,
            stage_completion=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
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
        classify_batch(
            config_path=args.config,
            method=args.method,
            stage_completion=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
            sample_id=args.sample_id,
        )
        return 0
    if args.action == "aggregate":
        aggregate_results(
            config_path=args.config,
            completion_path=_required_path(
                parser=parser,
                value=args.completion,
                option="--completion",
            ),
        )
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


if __name__ == "__main__":
    raise SystemExit(main())
