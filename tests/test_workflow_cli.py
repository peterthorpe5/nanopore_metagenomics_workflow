"""Tests for workflow adapters, failure tolerance, aggregation and CLI routing."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from helpers import build_test_project
from nanopore_realdata.cli import build_parser, main
from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.workflow import (
    _parse_classifier_report,
    _validate_fastq_prefix,
    aggregate_results,
    classify_batch,
    preflight,
    run_snakemake,
)


class TestPreflightAndAggregation(unittest.TestCase):
    """Validate public workflow boundaries with synthetic data."""

    def test_preflight_writes_resolved_inputs_and_resource_fingerprints(self) -> None:
        """Preflight should create auditable, non-recursive resource evidence."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            output = root / "preflight" / "preflight.json"
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow._version", return_value="test-version"),
            ):
                preflight(config_path=config_path, output_path=output)
            payload = json.loads(output.read_text(encoding="utf-8"))
            samples = (output.parent / "resolved_samples.tsv").read_text(encoding="utf-8")
        self.assertEqual(payload["status"], "success")
        self.assertEqual(payload["sample_count"], 1)
        self.assertIn("sample_1", samples)
        self.assertIn("kraken2_database", payload["resources"])

    def test_fastq_validation_rejects_malformed_record(self) -> None:
        """A broken first record must stop before host depletion."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "broken.fastq"
            path.write_text("@read\nACGT\n+\nIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "lengths differ"):
                _validate_fastq_prefix(path=path)

    def test_report_parser_normalises_kraken_style_rows(self) -> None:
        """Kraken2 and Metabuli report rows share one stable TSV schema."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "report.tsv"
            path.write_text("50.00\t10\t8\tS\t123\t  Species alpha\n", encoding="utf-8")
            rows = _parse_classifier_report(
                path=path,
                sample_id="sample_1",
                method="kraken2",
            )
        self.assertEqual(rows[0]["tax_id"], "123")
        self.assertEqual(rows[0]["taxon_name"], "Species alpha")

    def test_aggregate_survives_failed_kmersutra(self) -> None:
        """Kraken2 and Metabuli results remain usable when KmerSutra fails."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            host = workflow.output_directory / "01_host_depletion" / "sample_1"
            host.mkdir(parents=True)
            (host / "host_removal_summary.tsv").write_text(
                "sample_id\tinput_reads\tnon_host_reads\thost_reads_removed\thost_fraction\n"
                "sample_1\t100\t40\t60\t0.60000000\n",
                encoding="utf-8",
            )
            for method in ("kraken2", "metabuli"):
                directory = workflow.output_directory / "02_classification" / method / "sample_1"
                directory.mkdir(parents=True)
                (directory / "report.tsv").write_text(
                    "100.00\t40\t40\tS\t123\tSpecies alpha\n",
                    encoding="utf-8",
                )
                (directory / "complete.json").write_text(
                    json.dumps({"status": "success"}),
                    encoding="utf-8",
                )
            minimap = workflow.output_directory / "02_classification" / "minimap2" / "sample_1"
            minimap.mkdir(parents=True)
            (minimap / "taxon_report.tsv").write_text(
                "sample_id\tmethod\ttax_id\ttaxon_name\tbest_read_count\n"
                "sample_1\tminimap2\t123\tSpecies alpha\t4\n",
                encoding="utf-8",
            )
            (minimap / "complete.json").write_text(
                json.dumps({"status": "success"}),
                encoding="utf-8",
            )
            kmersutra = workflow.output_directory / "02_classification" / "kmersutra" / "sample_1"
            kmersutra.mkdir(parents=True)
            (kmersutra / "failure.json").write_text(
                json.dumps({"status": "failed"}),
                encoding="utf-8",
            )
            completion = workflow.output_directory / "03_final" / "workflow.complete.json"
            aggregate_results(config_path=config_path, completion_path=completion)
            sample_summary = (
                workflow.output_directory / "03_final" / "sample_summary.tsv"
            ).read_text(encoding="utf-8")
            with gzip.open(
                workflow.output_directory / "03_final" / "kmersutra_species_calls.tsv.gz",
                "rt",
                encoding="utf-8",
            ) as handle:
                calls = handle.read()
            completion_status = json.loads(completion.read_text(encoding="utf-8"))["status"]
            classifier_status = (
                workflow.output_directory / "03_final" / "classifier_status.tsv"
            ).read_text(encoding="utf-8")
        self.assertIn("failed", sample_summary)
        self.assertEqual(calls, "sample_id\tmethod\tspecies_name\tdetection_call\n")
        self.assertIn("kmersutra\tfailed", classifier_status)
        self.assertEqual(completion_status, "partial")


class TestKmerSutraFailurePolicy(unittest.TestCase):
    """Protect the explicitly non-blocking KmerSutra branch."""

    def test_disabled_kmersutra_writes_skipped_stage(self) -> None:
        """Disabling KmerSutra should require no executable or scratch staging."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, kmersutra_enabled=False)
            stage = root / "results" / "02_classification" / "kmersutra" / "stage.complete.json"
            classify_batch(
                config_path=config_path,
                method="kmersutra",
                stage_completion=stage,
            )
            payload = json.loads(stage.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "skipped")

    def test_kmersutra_failure_continues_and_records_sample(self) -> None:
        """A real KmerSutra error should produce a partial, successful stage job."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, failure_policy="continue")
            workflow = load_workflow_config(config_path=config_path)
            host = workflow.output_directory / "01_host_depletion" / "sample_1"
            host.mkdir(parents=True)
            with gzip.open(host / "non_host.fastq.gz", "wt", encoding="utf-8") as handle:
                handle.write("@read\nACGT\n+\nIIII\n")
            stage = (
                workflow.output_directory
                / "02_classification"
                / "kmersutra"
                / "stage.complete.json"
            )
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch(
                    "nanopore_realdata.workflow._classify_sample",
                    side_effect=RuntimeError("KmerSutra failed"),
                ),
            ):
                classify_batch(
                    config_path=config_path,
                    method="kmersutra",
                    stage_completion=stage,
                )
            payload = json.loads(stage.read_text(encoding="utf-8"))
            failure = stage.parent / "sample_1" / "failure.json"
            failure_exists = failure.is_file()
            self.assertEqual(payload["status"], "failed")
        self.assertEqual(payload["failed_samples"], ["sample_1"])
        self.assertTrue(failure_exists)

    def test_kmersutra_fail_policy_propagates_error(self) -> None:
        """Operators can choose fail-fast KmerSutra behaviour explicitly."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, failure_policy="fail")
            workflow = load_workflow_config(config_path=config_path)
            host = workflow.output_directory / "01_host_depletion" / "sample_1"
            host.mkdir(parents=True)
            with gzip.open(host / "non_host.fastq.gz", "wt", encoding="utf-8") as handle:
                handle.write("@read\nACGT\n+\nIIII\n")
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch(
                    "nanopore_realdata.workflow._classify_sample",
                    side_effect=RuntimeError("KmerSutra failed"),
                ),
            ):
                with self.assertRaisesRegex(RuntimeError, "KmerSutra failed"):
                    classify_batch(
                        config_path=config_path,
                        method="kmersutra",
                        stage_completion=root / "stage.json",
                    )


class TestIndependentClassifierFailurePolicies(unittest.TestCase):
    """Ensure every classifier can terminate without blocking final reporting."""

    def test_every_classifier_records_failure_under_continue_policy(self) -> None:
        """Kraken2, Metabuli, minimap2 and KmerSutra share the status contract."""
        for method in ("kraken2", "metabuli", "minimap2", "kmersutra"):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                config_path = build_test_project(root=root, input_read_state="host_removed")
                workflow = load_workflow_config(config_path=config_path)
                prepared = (
                    workflow.output_directory
                    / "01_host_depletion"
                    / "sample_1"
                    / "non_host.fastq.gz"
                )
                prepared.parent.mkdir(parents=True)
                with gzip.open(prepared, "wt", encoding="utf-8") as handle:
                    handle.write("@read\nACGT\n+\nIIII\n")
                stage = (
                    workflow.output_directory / "02_classification" / method / "stage.complete.json"
                )
                with (
                    patch("nanopore_realdata.workflow._require_tools"),
                    patch("nanopore_realdata.workflow._validate_classifier_resource"),
                    patch("nanopore_realdata.workflow.validate_scratch"),
                    patch(
                        "nanopore_realdata.workflow.stage_resource",
                        side_effect=lambda **kwargs: kwargs["source"],
                    ),
                    patch(
                        "nanopore_realdata.workflow._classify_sample",
                        side_effect=RuntimeError(f"{method} synthetic failure"),
                    ),
                ):
                    classify_batch(
                        config_path=config_path,
                        method=method,
                        stage_completion=stage,
                        sample_id="sample_1" if method == "kmersutra" else None,
                    )
                stage_payload = json.loads(stage.read_text(encoding="utf-8"))
                failure = (
                    workflow.output_directory
                    / "02_classification"
                    / method
                    / "sample_1"
                    / "failure.json"
                )
                failure_payload = json.loads(failure.read_text(encoding="utf-8"))
                self.assertEqual(stage_payload["status"], "failed")
                self.assertEqual(failure_payload["status"], "failed")
                self.assertEqual(failure_payload["method"], method)


class TestCliAndSnakemake(unittest.TestCase):
    """Verify named CLI routing and the generated Snakemake command."""

    def test_parser_requires_named_action_and_config(self) -> None:
        """The public interface should not require positional arguments."""
        args = build_parser().parse_args(["--action", "validate", "--config", "config.yaml"])
        self.assertEqual(args.action, "validate")
        self.assertEqual(args.cores, "all")

    def test_cli_validate_and_preflight_defaults(self) -> None:
        """CLI validation succeeds and preflight uses its configured default path."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            self.assertEqual(
                main(["--action", "validate", "--config", str(config_path)]),
                0,
            )
            with patch("nanopore_realdata.cli.preflight") as preflight_action:
                self.assertEqual(
                    main(["--action", "preflight", "--config", str(config_path)]),
                    0,
                )
            self.assertEqual(
                preflight_action.call_args.kwargs["output_path"],
                Path(temporary) / "results" / "00_preflight" / "preflight.json",
            )

    def test_run_snakemake_serialises_database_staging(self) -> None:
        """The launcher must always cap large database staging at one job."""
        completed = SimpleNamespace(returncode=0)
        with patch("nanopore_realdata.workflow.subprocess.run", return_value=completed) as mocked:
            code = run_snakemake(
                config_path=Path("config.yaml"),
                snakefile=Path("Snakefile"),
                profile=Path("profile"),
                cores="all",
                jobs=10,
                dry_run=True,
                unlock=False,
                extra_arguments=(
                    "--set-resources",
                    "classify_kmersutra:slurm_qos=4week",
                ),
            )
        command = mocked.call_args.args[0]
        self.assertEqual(code, 0)
        self.assertIn("database_stage=1", command)
        self.assertIn("kmersutra_worker=1", command)
        self.assertIn("--dry-run", command)
        self.assertIn("classify_kmersutra:slurm_qos=4week", command)

    def test_run_snakemake_rejects_bad_capacity_values(self) -> None:
        """Cores and job limits should fail before invoking Snakemake."""
        with self.assertRaisesRegex(ValueError, "cores"):
            run_snakemake(
                config_path=Path("config.yaml"),
                snakefile=Path("Snakefile"),
                profile=None,
                cores="zero",
                jobs=1,
                dry_run=False,
                unlock=False,
            )

    def test_snakefile_dry_run_builds_complete_dag(self) -> None:
        """The packaged Snakefile should parse and plan all real-data stages."""
        snakemake = Path(sys.executable).with_name("snakemake")
        if not snakemake.is_file():
            self.skipTest("Snakemake is not installed beside the test interpreter")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, kmersutra_enabled=False)
            snakefile = (
                Path(__file__).resolve().parents[1] / "src" / "nanopore_realdata" / "Snakefile"
            )
            snakefile_text = snakefile.read_text(encoding="utf-8")
            self.assertIn("rule classify_kmersutra:", snakefile_text)
            self.assertIn('method="kmersutra"', snakefile_text)
            kmersutra_rule = snakefile_text.split("rule classify_kmersutra:", maxsplit=1)[1]
            kmersutra_rule = kmersutra_rule.split("rule aggregate:", maxsplit=1)[0]
            self.assertNotIn("MINIMAP_INDEX_COMPLETION", kmersutra_rule)
            command = [
                str(snakemake),
                "--snakefile",
                str(snakefile),
                "--configfile",
                str(config_path),
                "--config",
                f"workflow_config_path={config_path}",
                "--dry-run",
                "--cores",
                "1",
                "--resources",
                "database_stage=1",
                "kmersutra_worker=1",
            ]
            completed = subprocess.run(
                command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env={**os.environ, "XDG_CACHE_HOME": str(root / "cache")},
            )
        self.assertEqual(completed.returncode, 0, msg=completed.stdout)
        self.assertIn("aggregate", completed.stdout)


if __name__ == "__main__":
    unittest.main()
