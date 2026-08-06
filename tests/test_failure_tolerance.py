"""Focused tests for classifier independence and partial-result reporting."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from helpers import build_test_project, write_fastq
from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.runtime import CommandTimeoutError
from nanopore_realdata.workflow import (
    _batch_status,
    _classifier_readiness,
    _failure_status,
    _method_status,
    _method_status_record,
    aggregate_results,
    build_minimap_index,
    classify_batch,
    preflight,
)


def _write_host_summary(*, root: Path, sample_id: str) -> None:
    """Write the core input summary required by aggregation tests."""
    directory = root / "01_host_depletion" / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "host_removal_summary.tsv").write_text(
        "sample_id\tinput_reads\tnon_host_reads\thost_reads_removed\thost_fraction\n"
        f"{sample_id}\t10\t10\t0\t0.00000000\n",
        encoding="utf-8",
    )


class TestConfigurationAndReadiness(unittest.TestCase):
    """Verify syntactic config loading and per-classifier readiness separation."""

    def test_reporting_configuration_is_optional_validated_and_deduplicated(self) -> None:
        """Old configs receive defaults while explicit focus terms remain clean."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["reporting"] = {
                "focus_taxa": ["Plasmodium", "plasmodium", "Babesia"],
                "top_n": 12,
                "max_table_rows": 250,
            }
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
            self.assertEqual(workflow.report_focus_taxa, ("Plasmodium", "Babesia"))
            self.assertEqual(workflow.report_top_n, 12)
            self.assertEqual(workflow.report_max_table_rows, 250)

            config["reporting"] = "not a mapping"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "reporting"):
                load_workflow_config(config_path=config_path)

    def test_reporting_configuration_rejects_bad_focus_and_limits(self) -> None:
        """Blank terms, non-lists and non-positive limits fail defensively."""
        cases = (
            ({"focus_taxa": "Plasmodium"}, "YAML list"),
            ({"focus_taxa": [""]}, "blank"),
            ({"focus_taxa": []}, "at least one"),
            ({"focus_taxa": ["Plasmodium"], "top_n": 0}, "positive integer"),
        )
        for reporting, message in cases:
            with self.subTest(reporting=reporting), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                config_path = build_test_project(root=root)
                config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
                config["reporting"] = reporting
                config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    load_workflow_config(config_path=config_path)

    def test_missing_classifier_resources_are_reported_not_global_config_errors(self) -> None:
        """All unavailable methods should be visible while core reads remain valid."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
                stage_resources=True,
            )
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["databases"] = {
                "kraken2": str(root / "missing_kraken"),
                "metabuli": str(root / "missing_metabuli"),
                "kmersutra_panel": str(root / "missing_panel.tsv.gz"),
            }
            config["minimap2"]["reference"] = str(root / "missing_reference.fa")
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
            output = root / "preflight" / "preflight.json"
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow._version", return_value="test"),
            ):
                preflight(config_path=config_path, output_path=output)
            payload = json.loads(output.read_text(encoding="utf-8"))
        statuses = {row["method"]: row["status"] for row in payload["classifier_readiness"]}
        self.assertEqual(statuses, {method: "unavailable" for method in statuses})
        self.assertIn("missing:", payload["resources"]["kraken2_database"])
        self.assertEqual(
            _classifier_readiness(workflow=workflow, method="minimap2")["status"],
            "unavailable",
        )


class TestRecoverableStages(unittest.TestCase):
    """Exercise index, staging and mixed-sample terminal paths."""

    def test_missing_minimap_reference_writes_unavailable_index_status(self) -> None:
        """An unbuildable minimap2 index should not block other classifiers."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["minimap2"]["reference"] = str(root / "missing.fa")
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            completion = root / "results" / "00_preflight" / "minimap.complete.json"
            index = root / "results" / "00_preflight" / "minimap.mmi"
            build_minimap_index(
                config_path=config_path,
                output_index=index,
                output_completion=completion,
            )
            payload = json.loads(completion.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "unavailable")
        self.assertFalse(index.exists())

    def test_strict_minimap_policy_raises_index_failure(self) -> None:
        """The opt-in strict policy remains available."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["minimap2"]["reference"] = str(root / "missing.fa")
            config["minimap2"]["failure_policy"] = "fail"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                build_minimap_index(
                    config_path=config_path,
                    output_index=root / "index.mmi",
                    output_completion=root / "index.complete.json",
                )

    def test_database_staging_failure_marks_every_sample(self) -> None:
        """A shared resource-copy error should become terminal sample failures."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
                stage_resources=True,
            )
            workflow = load_workflow_config(config_path=config_path)
            stage = workflow.output_directory / "02_classification" / "kraken2" / "stage.json"
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow._validate_classifier_resource"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch(
                    "nanopore_realdata.workflow.stage_resource",
                    side_effect=OSError("scratch copy failed"),
                ),
            ):
                classify_batch(
                    config_path=config_path,
                    method="kraken2",
                    stage_completion=stage,
                )
            payload = json.loads(stage.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "failed")
        self.assertEqual(payload["failed_samples"], ["sample_1"])

    def test_bundled_classifier_preserves_success_when_later_sample_fails(self) -> None:
        """One failed sample produces a partial stage rather than losing earlier work."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, input_read_state="host_removed")
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            second_fastq = root / "reads" / "second.fastq.gz"
            write_fastq(path=second_fastq, read_id="second")
            samples_path = Path(config["inputs"]["samples"])
            with samples_path.open("a", encoding="utf-8") as handle:
                handle.write(f"sample_2\t{second_fastq}\trun_1\tbarcode02\tSecond\n")
            workflow = load_workflow_config(config_path=config_path)
            stage = workflow.output_directory / "02_classification" / "kraken2" / "stage.json"
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
                    side_effect=(None, RuntimeError("second sample failed")),
                ),
            ):
                classify_batch(
                    config_path=config_path,
                    method="kraken2",
                    stage_completion=stage,
                )
            payload = json.loads(stage.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "partial")
        self.assertEqual(payload["successful_samples"], ["sample_1"])
        self.assertEqual(payload["failed_samples"], ["sample_2"])


class TestPartialAggregation(unittest.TestCase):
    """Verify final reporting remains useful with corrupt or absent methods."""

    def test_malformed_success_report_becomes_warning_and_html_still_exists(self) -> None:
        """A bad method table cannot suppress the final report suite."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, input_read_state="host_removed")
            workflow = load_workflow_config(config_path=config_path)
            _write_host_summary(root=workflow.output_directory, sample_id="sample_1")
            kraken = workflow.output_directory / "02_classification" / "kraken2" / "sample_1"
            kraken.mkdir(parents=True)
            (kraken / "complete.json").write_text(
                json.dumps({"status": "success"}), encoding="utf-8"
            )
            (kraken / "report.tsv").write_text("too\tfew\tfields\n", encoding="utf-8")
            for method in ("metabuli", "minimap2", "kmersutra"):
                directory = workflow.output_directory / "02_classification" / method / "sample_1"
                directory.mkdir(parents=True)
                (directory / "failure.json").write_text(
                    json.dumps({"status": "failed", "error": "synthetic"}),
                    encoding="utf-8",
                )
            completion = workflow.output_directory / "03_final" / "workflow.complete.json"
            aggregate_results(config_path=config_path, completion_path=completion)
            warning_text = (completion.parent / "report_warnings.tsv").read_text(encoding="utf-8")
            payload = json.loads(completion.read_text(encoding="utf-8"))
            final_report_exists = (completion.parent / "reports" / "index.html").is_file()
        self.assertIn("fewer than six fields", warning_text)
        self.assertEqual(payload["status"], "partial")
        self.assertTrue(payload["reports_generated"])
        self.assertTrue(final_report_exists)

    def test_all_classifiers_failed_still_generates_failure_report(self) -> None:
        """Zero successful classifiers should be an honest report, not no report."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, input_read_state="host_removed")
            workflow = load_workflow_config(config_path=config_path)
            _write_host_summary(root=workflow.output_directory, sample_id="sample_1")
            for method in ("kraken2", "metabuli", "minimap2", "kmersutra"):
                directory = workflow.output_directory / "02_classification" / method / "sample_1"
                directory.mkdir(parents=True)
                (directory / "failure.json").write_text(
                    json.dumps({"status": "timeout", "error": "time limit"}),
                    encoding="utf-8",
                )
            completion = workflow.output_directory / "03_final" / "workflow.complete.json"
            aggregate_results(config_path=config_path, completion_path=completion)
            payload = json.loads(completion.read_text(encoding="utf-8"))
            comparison_exists = (completion.parent / "reports" / "comparison.html").is_file()
        self.assertEqual(payload["status"], "failed")
        self.assertEqual(payload["reporting_status"], "success")
        self.assertTrue(comparison_exists)


class TestStatusHelpers(unittest.TestCase):
    """Cover terminal status edge cases used by the reports."""

    def test_status_classification_and_batch_summaries(self) -> None:
        """Timeout, unavailability, mixed and successful batches remain distinct."""
        self.assertEqual(_failure_status(error=CommandTimeoutError("late")), "timeout")
        self.assertEqual(_failure_status(error=FileNotFoundError("missing")), "unavailable")
        self.assertEqual(_failure_status(error=RuntimeError("broken")), "failed")
        self.assertEqual(_batch_status(successes=["a"], failures=[]), "success")
        self.assertEqual(
            _batch_status(successes=["a"], failures=[{"status": "failed"}]),
            "partial",
        )
        self.assertEqual(
            _batch_status(
                successes=[],
                failures=[{"status": "timeout"}, {"status": "failed"}],
            ),
            "failed",
        )

    def test_invalid_missing_and_disabled_status_records(self) -> None:
        """Bad metadata is visible and disabled KmerSutra is not treated as failed."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
                kmersutra_enabled=False,
            )
            workflow = load_workflow_config(config_path=config_path)
            sample = workflow.samples[0]
            kraken = workflow.output_directory / "02_classification" / "kraken2" / "sample_1"
            kraken.mkdir(parents=True)
            (kraken / "complete.json").write_text("not json", encoding="utf-8")
            self.assertEqual(
                _method_status(workflow=workflow, sample=sample, method="kraken2"),
                "invalid",
            )
            missing = _method_status_record(
                workflow=workflow,
                sample=sample,
                method="metabuli",
            )
            disabled = _method_status_record(
                workflow=workflow,
                sample=sample,
                method="kmersutra",
            )
        self.assertEqual(missing["status"], "missing")
        self.assertEqual(disabled["status"], "disabled")


if __name__ == "__main__":
    unittest.main()
