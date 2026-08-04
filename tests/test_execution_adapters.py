"""Direct unit tests for host and classifier execution adapters."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from helpers import build_test_project
from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.workflow import (
    _classify_sample,
    _handle_per_read_output,
    _host_deplete_sample,
    _method_expected_outputs,
    _method_status,
    _method_task_settings,
    _publish_failed_attempt,
    _publish_file,
    _read_single_tsv,
    _run_kmersutra,
    _run_kraken2,
    _run_metabuli,
    build_host_index,
    classify_batch,
    host_deplete_batch,
)


def copy_publication(*, source: Path, destination: Path, log_path: Path) -> None:
    """Test double for atomic cross-filesystem directory publication."""
    del log_path
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


class TestStageAdapters(unittest.TestCase):
    """Exercise stage-level staging and completion behaviour."""

    def test_build_host_index_stages_builds_publishes_and_resumes(self) -> None:
        """The index action should publish once and then honour its token."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, stage_resources=True)
            output = root / "results" / "00_preflight" / "host.mmi"
            completion = output.with_suffix(".complete.json")

            def fake_run(*, command, log_path, stdout_path=None, timeout_seconds=None):
                del log_path, stdout_path, timeout_seconds
                Path(command[command.index("-d") + 1]).write_text("index", encoding="utf-8")

            def fake_publish(*, source, destination, log_path):
                del log_path
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)

            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch(
                    "nanopore_realdata.workflow.stage_resource",
                    side_effect=lambda **kwargs: kwargs["source"],
                ) as staged,
                patch("nanopore_realdata.workflow.run_command", side_effect=fake_run),
                patch("nanopore_realdata.workflow._publish_file", side_effect=fake_publish),
            ):
                build_host_index(
                    config_path=config_path,
                    output_index=output,
                    output_completion=completion,
                )
                build_host_index(
                    config_path=config_path,
                    output_index=output,
                    output_completion=completion,
                )
            self.assertEqual(output.read_text(encoding="utf-8"), "index")
            self.assertEqual(staged.call_count, 1)
            self.assertEqual(json.loads(completion.read_text(encoding="utf-8"))["status"], "success")

    def test_host_batch_stages_index_and_fastq_once(self) -> None:
        """Host input data should be copied to scratch before mapping."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, stage_resources=True)
            host_index = root / "host.mmi"
            host_index.write_text("index", encoding="utf-8")
            completion = root / "results" / "01_host_depletion" / "stage.complete.json"
            observed: list[Path] = []

            def fake_stage(*, source, destination_root, log_path):
                del destination_root, log_path
                observed.append(source)
                return source

            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch("nanopore_realdata.workflow.stage_resource", side_effect=fake_stage),
                patch("nanopore_realdata.workflow._host_deplete_sample") as worker,
            ):
                host_deplete_batch(
                    config_path=config_path,
                    host_index=host_index,
                    stage_completion=completion,
                )
            payload = json.loads(completion.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "success")
            self.assertIn(host_index, observed)
            self.assertEqual(worker.call_count, 1)
            runtime_fastqs = worker.call_args.kwargs["runtime_fastqs"]
            self.assertEqual(runtime_fastqs, (load_workflow_config(config_path=config_path).samples[0].fastq_paths[0],))

    def test_classifier_batch_stages_database_and_completes(self) -> None:
        """Kraken2 should stage one database and dispatch each sample."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, stage_resources=True)
            workflow = load_workflow_config(config_path=config_path)
            host = workflow.output_directory / "01_host_depletion" / "sample_1"
            host.mkdir(parents=True)
            with gzip.open(host / "non_host.fastq.gz", "wt", encoding="utf-8") as handle:
                handle.write("@read\nACGT\n+\nIIII\n")
            completion = workflow.output_directory / "02_classification" / "kraken2" / "stage.complete.json"
            with (
                patch("nanopore_realdata.workflow._require_tools"),
                patch("nanopore_realdata.workflow.validate_scratch"),
                patch(
                    "nanopore_realdata.workflow.stage_resource",
                    side_effect=lambda **kwargs: kwargs["source"],
                ) as staged,
                patch("nanopore_realdata.workflow._classify_sample") as worker,
            ):
                classify_batch(
                    config_path=config_path,
                    method="kraken2",
                    stage_completion=completion,
                )
            self.assertEqual(json.loads(completion.read_text(encoding="utf-8"))["status"], "success")
            self.assertEqual(staged.call_count, 1)
            self.assertEqual(worker.call_count, 1)

    def test_classifier_batch_rejects_unknown_method(self) -> None:
        """An undeclared classifier cannot enter the DAG adapter."""
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            classify_batch(
                config_path=Path("unused.yaml"),
                method="unknown",
                stage_completion=Path("unused.json"),
            )

    def test_classifier_batch_rejects_unknown_sample(self) -> None:
        """Per-sample KmerSutra jobs must name a configured sample."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            with self.assertRaisesRegex(ValueError, "Unknown sample_id"):
                classify_batch(
                    config_path=config_path,
                    method="kmersutra",
                    stage_completion=Path(temporary) / "stage.json",
                    sample_id="absent_sample",
                )


class TestHostSampleAdapter(unittest.TestCase):
    """Exercise per-sample host depletion with mocked bioinformatics tools."""

    def test_host_sample_writes_counts_metadata_and_completion(self) -> None:
        """Primary mapped and unmapped counts should produce exact host metrics."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            sample = workflow.samples[0]
            index = root / "host.mmi"
            index.write_text("index", encoding="utf-8")
            workspace = root / "workspace"
            workspace.mkdir()

            pipeline_calls = 0

            def fake_pipeline(*, commands, log_path, stdout_path=None):
                nonlocal pipeline_calls
                del commands, log_path
                pipeline_calls += 1
                if pipeline_calls == 1:
                    bam = workspace / "samples" / sample.sample_id / "host_alignment.bam"
                    bam.write_text("bam", encoding="utf-8")
                else:
                    assert stdout_path is not None
                    with gzip.open(stdout_path, "wt", encoding="utf-8") as handle:
                        handle.write("@read\nACGT\n+\nIIII\n")

            with (
                patch("nanopore_realdata.workflow.run_pipeline", side_effect=fake_pipeline),
                patch("nanopore_realdata.workflow._capture_integer", side_effect=[100, 40]),
                patch("nanopore_realdata.workflow._version", return_value="test"),
                patch("nanopore_realdata.workflow.publish_directory", side_effect=copy_publication),
            ):
                _host_deplete_sample(
                    workflow=workflow,
                    sample=sample,
                    host_index=index,
                    signature_index=index,
                    runtime_fastqs=sample.fastq_paths,
                    workspace=workspace,
                )
            final = workflow.output_directory / "01_host_depletion" / sample.sample_id
            with (final / "host_removal_summary.tsv").open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["host_reads_removed"], "60")
            self.assertEqual(row["host_fraction"], "0.60000000")
            self.assertTrue((final / "complete.json").is_file())
            self.assertFalse((final / "host_alignment.bam").exists())

    def test_host_sample_resume_and_failure_paths(self) -> None:
        """Valid prior work is skipped, while failures retain a diagnostic attempt."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            sample = workflow.samples[0]
            index = root / "host.mmi"
            index.write_text("index", encoding="utf-8")
            workspace = root / "workspace"
            workspace.mkdir()
            with (
                patch("nanopore_realdata.workflow.completion_is_valid", return_value=True),
                patch("nanopore_realdata.workflow.run_pipeline") as pipeline,
            ):
                _host_deplete_sample(
                    workflow=workflow,
                    sample=sample,
                    host_index=index,
                    signature_index=index,
                    runtime_fastqs=sample.fastq_paths,
                    workspace=workspace,
                )
            pipeline.assert_not_called()
            with (
                patch("nanopore_realdata.workflow.completion_is_valid", return_value=False),
                patch("nanopore_realdata.workflow.run_pipeline", side_effect=RuntimeError("map failed")),
                patch("nanopore_realdata.workflow._publish_failed_attempt") as failed,
            ):
                with self.assertRaisesRegex(RuntimeError, "map failed"):
                    _host_deplete_sample(
                        workflow=workflow,
                        sample=sample,
                        host_index=index,
                        signature_index=index,
                        runtime_fastqs=sample.fastq_paths,
                        workspace=workspace,
                    )
            failed.assert_called_once()


class TestClassifierSampleAdapter(unittest.TestCase):
    """Exercise per-sample dispatch, completion and failure paths."""

    def test_each_classifier_branch_publishes_expected_outputs(self) -> None:
        """Kraken2, Metabuli and KmerSutra should share restart semantics."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            sample = workflow.samples[0]
            host = workflow.output_directory / "01_host_depletion" / sample.sample_id
            host.mkdir(parents=True)
            with gzip.open(host / "non_host.fastq.gz", "wt", encoding="utf-8") as handle:
                handle.write("@read\nACGT\n+\nIIII\n")

            for method in ("kraken2", "metabuli", "kmersutra"):
                with self.subTest(method=method):
                    workspace = root / f"workspace_{method}"
                    workspace.mkdir()

                    def fake_tool(**kwargs):
                        directory = kwargs["output_directory"]
                        directory.joinpath("metadata.placeholder").write_text("x", encoding="utf-8")
                        if method in {"kraken2", "metabuli"}:
                            directory.joinpath("report.tsv").write_text(
                                "100\t1\t1\tS\t1\tSpecies\n",
                                encoding="utf-8",
                            )
                            name = (
                                "classifications.tsv.gz"
                                if method == "kraken2"
                                else "metabuli_classifications.tsv.gz"
                            )
                            directory.joinpath(name).write_text("assignments", encoding="utf-8")
                        else:
                            for name in (
                                "species_detection_calls.tsv",
                                "sample_species_kmer_evidence.tsv",
                                "sample_lineage_interpretation.tsv",
                            ):
                                directory.joinpath(name).write_text("header\n", encoding="utf-8")

                    target = f"nanopore_realdata.workflow._run_{method}"
                    with (
                        patch("nanopore_realdata.workflow.stage_resource", return_value=host / "non_host.fastq.gz"),
                        patch(target, side_effect=fake_tool),
                        patch("nanopore_realdata.workflow._method_version", return_value="test"),
                        patch("nanopore_realdata.workflow.publish_directory", side_effect=copy_publication),
                    ):
                        _classify_sample(
                            workflow=workflow,
                            sample=sample,
                            method=method,
                            resource=Path("resource"),
                            signature_resource=workflow.kmersutra_panel,
                            workspace=workspace,
                        )
                    final = workflow.output_directory / "02_classification" / method / sample.sample_id
                    self.assertTrue((final / "complete.json").is_file())

    def test_classifier_sample_rejects_missing_input_skips_and_retains_failure(self) -> None:
        """Input, resume and exception branches should remain explicit."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            sample = workflow.samples[0]
            workspace = root / "workspace"
            workspace.mkdir()
            with self.assertRaisesRegex(FileNotFoundError, "Non-host FASTQ"):
                _classify_sample(
                    workflow=workflow,
                    sample=sample,
                    method="kraken2",
                    resource=workflow.kraken_database,
                    signature_resource=workflow.kraken_database,
                    workspace=workspace,
                )
            host = workflow.output_directory / "01_host_depletion" / sample.sample_id
            host.mkdir(parents=True)
            (host / "non_host.fastq.gz").write_text("fastq", encoding="utf-8")
            with (
                patch("nanopore_realdata.workflow.completion_is_valid", return_value=True),
                patch("nanopore_realdata.workflow.stage_resource") as stage,
            ):
                _classify_sample(
                    workflow=workflow,
                    sample=sample,
                    method="kraken2",
                    resource=workflow.kraken_database,
                    signature_resource=workflow.kraken_database,
                    workspace=workspace,
                )
            stage.assert_not_called()
            with (
                patch("nanopore_realdata.workflow.completion_is_valid", return_value=False),
                patch("nanopore_realdata.workflow.stage_resource", return_value=host / "non_host.fastq.gz"),
                patch("nanopore_realdata.workflow._run_kraken2", side_effect=RuntimeError("failed")),
                patch("nanopore_realdata.workflow._publish_failed_attempt") as failed,
            ):
                with self.assertRaisesRegex(RuntimeError, "failed"):
                    _classify_sample(
                        workflow=workflow,
                        sample=sample,
                        method="kraken2",
                        resource=workflow.kraken_database,
                        signature_resource=workflow.kraken_database,
                        workspace=workspace,
                    )
            failed.assert_called_once()


class TestToolAdaptersAndHelpers(unittest.TestCase):
    """Cover tool-specific output contracts and compact helper branches."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.config_path = build_test_project(root=self.root)
        self.workflow = load_workflow_config(config_path=self.config_path)
        self.sample = self.workflow.samples[0]

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_kraken_metabuli_and_kmersutra_output_contracts(self) -> None:
        """Each adapter should accept expected outputs and reject absent reports."""
        for method in ("kraken2", "metabuli", "kmersutra"):
            with self.subTest(method=method):
                output = self.root / f"tool_{method}"
                output.mkdir()

                def fake_run(*, command, log_path, stdout_path=None, timeout_seconds=None):
                    del command, log_path, stdout_path
                    if method == "kraken2":
                        (output / "classifications.tsv").write_text("raw", encoding="utf-8")
                        (output / "report.tsv").write_text("report", encoding="utf-8")
                    elif method == "metabuli":
                        (output / "metabuli_classifications.tsv").write_text("raw", encoding="utf-8")
                        (output / "metabuli_report.tsv").write_text("report", encoding="utf-8")
                    else:
                        self.assertEqual(
                            timeout_seconds,
                            self.workflow.kmersutra_timeout_minutes * 60,
                        )
                        (output / "species_detection_calls.tsv").write_text("calls", encoding="utf-8")

                with (
                    patch("nanopore_realdata.workflow.run_command", side_effect=fake_run),
                    patch("nanopore_realdata.workflow._handle_per_read_output"),
                ):
                    if method == "kraken2":
                        _run_kraken2(
                            workflow=self.workflow,
                            sample=self.sample,
                            input_fastq=Path("input.fastq.gz"),
                            database=Path("db"),
                            output_directory=output,
                            log_path=output / "log",
                        )
                    elif method == "metabuli":
                        _run_metabuli(
                            workflow=self.workflow,
                            sample=self.sample,
                            input_fastq=Path("input.fastq.gz"),
                            database=Path("db"),
                            output_directory=output,
                            log_path=output / "log",
                        )
                        self.assertTrue((output / "report.tsv").is_file())
                    else:
                        _run_kmersutra(
                            workflow=self.workflow,
                            sample=self.sample,
                            input_fastq=Path("input.fastq.gz"),
                            panel=Path("panel.tsv.gz"),
                            output_directory=output,
                            log_path=output / "log",
                        )

    def test_tool_adapters_reject_missing_declared_output(self) -> None:
        """Successful exit status alone is not sufficient evidence of success."""
        for method, function in (
            ("kraken2", _run_kraken2),
            ("metabuli", _run_metabuli),
            ("kmersutra", _run_kmersutra),
        ):
            with self.subTest(method=method):
                output = self.root / f"missing_{method}"
                output.mkdir()
                kwargs = {
                    "workflow": self.workflow,
                    "sample": self.sample,
                    "input_fastq": Path("input.fastq.gz"),
                    "output_directory": output,
                    "log_path": output / "log",
                }
                kwargs["panel" if method == "kmersutra" else "database"] = Path("resource")
                with (
                    patch("nanopore_realdata.workflow.run_command"),
                    patch("nanopore_realdata.workflow._handle_per_read_output"),
                ):
                    with self.assertRaisesRegex(RuntimeError, "produced no"):
                        function(**kwargs)

    def test_per_read_output_is_compressed_or_deleted(self) -> None:
        """Large assignment tables should never remain uncompressed."""
        raw = self.root / "classifications.tsv"
        raw.write_text("assignments", encoding="utf-8")

        def fake_compress(*, command, log_path, stdout_path=None, timeout_seconds=None):
            del command, log_path, timeout_seconds
            assert stdout_path is not None
            stdout_path.write_text("compressed", encoding="utf-8")

        with patch("nanopore_realdata.workflow.run_command", side_effect=fake_compress):
            _handle_per_read_output(
                raw_path=raw,
                keep=True,
                threads=2,
                log_path=self.root / "log",
            )
        self.assertFalse(raw.exists())
        self.assertTrue((self.root / "classifications.tsv.gz").is_file())
        raw.write_text("assignments", encoding="utf-8")
        _handle_per_read_output(
            raw_path=raw,
            keep=False,
            threads=2,
            log_path=self.root / "log",
        )
        self.assertFalse(raw.exists())
        with self.assertRaisesRegex(RuntimeError, "did not create"):
            _handle_per_read_output(
                raw_path=raw,
                keep=True,
                threads=2,
                log_path=self.root / "log",
            )

    def test_method_helpers_cover_all_methods_and_statuses(self) -> None:
        """Expected outputs, settings and statuses should remain method-specific."""
        for method in ("kraken2", "metabuli", "kmersutra"):
            outputs = _method_expected_outputs(
                directory=Path("out"),
                method=method,
                keep_per_read=True,
            )
            settings = _method_task_settings(
                workflow=self.workflow,
                sample=self.sample,
                method=method,
            )
            self.assertTrue(outputs)
            self.assertEqual(settings["sample_id"], "sample_1")
        self.assertEqual(
            _method_expected_outputs(
                directory=Path("out"),
                method="kraken2",
                keep_per_read=False,
            ),
            [Path("out/metadata.json"), Path("out/report.tsv")],
        )
        self.assertEqual(
            _method_status(workflow=self.workflow, sample=self.sample, method="kraken2"),
            "missing",
        )
        directory = self.workflow.output_directory / "02_classification" / "kraken2" / "sample_1"
        directory.mkdir(parents=True)
        (directory / "complete.json").write_text("not-json", encoding="utf-8")
        self.assertEqual(
            _method_status(workflow=self.workflow, sample=self.sample, method="kraken2"),
            "invalid",
        )

    def test_publish_file_failed_attempt_and_single_tsv_validation(self) -> None:
        """Small publication and diagnostic helpers should remain defensive."""
        source = self.root / "source.txt"
        source.write_text("new", encoding="utf-8")
        destination = self.root / "published" / "result.txt"

        def fake_rsync(*, command, log_path, stdout_path=None, timeout_seconds=None):
            del log_path, stdout_path, timeout_seconds
            shutil.copy2(command[-2], command[-1])

        with patch("nanopore_realdata.workflow.run_command", side_effect=fake_rsync):
            _publish_file(
                source=source,
                destination=destination,
                log_path=self.root / "publish.log",
            )
        self.assertEqual(destination.read_text(encoding="utf-8"), "new")

        failed_source = self.root / "failed_source"
        failed_source.mkdir()
        (failed_source / "tool.log").write_text("error", encoding="utf-8")
        (failed_source / "large.partial").write_text("do not retain", encoding="utf-8")
        failed_destination = self.root / "method" / "sample"
        _publish_failed_attempt(source=failed_source, destination=failed_destination)
        attempts = list((failed_destination.parent / "failed_attempts").glob("sample.*"))
        self.assertTrue(attempts)
        self.assertTrue((attempts[0] / "tool.log").is_file())
        self.assertFalse((attempts[0] / "large.partial").exists())

        empty = self.root / "empty.tsv"
        empty.write_text("column\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            _read_single_tsv(path=empty)


if __name__ == "__main__":
    unittest.main()
