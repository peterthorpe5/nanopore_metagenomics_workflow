"""Branch-focused regression tests for defensive and failure behaviour."""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import yaml

from helpers import build_test_project, write_fastq
from nanopore_realdata.cli import main
from nanopore_realdata.commands import (
    kmersutra_command,
    kraken2_command,
    metabuli_command,
    pigz_command,
    samtools_count_command,
)
from nanopore_realdata.config import load_samples, load_workflow_config
from nanopore_realdata.runtime import (
    capture_output,
    filesystem_type,
    metadata_fingerprint,
    publish_directory,
    run_command,
    run_pipeline,
    stage_resource,
    validate_scratch,
)
from nanopore_realdata.workflow import (
    _capture_integer,
    _method_status,
    _method_version,
    _parse_classifier_report,
    _publish_failed_attempt,
    _publish_file,
    _require_tools,
    _validate_fastq_prefix,
    _validate_kraken_database,
    _validate_non_empty_directory,
    _validate_panel,
    _version,
    aggregate_results,
    run_snakemake,
)


class TestCliRouting(unittest.TestCase):
    """Cover every named public CLI action."""

    def test_internal_actions_route_named_paths(self) -> None:
        """Snakemake-facing actions should call exactly one adapter."""
        cases = [
            (
                ["--action", "preflight", "--config", "c", "--output", "o"],
                "nanopore_realdata.cli.preflight",
            ),
            (
                [
                    "--action",
                    "build-host-index",
                    "--config",
                    "c",
                    "--output",
                    "o",
                    "--completion",
                    "done",
                ],
                "nanopore_realdata.cli.build_host_index",
            ),
            (
                [
                    "--action",
                    "build-minimap-index",
                    "--config",
                    "c",
                    "--output",
                    "o",
                    "--completion",
                    "done",
                ],
                "nanopore_realdata.cli.build_minimap_index",
            ),
            (
                [
                    "--action",
                    "accept-host-removed",
                    "--config",
                    "c",
                    "--completion",
                    "done",
                ],
                "nanopore_realdata.cli.accept_host_removed_batch",
            ),
            (
                [
                    "--action",
                    "host-deplete",
                    "--config",
                    "c",
                    "--host-index",
                    "host.mmi",
                    "--completion",
                    "done",
                ],
                "nanopore_realdata.cli.host_deplete_batch",
            ),
            (
                [
                    "--action",
                    "classify",
                    "--config",
                    "c",
                    "--method",
                    "kraken2",
                    "--completion",
                    "done",
                ],
                "nanopore_realdata.cli.classify_batch",
            ),
            (
                ["--action", "aggregate", "--config", "c", "--completion", "done"],
                "nanopore_realdata.cli.aggregate_results",
            ),
        ]
        for arguments, target in cases:
            with self.subTest(action=arguments[1]), patch(target) as mocked:
                self.assertEqual(main(arguments), 0)
                mocked.assert_called_once()

    def test_classify_requires_method(self) -> None:
        """A classifier action without a method is invalid."""
        with self.assertRaises(SystemExit):
            main(
                [
                    "--action",
                    "classify",
                    "--config",
                    "config.yaml",
                    "--completion",
                    "done",
                ]
            )

    def test_run_action_passes_profile_targets_resources_and_unlock(self) -> None:
        """Operator overrides should reach the Snakemake launcher unchanged."""
        with patch("nanopore_realdata.cli.run_snakemake", return_value=7) as mocked:
            code = main(
                [
                    "--action",
                    "run",
                    "--config",
                    "config.yaml",
                    "--snakefile",
                    "Snakefile",
                    "--profile",
                    "profile",
                    "--cores",
                    "8",
                    "--jobs",
                    "2",
                    "--target",
                    "aggregate",
                    "--set-resource",
                    "classify_kmersutra:runtime=60",
                    "--dry-run",
                    "--unlock",
                ]
            )
        self.assertEqual(code, 7)
        kwargs = mocked.call_args.kwargs
        self.assertEqual(kwargs["cores"], "8")
        self.assertIn("--set-resources", kwargs["extra_arguments"])
        self.assertTrue(kwargs["unlock"])

    def test_report_action_uses_the_configured_final_directory(self) -> None:
        """Salvage reporting should need only a validated workflow config."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            with patch("nanopore_realdata.cli.aggregate_results") as aggregate:
                self.assertEqual(
                    main(["--action", "report", "--config", str(config_path)]),
                    0,
                )
            completion = aggregate.call_args.kwargs["completion_path"]
        self.assertEqual(
            completion,
            root / "results" / "03_final" / "workflow.complete.json",
        )

    def test_detached_actions_resolve_defaults_and_sample_indices(self) -> None:
        """Detached CLI jobs should derive durable outputs from one config."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
            )
            workflow = load_workflow_config(config_path=config_path)
            cases = (
                (
                    ["--action", "build-minimap-index", "--config", str(config_path)],
                    "nanopore_realdata.cli.build_minimap_index",
                    "output_completion",
                    workflow.output_directory
                    / "00_preflight"
                    / "classification_reference_index.complete.json",
                ),
                (
                    ["--action", "accept-host-removed", "--config", str(config_path)],
                    "nanopore_realdata.cli.accept_host_removed_batch",
                    "stage_completion",
                    workflow.output_directory / "01_host_depletion" / "stage.complete.json",
                ),
                (
                    [
                        "--action",
                        "classify",
                        "--config",
                        str(config_path),
                        "--method",
                        "kraken2",
                        "--sample-index",
                        "0",
                    ],
                    "nanopore_realdata.cli.classify_batch",
                    "stage_completion",
                    workflow.output_directory
                    / "02_classification"
                    / "kraken2"
                    / "sample_1"
                    / "stage.complete.json",
                ),
                (
                    ["--action", "aggregate", "--config", str(config_path)],
                    "nanopore_realdata.cli.aggregate_results",
                    "completion_path",
                    workflow.output_directory / "03_final" / "workflow.complete.json",
                ),
            )
            for arguments, target, output_key, expected in cases:
                with self.subTest(action=arguments[1]), patch(target) as mocked:
                    self.assertEqual(main(arguments), 0)
                    self.assertEqual(mocked.call_args.kwargs[output_key], expected)

            with patch("nanopore_realdata.cli.record_scheduler_failure") as recorder:
                self.assertEqual(
                    main(
                        [
                            "--action",
                            "record-scheduler-failure",
                            "--config",
                            str(config_path),
                            "--method",
                            "minimap2",
                            "--sample-index",
                            "0",
                            "--message",
                            "Slurm hard kill",
                        ]
                    ),
                    0,
                )
            self.assertEqual(recorder.call_args.kwargs["sample_id"], "sample_1")

            with (
                patch("nanopore_realdata.cli.planned_commands", return_value=[{"key": "x"}]),
                patch("builtins.print") as printed,
            ):
                self.assertEqual(
                    main(["--action", "plan-slurm", "--config", str(config_path)]),
                    0,
                )
            self.assertTrue(printed.called)

            journal = root / "journal.json"
            with patch("nanopore_realdata.cli.submit_workflow", return_value=journal) as submit:
                self.assertEqual(
                    main(
                        [
                            "--action",
                            "submit-slurm",
                            "--config",
                            str(config_path),
                            "--new-attempt",
                        ]
                    ),
                    0,
                )
            self.assertTrue(submit.call_args.kwargs["new_attempt"])

    def test_generated_reference_action_uses_configured_defaults(self) -> None:
        """The optional generic reference builder should remain CLI-reachable."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            genome_config = root / "genomes.tsv"
            genome_config.write_text("placeholder\n", encoding="utf-8")
            data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            data["minimap2"]["reference"] = ""
            data["minimap2"]["genome_config"] = str(genome_config)
            config_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
            with patch("nanopore_realdata.cli.build_minimap_reference") as build:
                self.assertEqual(
                    main(
                        [
                            "--action",
                            "build-minimap-reference",
                            "--config",
                            str(config_path),
                        ]
                    ),
                    0,
                )
            kwargs = build.call_args.kwargs
        self.assertEqual(kwargs["output_reference"], workflow.minimap_reference)
        self.assertEqual(
            kwargs["output_manifest"],
            workflow.output_directory
            / "00_preflight"
            / "controlled_minimap_reference.manifest.tsv",
        )


class TestConfigurationEdges(unittest.TestCase):
    """Cover validation branches that protect real input data."""

    def test_missing_and_non_mapping_configuration_fail(self) -> None:
        """Configuration loading requires an existing YAML mapping."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(FileNotFoundError):
                load_workflow_config(config_path=root / "missing.yaml")
            config = root / "config.yaml"
            config.write_text("- list\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "YAML mapping"):
                load_workflow_config(config_path=config)

    def test_missing_section_and_wrong_scalar_types_fail(self) -> None:
        """Required sections, booleans, integers, floats and choices are strict."""
        mutations = [
            (lambda data: data.pop("run"), "section 'run'"),
            (lambda data: data["execution"].__setitem__("stage_resources", "yes"), "true or false"),
            (lambda data: data["resources"].__setitem__("host_threads", 0), "positive integer"),
            (lambda data: data["kraken2"].__setitem__("confidence", "zero"), "numeric"),
            (
                lambda data: data["kmersutra"].__setitem__("call_preset", "unknown"),
                "must be one of",
            ),
        ]
        for mutation, message in mutations:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as temporary:
                config_path = build_test_project(root=Path(temporary))
                data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
                mutation(data)
                config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    load_workflow_config(config_path=config_path)

    def test_resource_path_types_and_optional_index_are_validated(self) -> None:
        """Files and directories cannot be interchanged silently."""
        path_cases = (
            ("kraken2", "file", "Kraken2 database"),
            ("metabuli", "file", "Metabuli database"),
            ("kmersutra_panel", "directory", "KmerSutra panel"),
        )
        for key, kind, message in path_cases:
            with self.subTest(key=key), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                config_path = build_test_project(root=root)
                data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
                replacement = root / f"wrong_{key}"
                if kind == "file":
                    replacement.write_text("wrong", encoding="utf-8")
                else:
                    replacement.mkdir()
                data["databases"][key] = str(replacement)
                config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    load_workflow_config(config_path=config_path)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            index = root / "host.mmi"
            index.write_text("index", encoding="utf-8")
            data["host"]["index"] = str(index)
            config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
            self.assertEqual(workflow.host_index, index.resolve())

    def test_disabled_kmersutra_allows_a_blank_panel(self) -> None:
        """Kraken2 and Metabuli should not depend on a disabled KmerSutra panel."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, kmersutra_enabled=False)
            data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            data["databases"]["kmersutra_panel"] = ""
            config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
        self.assertIsNone(workflow.kmersutra_panel)

    def test_sample_manifest_structural_failures(self) -> None:
        """Missing headers, blank paths, bad suffixes and conflicting metadata fail."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = (
                ("name\tpath\n", "requires"),
                ("sample_id\tfastq\nsample_1\t\n", "Blank fastq"),
                ("sample_id\tfastq\nsample_1\treads.txt\n", "Unsupported FASTQ"),
                ("sample_id\tfastq\n", "no data rows"),
            )
            for index, (text, message) in enumerate(cases):
                path = root / f"samples_{index}.tsv"
                path.write_text(text, encoding="utf-8")
                with (
                    self.subTest(message=message),
                    self.assertRaisesRegex(
                        (ValueError, FileNotFoundError),
                        message,
                    ),
                ):
                    load_samples(samples_path=path)

            first = root / "first.fastq"
            second = root / "second.fastq"
            write_fastq(path=first)
            write_fastq(path=second)
            conflict = root / "conflict.tsv"
            conflict.write_text(
                "sample_id\tfastq\tbarcode\n"
                f"sample_1\t{first}\tbarcode01\n"
                f"sample_1\t{second}\tbarcode02\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Conflicting barcode"):
                load_samples(samples_path=conflict)


class TestRuntimeEdges(unittest.TestCase):
    """Cover remaining runtime safety branches."""

    def test_empty_commands_and_capture_failure_are_rejected(self) -> None:
        """Subprocess helpers require valid commands and successful exits."""
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "log"
            with self.assertRaisesRegex(ValueError, "must not be empty"):
                run_command(command=[], log_path=log)
            with self.assertRaisesRegex(ValueError, "must not be empty"):
                run_pipeline(commands=[], log_path=log)
            self.assertEqual(capture_output(command=["printf", "42"]), "42")
            with self.assertRaisesRegex(RuntimeError, "failed"):
                capture_output(command=["false"])

    def test_filesystem_detection_and_successful_scratch_probe(self) -> None:
        """Filesystem discovery supports missing, failed and successful findmnt."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch("nanopore_realdata.runtime.shutil.which", return_value=None):
                self.assertIsNone(filesystem_type(path=root))
            with (
                patch("nanopore_realdata.runtime.shutil.which", return_value="findmnt"),
                patch(
                    "nanopore_realdata.runtime.subprocess.run",
                    return_value=SimpleNamespace(returncode=1, stdout=""),
                ),
            ):
                self.assertIsNone(filesystem_type(path=root))
            with (
                patch("nanopore_realdata.runtime.shutil.which", return_value="findmnt"),
                patch(
                    "nanopore_realdata.runtime.subprocess.run",
                    return_value=SimpleNamespace(returncode=0, stdout="XFS\n"),
                ),
            ):
                self.assertEqual(filesystem_type(path=root), "xfs")
            with patch("nanopore_realdata.runtime.filesystem_type", return_value="xfs"):
                result = validate_scratch(scratch_root=root, minimum_gb=1)
            self.assertEqual(result["filesystem_type"], "xfs")
            with self.assertRaisesRegex(ValueError, "not a directory"):
                validate_scratch(scratch_root=root / "missing", minimum_gb=1)

    def test_missing_fingerprint_resource_and_invalid_stage_source_fail(self) -> None:
        """Missing resources must fail before a tool starts."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(FileNotFoundError):
                metadata_fingerprint(paths=[root / "missing"], checksum_files=False)
            with self.assertRaisesRegex(ValueError, "neither a file nor directory"):
                stage_resource(
                    source=root / "missing",
                    destination_root=root / "stage",
                    log_path=root / "log",
                )

    def test_publication_rejects_non_directory_empty_and_failed_copy(self) -> None:
        """Atomic publication must never replace results with invalid staging."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_file = root / "file"
            source_file.write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not a directory"):
                publish_directory(
                    source=source_file,
                    destination=root / "out",
                    log_path=root / "log",
                )
            empty = root / "empty"
            empty.mkdir()
            with patch("nanopore_realdata.runtime.run_command"):
                with self.assertRaisesRegex(RuntimeError, "empty directory"):
                    publish_directory(
                        source=empty,
                        destination=root / "out",
                        log_path=root / "log",
                    )
            source = root / "source"
            source.mkdir()
            (source / "new").write_text("new", encoding="utf-8")
            destination = root / "destination"
            destination.mkdir()
            (destination / "old").write_text("old", encoding="utf-8")
            with patch(
                "nanopore_realdata.runtime.run_command",
                side_effect=RuntimeError("rsync failed"),
            ):
                with self.assertRaisesRegex(RuntimeError, "rsync failed"):
                    publish_directory(
                        source=source,
                        destination=destination,
                        log_path=root / "log",
                    )
            self.assertTrue((destination / "old").is_file())


class TestWorkflowEdges(unittest.TestCase):
    """Cover compact workflow validation and reporting branches."""

    def test_kmersutra_calls_are_aggregated_when_present(self) -> None:
        """Successful KmerSutra calls should be retained, not replaced by status rows."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            host = workflow.output_directory / "01_host_depletion" / "sample_1"
            host.mkdir(parents=True)
            (host / "host_removal_summary.tsv").write_text(
                "sample_id\tinput_reads\tnon_host_reads\thost_reads_removed\thost_fraction\n"
                "sample_1\t10\t8\t2\t0.2\n",
                encoding="utf-8",
            )
            for method in ("kraken2", "metabuli"):
                directory = workflow.output_directory / "02_classification" / method / "sample_1"
                directory.mkdir(parents=True)
                (directory / "report.tsv").write_text(
                    "100\t8\t8\tS\t1\tSpecies alpha\n",
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
                "sample_1\tminimap2\t1\tSpecies alpha\t8\n",
                encoding="utf-8",
            )
            (minimap / "complete.json").write_text(
                json.dumps({"status": "success"}),
                encoding="utf-8",
            )
            kmer = workflow.output_directory / "02_classification" / "kmersutra" / "sample_1"
            kmer.mkdir(parents=True)
            (kmer / "species_detection_calls.tsv").write_text(
                "species_name\tcall\nSpecies alpha\tpresent\n",
                encoding="utf-8",
            )
            (kmer / "complete.json").write_text(
                json.dumps({"status": "success"}),
                encoding="utf-8",
            )
            completion = workflow.output_directory / "03_final" / "workflow.complete.json"
            aggregate_results(config_path=config_path, completion_path=completion)
            with gzip.open(
                workflow.output_directory / "03_final" / "kmersutra_species_calls.tsv.gz",
                "rt",
                encoding="utf-8",
            ) as handle:
                text = handle.read()
        self.assertIn("Species alpha", text)
        self.assertIn("present", text)

    def test_run_snakemake_numeric_and_unlock_branches(self) -> None:
        """Local numeric cores, unlock and capacity errors should be explicit."""
        with patch(
            "nanopore_realdata.workflow.subprocess.run",
            return_value=SimpleNamespace(returncode=0),
        ) as mocked:
            self.assertEqual(
                run_snakemake(
                    config_path=Path("config"),
                    snakefile=Path("Snakefile"),
                    profile=None,
                    cores="2",
                    jobs=1,
                    dry_run=False,
                    unlock=True,
                ),
                0,
            )
        self.assertIn("--unlock", mocked.call_args.args[0])
        with self.assertRaisesRegex(ValueError, "cores must be positive"):
            run_snakemake(
                config_path=Path("config"),
                snakefile=Path("Snakefile"),
                profile=None,
                cores="0",
                jobs=1,
                dry_run=False,
                unlock=False,
            )
        with self.assertRaisesRegex(ValueError, "jobs"):
            run_snakemake(
                config_path=Path("config"),
                snakefile=Path("Snakefile"),
                profile=None,
                cores="all",
                jobs=0,
                dry_run=False,
                unlock=False,
            )

    def test_status_report_and_validation_failures(self) -> None:
        """Disabled status, malformed reports and invalid resources are reported."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root, kmersutra_enabled=False)
            workflow = load_workflow_config(config_path=config_path)
            self.assertEqual(
                _method_status(
                    workflow=workflow,
                    sample=workflow.samples[0],
                    method="kmersutra",
                ),
                "disabled",
            )
            report = root / "report.tsv"
            report.write_text("# comment\n\ninvalid\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fewer than six"):
                _parse_classifier_report(
                    path=report,
                    sample_id="sample_1",
                    method="kraken2",
                )
            bad_fastq = root / "bad.fastq"
            bad_fastq.write_text("@read\nACGT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no complete"):
                _validate_fastq_prefix(path=bad_fastq)
            bad_fastq.write_text("read\nACGT\nplus\nIIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "malformed"):
                _validate_fastq_prefix(path=bad_fastq)
            bad_panel = root / "panel.tsv"
            bad_panel.write_text("wrong\theader\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "lacks required"):
                _validate_panel(panel=bad_panel)
            empty_db = root / "empty_db"
            empty_db.mkdir()
            with self.assertRaisesRegex(ValueError, "missing or empty"):
                _validate_non_empty_directory(label="Database", path=empty_db)
            with self.assertRaisesRegex(ValueError, "missing required"):
                _validate_kraken_database(database=empty_db)

    def test_tool_presence_counts_versions_and_failure_helpers(self) -> None:
        """Executable, count, version and failed-attempt helpers cover error paths."""
        with patch("nanopore_realdata.workflow.shutil.which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "not on PATH"):
                _require_tools(action="build-host-index")
        with patch("nanopore_realdata.workflow.shutil.which", return_value="tool"):
            _require_tools(action="build-host-index")
        with patch("nanopore_realdata.workflow.capture_output", return_value="12"):
            self.assertEqual(_capture_integer(command=["count"]), 12)
        with patch("nanopore_realdata.workflow.capture_output", return_value="bad"):
            with self.assertRaisesRegex(RuntimeError, "Expected integer"):
                _capture_integer(command=["count"])
        with patch("nanopore_realdata.workflow.capture_output", return_value="v1\nmore"):
            self.assertEqual(_version(command=["tool"]), "v1")
        with patch("nanopore_realdata.workflow.capture_output", side_effect=RuntimeError("bad")):
            self.assertEqual(_version(command=["tool"]), "unavailable")
        with patch("nanopore_realdata.workflow._version", return_value="tool-version"):
            self.assertEqual(_method_version(method="kraken2"), "tool-version")
            self.assertEqual(_method_version(method="metabuli"), "tool-version")
            self.assertTrue(_method_version(method="kmersutra"))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _publish_failed_attempt(source=root / "missing", destination=root / "out")
            empty = root / "empty"
            empty.mkdir()
            _publish_failed_attempt(source=empty, destination=root / "out")

    def test_publish_file_rejects_size_mismatch(self) -> None:
        """An incomplete rsync result cannot replace a declared file."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.write_text("source", encoding="utf-8")

            def wrong_copy(*, command, log_path, stdout_path=None, timeout_seconds=None):
                del log_path, stdout_path, timeout_seconds
                Path(command[-1]).write_text("x", encoding="utf-8")

            with patch("nanopore_realdata.workflow.run_command", side_effect=wrong_copy):
                with self.assertRaisesRegex(RuntimeError, "validation failed"):
                    _publish_file(
                        source=source,
                        destination=root / "destination",
                        log_path=root / "log",
                    )


class TestCommandEdges(unittest.TestCase):
    """Cover optional command flags and remaining numeric validation."""

    def test_plain_kraken_samtools_and_parquet_kmersutra_branches(self) -> None:
        """Plain FASTQ, total counts and Parquet output remain supported."""
        kraken = kraken2_command(
            input_fastq=Path("reads.fastq"),
            database=Path("db"),
            classifications=Path("out"),
            report=Path("report"),
            threads=1,
            confidence=0.0,
        )
        self.assertNotIn("--gzip-compressed", kraken)
        count = samtools_count_command(
            bam_path=Path("host.bam"),
            threads=1,
            unmapped_only=False,
        )
        self.assertNotIn("-f", count)
        kmer = kmersutra_command(
            input_fastq=Path("reads.fastq.gz"),
            panel=Path("panel.tsv.gz"),
            sample_id="sample",
            output_directory=Path("out"),
            threads=1,
            screen_preset="exact",
            call_preset="lineage_aware",
            same_genus_fraction=0.05,
            write_parquet=True,
        )
        self.assertIn("--write_parquet_outputs", kmer)

    def test_remaining_invalid_numeric_options(self) -> None:
        """Metabuli, KmerSutra and compression reject invalid resources."""
        base = {
            "input_fastq": Path("reads.fastq"),
            "database": Path("db"),
            "output_directory": Path("out"),
            "output_prefix": "metabuli",
            "threads": 1,
            "maximum_ram_gb": 1,
            "minimum_score": 0.1,
        }
        with self.assertRaisesRegex(ValueError, "maximum RAM"):
            metabuli_command(**{**base, "maximum_ram_gb": 0})
        with self.assertRaisesRegex(ValueError, "minimum score"):
            metabuli_command(**{**base, "minimum_score": 2.0})
        with self.assertRaisesRegex(ValueError, "reportability fraction"):
            kmersutra_command(
                input_fastq=Path("reads.fastq"),
                panel=Path("panel"),
                sample_id="sample",
                output_directory=Path("out"),
                threads=1,
                screen_preset="exact",
                call_preset="lineage_aware",
                same_genus_fraction=2.0,
                write_parquet=False,
            )
        with self.assertRaisesRegex(ValueError, "Thread count"):
            pigz_command(threads=0)


if __name__ == "__main__":
    unittest.main()
