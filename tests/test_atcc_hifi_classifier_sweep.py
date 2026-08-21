"""Unit tests for the ATCC HiFi classifier operating-point sweep."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY_ROOT = Path(__file__).parents[1]
SCRIPT_PATH = REPOSITORY_ROOT / "scripts" / "atcc_hifi_classifier_sweep.py"
MANIFEST_PATH = (
    REPOSITORY_ROOT
    / "config"
    / "atcc_msa1003_hifi_classifier_operating_points_20260821.tsv"
)
SPEC = importlib.util.spec_from_file_location("atcc_hifi_classifier_sweep", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {SCRIPT_PATH}")
SWEEP = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SWEEP
SPEC.loader.exec_module(SWEEP)


class OperatingPointManifestTests(unittest.TestCase):
    """Test the locked classifier operating-point manifest."""

    def test_manifest_has_expected_complete_grid(self) -> None:
        points = SWEEP.load_operating_points(MANIFEST_PATH)
        self.assertEqual(len(points), 11)
        self.assertEqual(len(SWEEP.pending_points(points, "kraken2")), 4)
        self.assertEqual(len(SWEEP.pending_points(points, "metabuli")), 5)

    def test_metabuli_hifi_precision_preset_is_present(self) -> None:
        points = SWEEP.load_operating_points(MANIFEST_PATH)
        recommended = [
            point
            for point in points
            if point.setting_context == "hifi_recommended_current"
        ]
        self.assertEqual(len(recommended), 1)
        self.assertEqual(recommended[0].method, "metabuli")
        self.assertEqual(recommended[0].precise, 2)
        self.assertIsNone(recommended[0].min_score)
        self.assertIsNone(recommended[0].min_sp_score)

    def test_manifest_preserves_existing_baselines(self) -> None:
        points = SWEEP.load_operating_points(MANIFEST_PATH)
        reused = {(point.method, point.setting_id) for point in points if point.reuse_existing}
        self.assertEqual(
            reused,
            {
                ("kraken2", "confidence_0p00"),
                ("metabuli", "previous_min0p008_sp0p000"),
            },
        )

    def test_duplicate_setting_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "duplicate.tsv"
            header = (
                "method\tsetting_id\tconfidence\tmin_score\tmin_sp_score\t"
                "precise\treuse_existing\tsetting_context\tdescription\n"
            )
            row = "kraken2\tpoint\t0.1\t\t\t\tfalse\ttest\ttest\n"
            path.write_text(header + row + row, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate operating point"):
                SWEEP.load_operating_points(path)

    def test_cross_method_parameters_are_rejected(self) -> None:
        point = SWEEP.OperatingPoint(
            method="kraken2",
            setting_id="invalid",
            confidence=0.1,
            min_score=0.1,
            min_sp_score=None,
            precise=None,
            reuse_existing=False,
            setting_context="test",
            description="test",
        )
        with self.assertRaisesRegex(ValueError, "Metabuli settings"):
            SWEEP.validate_operating_point(point)


class CommandConstructionTests(unittest.TestCase):
    """Test shell-free classifier command construction."""

    def setUp(self) -> None:
        self.input_fastq = Path("/data/reads.fastq.gz")
        self.database = Path("/data/database")
        self.output = Path("/data/output")

    def test_kraken2_command_contains_confidence(self) -> None:
        point = SWEEP.OperatingPoint(
            method="kraken2",
            setting_id="confidence_0p10",
            confidence=0.1,
            min_score=None,
            min_sp_score=None,
            precise=None,
            reuse_existing=False,
            setting_context="test",
            description="test",
        )
        command = SWEEP.build_classifier_command(
            point,
            input_fastq=self.input_fastq,
            database=self.database,
            staging_directory=self.output,
            threads=12,
            max_ram_gb=120,
        )
        self.assertEqual(command[0], "kraken2")
        self.assertEqual(command[command.index("--confidence") + 1], "0.1")
        self.assertIn("--gzip-compressed", command)

    def test_metabuli_command_contains_both_score_thresholds(self) -> None:
        point = SWEEP.OperatingPoint(
            method="metabuli",
            setting_id="hifi_precision",
            confidence=None,
            min_score=0.07,
            min_sp_score=0.3,
            precise=0,
            reuse_existing=False,
            setting_context="hifi_recommended",
            description="test",
        )
        command = SWEEP.build_classifier_command(
            point,
            input_fastq=self.input_fastq,
            database=self.database,
            staging_directory=self.output,
            threads=12,
            max_ram_gb=120,
        )
        self.assertEqual(command[:2], ["metabuli", "classify"])
        self.assertEqual(command[command.index("--seq-mode") + 1], "3")
        self.assertEqual(command[command.index("--min-score") + 1], "0.07")
        self.assertEqual(command[command.index("--min-sp-score") + 1], "0.3")

    def test_metabuli_hifi_preset_does_not_override_scores(self) -> None:
        point = SWEEP.OperatingPoint(
            method="metabuli",
            setting_id="hifi_precise_preset2",
            confidence=None,
            min_score=None,
            min_sp_score=None,
            precise=2,
            reuse_existing=False,
            setting_context="hifi_recommended_current",
            description="test",
        )
        command = SWEEP.build_classifier_command(
            point,
            input_fastq=self.input_fastq,
            database=self.database,
            staging_directory=self.output,
            threads=12,
            max_ram_gb=120,
        )
        self.assertEqual(command[command.index("--seq-mode") + 1], "3")
        self.assertEqual(command[command.index("--precise") + 1], "2")
        self.assertNotIn("--min-score", command)
        self.assertNotIn("--min-sp-score", command)


class ReportSummaryTests(unittest.TestCase):
    """Test direct species assignment parsing and threshold summaries."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.report = self.root / "report.tsv"
        self.report.write_text(
            "9.09\t10\t10\tU\t0\tunclassified\n"
            "90.91\t100\t0\tR\t1\troot\n"
            "50.00\t55\t100\tS\t101\tExpected alpha\n"
            "1.00\t1\t1\tspecies\t202\tExpected beta\n"
            "4.55\t5\t5\tS\t303\tAdditional gamma\n"
            "1.82\t2\t2\tS\t99158\tHammondia hammondi\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_report_parser_accepts_kraken_and_metabuli_species_ranks(self) -> None:
        rows = SWEEP.parse_classifier_report(self.report)
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[-1].taxid, 99158)

    def test_summary_counts_expected_and_additional_species(self) -> None:
        metrics, species = SWEEP.summarise_report(
            SWEEP.parse_classifier_report(self.report),
            expected_taxids={101, 202},
            focus_taxid=99158,
        )
        self.assertEqual(metrics["total_reads"], 110)
        self.assertEqual(metrics["expected_species_ge_1"], 2)
        self.assertEqual(metrics["additional_species_ge_1"], 2)
        self.assertEqual(metrics["expected_species_ge_2"], 1)
        self.assertEqual(metrics["additional_species_ge_2"], 2)
        self.assertEqual(metrics["expected_species_ge_10"], 1)
        self.assertEqual(metrics["additional_species_ge_10"], 0)
        self.assertEqual(metrics["expected_species_ge_100"], 1)
        self.assertEqual(metrics["focus_taxid_direct_reads"], 2)
        self.assertEqual(len(species), 4)

    def test_malformed_report_is_rejected(self) -> None:
        self.report.write_text("1\t2\t3\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Malformed classifier report"):
            SWEEP.parse_classifier_report(self.report)

    def test_gzip_output_is_deterministic(self) -> None:
        source = self.root / "classifications.tsv"
        source.write_text("one\ttwo\n", encoding="utf-8")
        output = SWEEP.gzip_file(source)
        with gzip.open(output, "rt", encoding="utf-8") as handle:
            self.assertEqual(handle.read(), "one\ttwo\n")
        self.assertFalse(source.exists())

    def test_tsv_writer_does_not_create_comma_separated_output(self) -> None:
        output = self.root / "summary.tsv"
        SWEEP.write_tsv(output, ["first", "second"], [{"first": 1, "second": 2}])
        with output.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(rows, [["first", "second"], ["1", "2"]])


class ClassifierExecutionTests(unittest.TestCase):
    """Test atomic classifier execution with mocked external programs."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.input_fastq = self.root / "reads.fastq.gz"
        self.input_fastq.write_bytes(b"reads")
        self.database = self.root / "database"
        self.database.mkdir()
        self.output_root = self.root / "output"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def report_text() -> str:
        """Return a minimal valid classifier report."""
        return (
            "10.0\t1\t1\tU\t0\tunclassified\n"
            "90.0\t9\t0\tR\t1\troot\n"
            "80.0\t8\t8\tS\t101\tExpected species\n"
        )

    def fake_successful_run(self, command, **kwargs):
        """Create the outputs expected from a successful fake classifier."""
        if "--version" in command or command[-1:] == ["version"]:
            return subprocess.CompletedProcess(command, 0, stdout="classifier 1.0\n", stderr="")
        if command[0] == "kraken2":
            report_path = Path(command[command.index("--report") + 1])
            output_path = Path(command[command.index("--output") + 1])
            report_path.write_text(self.report_text(), encoding="utf-8")
            output_path.write_text("C\tread\t101\n", encoding="utf-8")
        elif command[:2] == ["metabuli", "classify"]:
            input_index = command.index("--seq-mode") + 2
            output_directory = Path(command[input_index + 2])
            setting_id = command[input_index + 3]
            (output_directory / f"{setting_id}_report.tsv").write_text(
                self.report_text(),
                encoding="utf-8",
            )
            (output_directory / f"{setting_id}_classifications.tsv").write_text(
                "1\tread\t101\n",
                encoding="utf-8",
            )
            (output_directory / f"{setting_id}_krona.html").write_text(
                "<html></html>\n",
                encoding="utf-8",
            )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    def test_kraken2_result_is_published_atomically(self) -> None:
        point = SWEEP.OperatingPoint(
            method="kraken2",
            setting_id="confidence_0p10",
            confidence=0.1,
            min_score=None,
            min_sp_score=None,
            precise=None,
            reuse_existing=False,
            setting_context="test",
            description="test",
        )
        with (
            mock.patch.object(SWEEP.shutil, "which", return_value="/bin/kraken2"),
            mock.patch.object(
                SWEEP.subprocess,
                "run",
                side_effect=self.fake_successful_run,
            ),
        ):
            result = SWEEP.run_operating_point(
                point,
                input_fastq=self.input_fastq,
                database=self.database,
                output_root=self.output_root,
                threads=2,
                max_ram_gb=4,
            )
        self.assertTrue((result / "complete.json").is_file())
        self.assertTrue((result / "classifications.tsv.gz").is_file())
        self.assertTrue((result / "metadata.json").is_file())
        self.assertFalse((result / "classifications.tsv").exists())

    def test_metabuli_outputs_are_normalised(self) -> None:
        point = SWEEP.OperatingPoint(
            method="metabuli",
            setting_id="hifi_precise_preset2",
            confidence=None,
            min_score=None,
            min_sp_score=None,
            precise=2,
            reuse_existing=False,
            setting_context="test",
            description="test",
        )
        with (
            mock.patch.object(SWEEP.shutil, "which", return_value="/bin/metabuli"),
            mock.patch.object(
                SWEEP.subprocess,
                "run",
                side_effect=self.fake_successful_run,
            ),
        ):
            result = SWEEP.run_operating_point(
                point,
                input_fastq=self.input_fastq,
                database=self.database,
                output_root=self.output_root,
                threads=2,
                max_ram_gb=4,
            )
        self.assertTrue((result / "report.tsv").is_file())
        self.assertTrue((result / "classifications.tsv.gz").is_file())
        self.assertTrue((result / "krona.html").is_file())

    def test_classifier_failure_retains_diagnostics(self) -> None:
        point = SWEEP.OperatingPoint(
            method="kraken2",
            setting_id="confidence_0p20",
            confidence=0.2,
            min_score=None,
            min_sp_score=None,
            precise=None,
            reuse_existing=False,
            setting_context="test",
            description="test",
        )
        failure = subprocess.CalledProcessError(returncode=1, cmd=["kraken2"])
        with (
            mock.patch.object(SWEEP.shutil, "which", return_value="/bin/kraken2"),
            mock.patch.object(SWEEP.subprocess, "run", side_effect=failure),
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                SWEEP.run_operating_point(
                    point,
                    input_fastq=self.input_fastq,
                    database=self.database,
                    output_root=self.output_root,
                    threads=2,
                    max_ram_gb=4,
                )
        failures = list((self.output_root / "kraken2").glob("*.failed.*"))
        self.assertEqual(len(failures), 1)
        self.assertTrue((failures[0] / "failure.json").is_file())


class CommandLineTests(unittest.TestCase):
    """Test named command-line actions."""

    def test_task_count_action(self) -> None:
        exit_code = SWEEP.main(
            [
                "--action",
                "task-count",
                "--manifest",
                str(MANIFEST_PATH),
                "--method",
                "kraken2",
            ]
        )
        self.assertEqual(exit_code, 0)


if __name__ == "__main__":
    unittest.main()
