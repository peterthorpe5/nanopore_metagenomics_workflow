"""Tests for configuration validation and external command construction."""

from __future__ import annotations

import tempfile
import unittest
import csv
from pathlib import Path

import yaml

from helpers import build_test_project, write_fastq
from nanopore_realdata.commands import (
    kmersutra_command,
    kraken2_command,
    metabuli_command,
    minimap2_classification_command,
    minimap2_host_command,
    required_executables,
    samtools_count_command,
)
from nanopore_realdata.config import load_samples, load_workflow_config


class TestConfiguration(unittest.TestCase):
    """Validate YAML and repeated-row sample-sheet behaviour."""

    def test_valid_configuration_resolves_all_fields(self) -> None:
        """A complete project should load into typed settings."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            workflow = load_workflow_config(config_path=config_path)
        self.assertEqual(workflow.run_id, "test_run")
        self.assertEqual(workflow.samples[0].sample_id, "sample_1")
        self.assertEqual(workflow.threads_kmersutra, 5)
        self.assertEqual(workflow.runtime_kmersutra_minutes, 40)
        self.assertEqual(workflow.kmersutra_failure_policy, "continue")
        self.assertEqual(workflow.minimap_min_mapq, 15)
        self.assertEqual(workflow.minimap_min_alignment, 500)
        self.assertEqual(workflow.minimap_required_species, ("Species alpha",))
        self.assertIsNotNone(workflow.pcr_truth_path)

    def test_pcr_truth_is_optional_for_unlabelled_real_samples(self) -> None:
        """Classification should not require an independent benchmark truth set."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(
                root=Path(temporary),
                pcr_truth_enabled=False,
            )
            workflow = load_workflow_config(config_path=config_path)
        self.assertIsNone(workflow.pcr_truth_path)
        self.assertEqual(workflow.minimap_required_species, ("Species alpha",))

    def test_host_removed_configuration_does_not_require_a_host_reference(self) -> None:
        """Classification-ready reads must bypass host-reference requirements."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
            )
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["host"] = {"reference": "", "index": ""}
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            workflow = load_workflow_config(config_path=config_path)
        self.assertEqual(workflow.input_read_state, "host_removed")
        self.assertIsNone(workflow.host_reference)

    def test_repeated_sample_rows_preserve_part_order(self) -> None:
        """Nanopore chunks for one sample should remain ordered."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first.fastq"
            second = root / "second.fastq.gz"
            write_fastq(path=first, read_id="first")
            write_fastq(path=second, read_id="second")
            samples = root / "samples.tsv"
            samples.write_text(
                "sample_id\tfastq\tbarcode\n"
                f"sample-A\t{first}\tbarcode01\n"
                f"sample-A\t{second}\tbarcode01\n",
                encoding="utf-8",
            )
            loaded = load_samples(samples_path=samples)
        self.assertEqual(len(loaded), 1)
        self.assertEqual(loaded[0].fastq_paths, (first.resolve(), second.resolve()))

    def test_sample_sheet_rejects_duplicate_fastq(self) -> None:
        """The same physical input must not be classified twice."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq = root / "reads.fastq"
            write_fastq(path=fastq)
            samples = root / "samples.tsv"
            samples.write_text(
                f"sample_id\tfastq\nsample_1\t{fastq}\nsample_2\t{fastq}\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "listed more than once"):
                load_samples(samples_path=samples)

    def test_sample_sheet_rejects_unsafe_identifier_and_suffix(self) -> None:
        """Unsafe output paths and non-FASTQ inputs should fail early."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            file_path = root / "reads.txt"
            file_path.write_text("data", encoding="utf-8")
            samples = root / "samples.tsv"
            samples.write_text(
                f"sample_id\tfastq\nbad/sample\t{file_path}\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Invalid sample_id"):
                load_samples(samples_path=samples)

    def test_configuration_rejects_wrong_schema_and_bounds(self) -> None:
        """Schema and numeric bounds are strict rather than silently coerced."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["schema_version"] = 1
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "schema_version"):
                load_workflow_config(config_path=config_path)
            config["schema_version"] = 3
            config["kraken2"]["confidence"] = 1.2
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "between"):
                load_workflow_config(config_path=config_path)

    def test_configuration_rejects_unresolved_variable(self) -> None:
        """Unresolved environment placeholders must never become literal paths."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["run"]["output_directory"] = "${MISSING_REALDATA_TEST}/results"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unresolved environment"):
                load_workflow_config(config_path=config_path)

    def test_configuration_rejects_output_containing_inputs(self) -> None:
        """The output tree cannot enclose raw input or database paths."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            config["run"]["output_directory"] = str(root)
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must not contain"):
                load_workflow_config(config_path=config_path)

    def test_prepared_dundee_manifest_and_masked_reference_are_frozen(self) -> None:
        """Release inputs should retain all 11 samples and exact PCR exclusions."""
        repository = Path(__file__).resolve().parents[1]
        config_path = repository / "config" / "dundee_real_reads_v2_20260810.yaml"
        samples_path = config_path.with_name("dundee_real_reads_v2_20260810.samples.tsv")
        pcr_path = config_path.with_name("dundee_real_reads_v2_20260810.pcr_truth.tsv")
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        with samples_path.open("r", encoding="utf-8", newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        with pcr_path.open("r", encoding="utf-8", newline="") as handle:
            pcr = list(csv.DictReader(handle, delimiter="\t"))

        expected_reference = (
            "/home/pthorpe001/data/project_back_up_2024/Janet_genome_databases/"
            "genome_to_use/plas_outgrps_genomes_Hard_MASKED.fasta"
        )
        expected_repository = (
            "/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/"
            "2026_plasmodium_kraken_sensitivity/nanopore_metagenomics_workflow"
        )
        self.assertEqual(config["schema_version"], 3)
        self.assertEqual(
            config["deployment"]["expected_repository_root"],
            expected_repository,
        )
        self.assertEqual(config["deployment"]["expected_package_version"], "0.4.2")
        self.assertEqual(config["resources"]["kraken2_memory_mb"], 409600)
        self.assertEqual(config["minimap2"]["reference"], expected_reference)
        self.assertEqual(
            config["minimap2"]["required_species"],
            ["Plasmodium inui", "Plasmodium cynomolgi"],
        )
        self.assertEqual(config["minimap2"]["genome_config"], "")
        self.assertEqual(config["minimap2"]["index"], "")
        self.assertEqual(len(samples), 11)
        self.assertEqual(len({row["sample_id"] for row in samples}), 11)
        self.assertEqual({row["sample_id"] for row in pcr}, {row["sample_id"] for row in samples})
        included = [row for row in pcr if row["include_in_primary_comparison"] == "true"]
        excluded = [row for row in pcr if row["include_in_primary_comparison"] == "false"]
        self.assertEqual(len(included), 10)
        self.assertEqual([row["sample_id"] for row in excluded], ["MRC1123_WuH001WB"])
        self.assertTrue(all("P. inui" in row["pcr_species"] for row in included))
        self.assertEqual(sum("P. cynomolgi" in row["pcr_species"] for row in included), 6)
        self.assertEqual(config["slurm"]["kmersutra_qos"], "4week")


class TestCommands(unittest.TestCase):
    """Protect exact tool interfaces and shell-free argument vectors."""

    def test_minimap_and_samtools_commands_are_primary_only(self) -> None:
        """Host depletion should exclude secondary and supplementary records."""
        command = minimap2_host_command(
            host_index=Path("host.mmi"),
            fastq_paths=[Path("a.fastq.gz"), Path("b.fastq")],
            threads=8,
        )
        self.assertIn("--secondary=no", command)
        self.assertEqual(command[-2:], ["a.fastq.gz", "b.fastq"])
        count = samtools_count_command(
            bam_path=Path("host.bam"),
            threads=4,
            unmapped_only=True,
        )
        self.assertIn("2304", count)
        self.assertIn("4", count)

    def test_classification_minimap_command_retains_secondary_paf_hits(self) -> None:
        """Controlled-reference mapping should preserve ties for ambiguity reporting."""
        command = minimap2_classification_command(
            reference_index=Path("classification.mmi"),
            input_fastq=Path("reads.fastq.gz"),
            threads=12,
        )
        self.assertIn("--secondary=yes", command)
        self.assertIn("-c", command)
        self.assertEqual(command[-2:], ["classification.mmi", "reads.fastq.gz"])
        with self.assertRaisesRegex(ValueError, "Thread count"):
            minimap2_classification_command(
                reference_index=Path("classification.mmi"),
                input_fastq=Path("reads.fastq.gz"),
                threads=0,
            )

    def test_kraken_command_detects_gzip_and_confidence(self) -> None:
        """Compressed ONT input and explicit confidence should reach Kraken2."""
        command = kraken2_command(
            input_fastq=Path("reads.fastq.gz"),
            database=Path("db"),
            classifications=Path("classifications.tsv"),
            report=Path("report.tsv"),
            threads=12,
            confidence=0.1,
        )
        self.assertIn("--gzip-compressed", command)
        self.assertEqual(command[command.index("--confidence") + 1], "0.1")

    def test_metabuli_command_preserves_ont_mode(self) -> None:
        """Metabuli should use the established long-read sequence mode."""
        command = metabuli_command(
            input_fastq=Path("reads.fastq.gz"),
            database=Path("db"),
            output_directory=Path("out"),
            output_prefix="metabuli",
            threads=12,
            maximum_ram_gb=120,
            minimum_score=0.008,
        )
        self.assertEqual(command[:4], ["metabuli", "classify", "--seq-mode", "3"])
        self.assertEqual(command[-1], "metabuli")

    def test_kmersutra_command_is_conservative_and_compact(self) -> None:
        """The optional branch should omit bulky read-level hit tables."""
        command = kmersutra_command(
            input_fastq=Path("reads.fastq.gz"),
            panel=Path("panel.tsv.gz"),
            sample_id="sample_1",
            output_directory=Path("out"),
            threads=24,
            screen_preset="exact",
            call_preset="lineage_aware",
            same_genus_fraction=0.05,
            write_parquet=False,
        )
        self.assertIn("--no_read_level_hits", command)
        self.assertIn("--consolidate_species_calls", command)
        self.assertNotIn("--write_parquet_outputs", command)

    def test_invalid_command_values_fail(self) -> None:
        """Unsafe or nonsensical command settings must not reach subprocesses."""
        with self.assertRaisesRegex(ValueError, "At least one FASTQ"):
            minimap2_host_command(
                host_index=Path("host.mmi"),
                fastq_paths=[],
                threads=1,
            )
        with self.assertRaisesRegex(ValueError, "confidence"):
            kraken2_command(
                input_fastq=Path("reads.fastq"),
                database=Path("db"),
                classifications=Path("out"),
                report=Path("report"),
                threads=1,
                confidence=-0.1,
            )
        with self.assertRaisesRegex(ValueError, "prefix"):
            metabuli_command(
                input_fastq=Path("reads.fastq"),
                database=Path("db"),
                output_directory=Path("out"),
                output_prefix="bad/name",
                threads=1,
                maximum_ram_gb=1,
                minimum_score=0.1,
            )
        with self.assertRaisesRegex(ValueError, "Unknown workflow action"):
            required_executables(action="unknown")


if __name__ == "__main__":
    unittest.main()
