"""Tests for the held-out, source-matched Kraken2 reference builder."""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).parents[1] / "scripts" / "matched_kraken2_reference.py"
SPEC = importlib.util.spec_from_file_location("matched_kraken2_reference", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {SCRIPT_PATH}")
REFERENCE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REFERENCE
SPEC.loader.exec_module(REFERENCE)


class MatchedKraken2ReferenceTests(unittest.TestCase):
    """Exercise preparation, leakage rejection and database validation."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.first_fasta = self.root / "first.fna"
        self.first_fasta.write_text(">contig_a\nACGTACGT\n", encoding="utf-8")
        self.second_fasta = self.root / "second.fna.gz"
        with gzip.open(self.second_fasta, "wt", encoding="utf-8") as handle:
            handle.write(">contig_b description\nTTGGCCAA\n")
        self.genome_config = self.root / "genomes.tsv"
        self.genome_config.write_text(
            "genome_fasta\tspecies_name\ttaxid\tassembly_accession\trole\tsource\n"
            f"{self.first_fasta}\tSpecies alpha\t101\tGCF_000001.1\ttarget\ttest\n"
            f"{self.second_fasta}\tSpecies beta\t202\tGCF_000002.1\toff_target\ttest\n",
            encoding="utf-8",
        )
        self.truth_manifest = self.root / "truth.tsv"
        self.truth_manifest.write_text(
            "species_name\taccepted_species_names\texpected\tncbi_taxid\t"
            "truth_assembly_accessions\ttruth_sequence_accessions\n"
            "Species old alpha\tSpecies alpha\ttrue\t101\t\tNC_000001.1\n"
            "Species beta\t\ttrue\t202\t\tNC_000002.1\n",
            encoding="utf-8",
        )
        self.gate = self.root / "gate.tsv"
        self.gate.write_text("check\tstatus\nreference\tPASS\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def prepare(self) -> dict[str, object]:
        """Prepare the two-genome fixture and return its summary."""
        return REFERENCE.prepare_library(
            genome_config=self.genome_config,
            truth_manifest=self.truth_manifest,
            gate_summary=self.gate,
            output_fasta=self.root / "combined.fna",
            output_manifest=self.root / "assemblies.tsv",
            output_summary=self.root / "preparation.json",
            expected_genome_count=2,
            expected_source_species_count=2,
            expected_truth_species_count=2,
        )

    def test_prepare_library_writes_taxid_headers_and_counts(self) -> None:
        summary = self.prepare()
        combined = (self.root / "combined.fna").read_text(encoding="utf-8")
        self.assertIn(">kraken:taxid|101|GCF_000001.1|g1|r1|contig_a", combined)
        self.assertIn(">kraken:taxid|202|GCF_000002.1|g2|r1|contig_b", combined)
        self.assertEqual(summary["genome_count"], 2)
        self.assertEqual(summary["source_species_count"], 2)
        self.assertEqual(summary["truth_species_represented"], 2)
        self.assertEqual(summary["held_out_sequence_leakage_count"], 0)
        self.assertEqual(summary["sequence_record_count"], 2)
        self.assertEqual(summary["sequence_base_count"], 16)

    def test_prepare_library_rejects_truth_sequence_leakage(self) -> None:
        self.first_fasta.write_text(">NC_000001.1 leaked_truth\nACGT\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Held-out truth sequence accession"):
            self.prepare()

    def test_reference_gate_must_explicitly_pass(self) -> None:
        self.gate.write_text("check\tstatus\nreference\tFAIL\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "failing status"):
            self.prepare()

    def test_duplicate_assembly_is_rejected(self) -> None:
        self.genome_config.write_text(
            "genome_fasta\tspecies_name\ttaxid\tassembly_accession\trole\tsource\n"
            f"{self.first_fasta}\tSpecies alpha\t101\tGCF_000001.1\ttarget\ttest\n"
            f"{self.second_fasta}\tSpecies beta\t202\tGCF_000001.2\toff_target\ttest\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "Assembly accession is listed more than once"):
            REFERENCE.load_genome_config(self.genome_config)

    def test_validate_database_requires_every_truth_taxid(self) -> None:
        database = self.root / "database"
        database.mkdir()
        for name in REFERENCE.KRAKEN_DATABASE_FILES:
            (database / name).write_bytes(b"database")
        inspect_report = self.root / "inspect.tsv"
        inspect_report.write_text(
            "50.00\t100\t100\tS\t101\t  Species alpha\n"
            "50.00\t100\t100\tS\t202\t  Species beta\n",
            encoding="utf-8",
        )
        preparation = self.root / "production_preparation.json"
        preparation.write_text(
            json.dumps(
                {
                    "genome_count": 475,
                    "source_species_count": 396,
                    "truth_species_count": 20,
                    "held_out_assembly_leakage_count": 0,
                    "held_out_sequence_leakage_count": 0,
                }
            ),
            encoding="utf-8",
        )
        summary = REFERENCE.validate_database(
            database=database,
            inspect_report=inspect_report,
            truth_manifest=self.truth_manifest,
            preparation_summary=preparation,
            output_summary=self.root / "complete.json",
        )
        self.assertEqual(summary["truth_taxids_represented"], 2)
        self.assertEqual(summary["status"], "complete")

        inspect_report.write_text(
            "100.00\t100\t100\tS\t101\t  Species alpha\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, r"Species beta \(202\)"):
            REFERENCE.validate_database(
                database=database,
                inspect_report=inspect_report,
                truth_manifest=self.truth_manifest,
                preparation_summary=preparation,
                output_summary=None,
            )

    def test_cli_requires_named_action_arguments(self) -> None:
        exit_code = REFERENCE.main(
            [
                "--action",
                "prepare-library",
                "--truth-manifest",
                str(self.truth_manifest),
            ]
        )
        self.assertEqual(exit_code, 2)


if __name__ == "__main__":
    unittest.main()
