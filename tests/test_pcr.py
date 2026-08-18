"""Tests for independent PCR truth validation and concordance."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from nanopore_realdata.config import Sample
from nanopore_realdata.pcr import build_pcr_concordance, load_pcr_truth


class TestPcrTruth(unittest.TestCase):
    """Validate the frozen PCR comparison contract."""

    def test_truth_requires_exact_manifest_coverage(self) -> None:
        """Every FASTQ sample must have one explicit PCR interpretation."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "truth.tsv"
            truth.write_text(
                "sample_id\tpcr_status\tpcr_species\tpcr_assay_or_source\t"
                "pcr_notes\tinclude_in_primary_comparison\n"
                "sample_1\tpositive\tP. inui\tPCR\t\ttrue\n",
                encoding="utf-8",
            )
            samples = (
                Sample("sample_1", (root / "one.fastq.gz",)),
                Sample("sample_2", (root / "two.fastq.gz",)),
            )
            with self.assertRaisesRegex(ValueError, r"missing=\['sample_2'\]"):
                load_pcr_truth(path=truth, samples=samples)

    def test_unknown_truth_is_explicitly_excluded(self) -> None:
        """An untested sample remains in the run but not the primary comparison."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "truth.tsv"
            truth.write_text(
                "sample_id\tpcr_status\tpcr_species\tpcr_assay_or_source\t"
                "pcr_notes\tinclude_in_primary_comparison\n"
                "sample_1\tunknown\t\t\tNo PCR\tfalse\n",
                encoding="utf-8",
            )
            records = load_pcr_truth(
                path=truth,
                samples=(Sample("sample_1", (root / "one.fastq.gz",)),),
            )
        self.assertEqual(records[0].pcr_status, "unknown")
        self.assertFalse(records[0].include_in_primary_comparison)


class TestPcrConcordance(unittest.TestCase):
    """Keep biological comparison separate from classifier availability."""

    def test_species_matching_is_case_insensitive_and_exact(self) -> None:
        """Capitalisation differences must not create false PCR discordance."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "truth.tsv"
            truth.write_text(
                "sample_id\tpcr_status\tpcr_species\tpcr_assay_or_source\t"
                "pcr_notes\tinclude_in_primary_comparison\n"
                "sample_1\tpositive\tP. inui; P. cynomolgi\tPCR\t\ttrue\n",
                encoding="utf-8",
            )
            records = load_pcr_truth(
                path=truth,
                samples=(Sample("sample_1", (root / "one.fastq.gz",)),),
            )
        statuses = [
            {"sample_id": "sample_1", "method": "minimap2", "status": "success"},
            {"sample_id": "sample_1", "method": "kmersutra", "status": "failed"},
        ]
        evidence = [
            {
                "sample_id": "sample_1",
                "method": "minimap2",
                "taxon_name": "Plasmodium Inui",
                "rank": "controlled_reference",
                "detected": True,
                "evidence_count": 7,
            },
            {
                "sample_id": "sample_1",
                "method": "minimap2",
                "taxon_name": "Plasmodium cynomolgi",
                "rank": "controlled_reference",
                "detected": True,
                "evidence_count": 3,
            },
        ]
        rows, summaries = build_pcr_concordance(
            truth_records=records,
            status_rows=statuses,
            evidence_rows=evidence,
            methods=("minimap2", "kmersutra"),
        )
        self.assertEqual(rows[0]["comparison_status"], "exact_species_match")
        self.assertEqual(rows[0]["expected_species_evidence_count"], "10")
        self.assertEqual(rows[1]["comparison_status"], "classifier_unavailable")
        self.assertEqual(summaries[1]["unavailable_sample_count"], "1")

    def test_negative_truth_formats_zero_evidence_without_type_failure(self) -> None:
        """PCR-negative rows must accept the integer zero returned by empty sums."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "truth.tsv"
            truth.write_text(
                "sample_id\tpcr_status\tpcr_species\tpcr_assay_or_source\t"
                "pcr_notes\tinclude_in_primary_comparison\n"
                "sample_1\tnegative\t\tPCR\tNo Plasmodium detected\ttrue\n",
                encoding="utf-8",
            )
            records = load_pcr_truth(
                path=truth,
                samples=(Sample("sample_1", (root / "one.fastq.gz",)),),
            )

        rows, summaries = build_pcr_concordance(
            truth_records=records,
            status_rows=[{"sample_id": "sample_1", "method": "minimap2", "status": "success"}],
            evidence_rows=[],
            methods=("minimap2",),
        )

        self.assertEqual(rows[0]["expected_species_evidence_count"], "0")
        self.assertEqual(rows[0]["comparison_status"], "concordant_negative")
        self.assertEqual(summaries[0]["concordant_negative_count"], "1")

    def test_bacterial_truth_uses_exact_species_ranks_only(self) -> None:
        """Evaluate bacterial truth without promoting strain descendants."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "truth.tsv"
            truth.write_text(
                "sample_id\tpcr_status\tpcr_species\tpcr_assay_or_source\t"
                "pcr_notes\tinclude_in_primary_comparison\n"
                "sample_1\tpositive\tEscherichia coli; Cereibacter sphaeroides\t"
                "reference composition\t\ttrue\n",
                encoding="utf-8",
            )
            records = load_pcr_truth(
                path=truth,
                samples=(Sample("sample_1", (root / "one.fastq.gz",)),),
            )

        evidence = [
            {
                "sample_id": "sample_1",
                "method": "metabuli",
                "taxon_name": "Escherichia coli",
                "rank": "species",
                "detected": True,
                "evidence_count": 100,
            },
            {
                "sample_id": "sample_1",
                "method": "metabuli",
                "taxon_name": "Escherichia coli K-12",
                "rank": "strain",
                "detected": True,
                "evidence_count": 10,
            },
            {
                "sample_id": "sample_1",
                "method": "metabuli",
                "taxon_name": "Shigella flexneri",
                "rank": "species",
                "detected": True,
                "evidence_count": 3,
            },
        ]
        rows, _ = build_pcr_concordance(
            truth_records=records,
            status_rows=[
                {
                    "sample_id": "sample_1",
                    "method": "metabuli",
                    "status": "success",
                }
            ],
            evidence_rows=evidence,
            methods=("metabuli",),
        )

        self.assertEqual(rows[0]["detected_expected_species"], "Escherichia coli")
        self.assertEqual(rows[0]["missed_expected_species"], "Cereibacter sphaeroides")
        self.assertEqual(rows[0]["additional_species"], "Shigella flexneri")
        self.assertEqual(rows[0]["comparison_status"], "partial_expected_detection")


if __name__ == "__main__":
    unittest.main()
