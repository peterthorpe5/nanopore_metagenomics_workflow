"""Unit tests for offline, failure-aware classifier HTML reporting."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from helpers import build_test_project
from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.reporting import (
    METHODS,
    _call_is_positive,
    _comparison_matrix,
    _pairwise_overlap,
    build_normalised_evidence,
    generate_html_reports,
    serialisable_report_data,
)


class TestEvidenceNormalisation(unittest.TestCase):
    """Protect conservative cross-method evidence semantics."""

    def test_all_method_rows_are_normalised_without_forced_consensus(self) -> None:
        """Each tool retains its native metric and focus-taxon flag."""
        evidence = build_normalised_evidence(
            classifier_rows=[
                {
                    "sample_id": "sample_1",
                    "method": "kraken2",
                    "taxon_name": "Plasmodium knowlesi",
                    "tax_id": "5850",
                    "rank_code": "S",
                    "direct_reads": "7",
                    "clade_reads": "9",
                    "fraction_percent": "0.7",
                }
            ],
            minimap_rows=[
                {
                    "sample_id": "sample_1",
                    "method": "minimap2",
                    "taxon_name": "Plasmodium knowlesi",
                    "tax_id": "5850",
                    "best_read_count": "3",
                    "alignment_count": "5",
                }
            ],
            kmersutra_rows=[
                {
                    "sample_id": "sample_1",
                    "method": "kmersutra",
                    "species_name": "Plasmodium knowlesi",
                    "detection_call": "detected",
                    "supporting_read_count": "2",
                }
            ],
            focus_taxa=("Plasmodium",),
        )
        self.assertEqual([row["method"] for row in evidence], ["kraken2", "minimap2", "kmersutra"])
        self.assertTrue(all(row["is_focus"] for row in evidence))
        self.assertEqual(
            {row["metric"] for row in evidence},
            {"direct reads", "best-aligned reads", "KmerSutra evidence"},
        )
        self.assertNotIn("consensus", evidence[0])

    def test_negative_kmersutra_calls_remain_non_detected(self) -> None:
        """Explicit negative calls cannot become positive from a non-zero count."""
        evidence = build_normalised_evidence(
            classifier_rows=[],
            minimap_rows=[],
            kmersutra_rows=[
                {
                    "sample_id": "sample_1",
                    "method": "kmersutra",
                    "species_name": "Species alpha",
                    "detection_call": "not detected",
                    "supporting_read_count": "4",
                }
            ],
            focus_taxa=("Plasmodium",),
        )
        self.assertFalse(evidence[0]["detected"])
        self.assertFalse(_call_is_positive(value="no_call"))
        self.assertTrue(_call_is_positive(value="reportable"))
        self.assertFalse(_call_is_positive(value="neighbour_lineage_evidence"))
        self.assertFalse(_call_is_positive(value="background_candidate_signal"))

    def test_comparison_matrix_and_overlap_are_descriptive(self) -> None:
        """Comparison helpers expose method counts and bounded Jaccard values."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            workflow = load_workflow_config(config_path=config_path)
            evidence = [
                {
                    "sample_id": "sample_1",
                    "method": method,
                    "taxon_name": "Species alpha",
                    "evidence_count": 1,
                    "detected": True,
                    "comparable": True,
                    "is_focus": False,
                }
                for method in ("kraken2", "metabuli")
            ]
            matrix = _comparison_matrix(evidence_rows=evidence, max_rows=10)
            pairwise = _pairwise_overlap(workflow=workflow, evidence_rows=evidence)
        self.assertEqual(matrix[0]["classifier_count"], 2)
        kraken_metabuli = next(
            row
            for row in pairwise
            if row["classifier_a"] == "Kraken2" and row["classifier_b"] == "Metabuli"
        )
        self.assertEqual(kraken_metabuli["mean_jaccard"], "1.000")


class TestHtmlReports(unittest.TestCase):
    """Verify complete, offline and escaped report generation."""

    def test_report_suite_contains_final_comparison_classifier_and_sample_pages(self) -> None:
        """Every requested reporting layer should be generated and navigable."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(root=root)
            workflow = load_workflow_config(config_path=config_path)
            status_rows = [
                {
                    "sample_id": "sample_1",
                    "method": method,
                    "status": "success" if method in {"kraken2", "metabuli"} else "failed",
                    "message": "<script>alert('unsafe')</script>" if method == "minimap2" else "ok",
                    "completed_at_utc": "2026-08-06T10:00:00Z",
                }
                for method in METHODS
            ]
            evidence_rows = [
                {
                    "sample_id": "sample_1",
                    "method": "kraken2",
                    "taxon_name": "Plasmodium knowlesi",
                    "tax_id": "5850",
                    "rank": "S",
                    "evidence_count": 7,
                    "supporting_count": 9,
                    "fraction": 0.7,
                    "metric": "direct reads",
                    "detected": True,
                    "comparable": True,
                    "is_focus": True,
                }
            ]
            final_root = root / "final"
            paths = generate_html_reports(
                workflow=workflow,
                sample_rows=[{"sample_id": "sample_1"}],
                status_rows=status_rows,
                evidence_rows=evidence_rows,
                warning_rows=[{"source": "minimap2", "message": "non-fatal"}],
                final_root=final_root,
            )
            final_html = (final_root / "reports" / "index.html").read_text(encoding="utf-8")
            comparison_html = (final_root / "reports" / "comparison.html").read_text(
                encoding="utf-8"
            )
            manifest = json.loads((final_root / "report_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(paths), 8)
        self.assertIn("Classifier health", final_html)
        self.assertIn("Agreement without forced consensus", comparison_html)
        self.assertEqual(final_html.count("<section"), final_html.count("</section>"))
        self.assertEqual(comparison_html.count("<section"), comparison_html.count("</section>"))
        self.assertNotIn("<script>alert('unsafe')</script>", final_html)
        self.assertNotIn("https://", final_html)
        self.assertTrue(manifest["offline"])
        self.assertFalse(manifest["external_assets"])

    def test_serialisable_payload_accepts_non_json_scalar_values(self) -> None:
        """Durable report data converts unusual values instead of failing late."""
        payload = serialisable_report_data(
            sample_rows=[{"sample_id": Path("sample_1")}],
            status_rows=[],
            evidence_rows=[],
            warning_rows=[],
        )
        self.assertEqual(payload["sample_summary"][0]["sample_id"], "sample_1")


if __name__ == "__main__":
    unittest.main()
