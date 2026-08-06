"""Unit tests for controlled-reference minimap2 summarisation."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from nanopore_realdata.minimap import (
    PafHit,
    parse_paf,
    parse_reference_fasta,
    summarise_minimap_paf,
)


def paf_row(
    *,
    query: str,
    target: str,
    matches: int,
    block_length: int,
    mapq: int,
) -> str:
    """Return one minimal valid PAF row."""
    return (
        f"{query}\t1000\t0\t{block_length}\t+\t{target}\t2000\t0\t"
        f"{block_length}\t{matches}\t{block_length}\t{mapq}\n"
    )


class TestReferenceAndPafParsing(unittest.TestCase):
    """Validate reference and alignment parsing boundaries."""

    def test_reference_metadata_handles_taxids_names_and_unknowns(self) -> None:
        """FASTA metadata should preserve exact IDs and avoid merged unknown taxa."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.fa.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(">kraken:taxid|123|ref_a\nACGT\n>reference_B.species_token\nACGT\n")
            records = parse_reference_fasta(path=path)
        self.assertEqual(records["kraken:taxid|123|ref_a"].tax_id, "123")
        self.assertEqual(records["kraken:taxid|123|ref_a"].taxon_name, "taxid_123")
        self.assertEqual(
            records["reference_B.species_token"].tax_id,
            "unknown:reference_B.species_token",
        )

    def test_reference_parser_rejects_empty_blank_and_duplicate_headers(self) -> None:
        """Invalid reference inventories should fail before minimap2 reporting."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = (
                ("empty.fa", "ACGT\n", "no records"),
                ("blank.fa", ">\nACGT\n", "Blank FASTA"),
                ("duplicate.fa", ">ref\nACGT\n>ref\nTGCA\n", "Duplicate"),
            )
            for name, text, message in cases:
                path = root / name
                path.write_text(text, encoding="utf-8")
                with self.subTest(name=name), self.assertRaisesRegex(ValueError, message):
                    parse_reference_fasta(path=path)

    def test_paf_parser_accepts_gzip_and_rejects_malformed_rows(self) -> None:
        """PAF numeric and column validation should be explicit."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid = root / "valid.paf.gz"
            with gzip.open(valid, "wt", encoding="utf-8") as handle:
                handle.write(
                    paf_row(
                        query="read_1",
                        target="ref_1",
                        matches=900,
                        block_length=950,
                        mapq=60,
                    )
                )
            hit = next(parse_paf(path=valid))
            self.assertEqual(hit, PafHit("read_1", "ref_1", 900, 950, 60))

            short = root / "short.paf"
            short.write_text("read\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fewer than 12"):
                list(parse_paf(path=short))
            numeric = root / "numeric.paf"
            numeric.write_text(
                paf_row(
                    query="read",
                    target="ref",
                    matches=1,
                    block_length=1,
                    mapq=1,
                ).replace("\t1\t1\t1\n", "\tbad\t1\t1\n"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "invalid numeric"):
                list(parse_paf(path=numeric))


class TestMinimapSummarisation(unittest.TestCase):
    """Exercise filtered streaming best-hit reporting."""

    def test_summary_counts_unique_and_ambiguous_best_taxa(self) -> None:
        """Best-hit ties across taxa should be retained as ambiguity, not certainty."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            first = "kraken:taxid|1|ref_a"
            second = "kraken:taxid|2|ref_b"
            reference.write_text(
                f">{first}\nACGT\n>{second}\nACGT\n",
                encoding="utf-8",
            )
            paf = root / "alignments.paf.gz"
            with gzip.open(paf, "wt", encoding="utf-8") as handle:
                handle.write(
                    paf_row(
                        query="read_1",
                        target=first,
                        matches=900,
                        block_length=900,
                        mapq=60,
                    )
                )
                handle.write(
                    paf_row(
                        query="read_1",
                        target=second,
                        matches=700,
                        block_length=900,
                        mapq=40,
                    )
                )
                for target in (first, second):
                    handle.write(
                        paf_row(
                            query="read_2",
                            target=target,
                            matches=800,
                            block_length=800,
                            mapq=50,
                        )
                    )
                handle.write(
                    paf_row(
                        query="read_low",
                        target=first,
                        matches=400,
                        block_length=400,
                        mapq=10,
                    )
                )
            report = root / "taxon_report.tsv"
            summary = root / "mapping_summary.tsv"
            summarise_minimap_paf(
                paf_path=paf,
                reference_fasta=reference,
                taxon_report_path=report,
                mapping_summary_path=summary,
                sample_id="sample_1",
                minimum_mapq=15,
                minimum_alignment=500,
            )
            with report.open("r", encoding="utf-8", newline="") as handle:
                rows = {row["tax_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            with summary.open("r", encoding="utf-8", newline="") as handle:
                summary_row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows["1"]["best_read_count"], "1")
        self.assertEqual(rows["1"]["ambiguous_best_read_count"], "1")
        self.assertEqual(rows["2"]["best_read_count"], "0")
        self.assertEqual(rows["2"]["ambiguous_best_read_count"], "1")
        self.assertEqual(summary_row["mapped_read_count"], "2")
        self.assertEqual(summary_row["retained_alignment_count"], "4")

    def test_summary_rejects_thresholds_and_mismatched_index_targets(self) -> None:
        """Invalid thresholds and FASTA/index mismatches should never be hidden."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            reference.write_text(">ref_a\nACGT\n", encoding="utf-8")
            paf = root / "alignments.paf"
            paf.write_text(
                paf_row(
                    query="read",
                    target="different_ref",
                    matches=900,
                    block_length=900,
                    mapq=60,
                ),
                encoding="utf-8",
            )
            arguments = {
                "paf_path": paf,
                "reference_fasta": reference,
                "taxon_report_path": root / "report.tsv",
                "mapping_summary_path": root / "summary.tsv",
                "sample_id": "sample",
                "minimum_mapq": 15,
                "minimum_alignment": 500,
            }
            with self.assertRaisesRegex(ValueError, "may not match"):
                summarise_minimap_paf(**arguments)
            with self.assertRaisesRegex(ValueError, "non-negative"):
                summarise_minimap_paf(**{**arguments, "minimum_mapq": -1})
            with self.assertRaisesRegex(ValueError, "positive"):
                summarise_minimap_paf(**{**arguments, "minimum_alignment": 0})


if __name__ == "__main__":
    unittest.main()
