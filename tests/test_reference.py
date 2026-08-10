"""Tests for controlled minimap2 reference construction."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from nanopore_realdata.minimap import parse_reference_fasta
from nanopore_realdata.reference import (
    build_controlled_reference,
    load_genome_sources,
    species_name_from_header,
    validate_required_reference_species,
    validate_required_species,
    validate_single_part_index_log,
)


class TestControlledReference(unittest.TestCase):
    """Ensure the minimap2 organism set is bounded and auditable."""

    def _write_sources(self, *, root: Path) -> Path:
        first = root / "inui.fna"
        first.write_text(">contig 1\nACGTURYN\n", encoding="utf-8")
        second = root / "cynomolgi.fna"
        second.write_text(">chromosome\nAACCGGTT\n", encoding="utf-8")
        table = root / "genomes.tsv"
        table.write_text(
            "genome_fasta\tspecies_name\ttaxid\tassembly_accession\trole\tclade\n"
            f"{first.name}\tP. inui\t5854\tASM_INUI\ttarget_species\tPlasmodium\n"
            f"{second.name}\tPlasmodium cynomolgi\t5827\tASM_CYNO\t"
            "near_neighbour\tPlasmodium\n"
            "missing.fna\tExcluded species\t1\tEXCLUDED\texclude\tExcluded\n",
            encoding="utf-8",
        )
        return table

    def test_reference_uses_non_excluded_kmersutra_sources(self) -> None:
        """Excluded rows are ignored and sequence headers retain exact taxonomy."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = load_genome_sources(config_path=self._write_sources(root=root))
            validate_required_species(
                sources=sources,
                required_species=("Plasmodium inui", "Plasmodium cynomolgi"),
            )
            reference = root / "controlled.fa"
            manifest = root / "controlled.manifest.tsv"
            stats = build_controlled_reference(
                sources=sources,
                output_fasta=reference,
                output_manifest=manifest,
                maximum_reference_bases=100,
            )
            parsed = parse_reference_fasta(path=reference)
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(stats.genome_count, 2)
        self.assertEqual(stats.reference_record_count, 2)
        self.assertEqual(stats.total_bases, 16)
        self.assertEqual(
            {record.taxon_name for record in parsed.values()},
            {"Plasmodium inui", "Plasmodium cynomolgi"},
        )
        self.assertEqual({row["source_record_id"] for row in rows}, {"contig", "chromosome"})

    def test_reference_limit_fails_without_partial_publication(self) -> None:
        """A runaway reference must not replace an existing valid destination."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = load_genome_sources(config_path=self._write_sources(root=root))
            reference = root / "controlled.fa"
            manifest = root / "controlled.manifest.tsv"
            reference.write_text("previous\n", encoding="utf-8")
            manifest.write_text("previous\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hard limit"):
                build_controlled_reference(
                    sources=sources,
                    output_fasta=reference,
                    output_manifest=manifest,
                    maximum_reference_bases=10,
                )
            self.assertEqual(reference.read_text(encoding="utf-8"), "previous\n")
            self.assertEqual(manifest.read_text(encoding="utf-8"), "previous\n")

    def test_index_log_rejects_multipart_and_oversized_builds(self) -> None:
        """The old count-inflating multipart index signature is a hard failure."""
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "minimap.log"
            log.write_text(
                "[M::mm_idx_stat] total length: 100\n[M::mm_idx_stat] total length: 100\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "exactly one part"):
                validate_single_part_index_log(
                    log_path=log,
                    maximum_reference_bases=1000,
                )
            log.write_text("[M::mm_idx_stat] total length: 1001\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "outside the allowed range"):
                validate_single_part_index_log(
                    log_path=log,
                    maximum_reference_bases=1000,
                )

    def test_masked_reference_header_variants_are_species_auditable(self) -> None:
        """Legacy focused-reference labels should support exact PCR species checks."""
        with tempfile.TemporaryDirectory() as temporary:
            reference = Path(temporary) / "masked.fa"
            reference.write_text(
                ">GCF_000524495.1_Plas_inui_San_Antonio chromosome_1\nACGT\n"
                ">contig_2 [organism=Plasmodium cynomolgi]\nACGT\n"
                ">P.knowlesi_reference\nACGT\n",
                encoding="utf-8",
            )
            species = validate_required_reference_species(
                reference_fasta=reference,
                required_species=("Plasmodium inui", "Plasmodium cynomolgi"),
            )
            with self.assertRaisesRegex(ValueError, "PCR-expected"):
                validate_required_reference_species(
                    reference_fasta=reference,
                    required_species=("Plasmodium vivax",),
                )
        self.assertEqual(
            species,
            ("Plasmodium inui", "Plasmodium cynomolgi", "Plasmodium knowlesi"),
        )
        self.assertEqual(
            species_name_from_header(header="ref taxon_name=P._inui"),
            "Plasmodium inui",
        )


if __name__ == "__main__":
    unittest.main()
