"""Validate the ATCC MSA-1003 HiFi reference-truth recovery preset."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

import yaml


class TestAtccHifiTruthPreset(unittest.TestCase):
    """Keep the reference truth aligned with the minimap2 species contract."""

    def test_truth_and_required_species_match(self) -> None:
        """Require the same 20 current taxonomy names in truth and minimap2."""
        repository_root = Path(__file__).resolve().parents[1]
        config_path = (
            repository_root
            / "config"
            / "atcc_msa1003_hifi_srr9328980_20260818.yaml"
        )
        truth_path = config_path.with_name(
            "atcc_msa1003_hifi_srr9328980_20260818.reference_truth.tsv"
        )
        with config_path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        with truth_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(config["deployment"]["expected_package_version"], "0.4.3")
        self.assertEqual(len(rows), 1)
        truth_species = {
            item.strip()
            for item in rows[0]["pcr_species"].split(";")
            if item.strip()
        }
        required_species = set(config["minimap2"]["required_species"])
        self.assertEqual(len(truth_species), 20)
        self.assertEqual(truth_species, required_species)
        for current_name in (
            "Schaalia odontolytica",
            "Phocaeicola vulgatus",
            "Cereibacter sphaeroides",
        ):
            self.assertIn(current_name, truth_species)
        self.assertIn("not PCR", rows[0]["pcr_notes"])


if __name__ == "__main__":
    unittest.main()
