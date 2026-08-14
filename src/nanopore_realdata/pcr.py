"""PCR truth-table validation and classifier concordance summaries."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from nanopore_realdata.config import SAMPLE_ID_PATTERN, Sample


PCR_STATUSES = frozenset({"positive", "negative", "unresolved", "unknown"})
BOOLEAN_TEXT = {"true": True, "false": False}
METHODS = ("kraken2", "metabuli", "minimap2", "kmersutra")
PLASMODIUM_ABBREVIATION = re.compile(r"^P\.\s*", flags=re.IGNORECASE)


@dataclass(frozen=True)
class PcrTruth:
    """Independent PCR interpretation for one manifest sample."""

    sample_id: str
    pcr_status: str
    pcr_species: tuple[str, ...]
    pcr_species_source_text: str
    pcr_assay_or_source: str
    pcr_notes: str
    include_in_primary_comparison: bool


def load_pcr_truth(
    *,
    path: Path,
    samples: Sequence[Sample],
) -> tuple[PcrTruth, ...]:
    """Load a complete, sample-keyed PCR truth table.

    Args:
        path: Tab-separated PCR truth table.
        samples: Validated workflow samples in manifest order.

    Returns:
        PCR records ordered to match ``samples``.

    Raises:
        FileNotFoundError: If the truth table is unavailable.
        ValueError: If the table is incomplete, duplicated or inconsistent.
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"PCR truth table is missing or empty: {resolved}")
    required = {
        "sample_id",
        "pcr_status",
        "pcr_species",
        "pcr_assay_or_source",
        "pcr_notes",
        "include_in_primary_comparison",
    }
    records: dict[str, PcrTruth] = {}
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValueError(f"PCR truth table is missing required columns: {missing}")
        for row_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not SAMPLE_ID_PATTERN.fullmatch(sample_id):
                raise ValueError(f"Invalid PCR sample_id on row {row_number}: {sample_id!r}")
            if sample_id in records:
                raise ValueError(f"Duplicate PCR truth row for sample: {sample_id}")
            status = (row.get("pcr_status") or "").strip().casefold()
            if status not in PCR_STATUSES:
                raise ValueError(
                    f"Invalid pcr_status on row {row_number}: {status!r}; "
                    f"expected one of {sorted(PCR_STATUSES)}"
                )
            source_species = (row.get("pcr_species") or "").strip()
            species = _split_species(value=source_species)
            if status == "positive" and not species:
                raise ValueError(f"PCR-positive sample has no pcr_species: {sample_id}")
            if status == "negative" and species:
                raise ValueError(f"PCR-negative sample must not list pcr_species: {sample_id}")
            include_text = (row.get("include_in_primary_comparison") or "").strip().casefold()
            if include_text not in BOOLEAN_TEXT:
                raise ValueError(
                    f"include_in_primary_comparison must be true or false on row {row_number}"
                )
            include = BOOLEAN_TEXT[include_text]
            if include and status not in {"positive", "negative"}:
                raise ValueError(
                    f"Non-terminal PCR status cannot enter the primary comparison: {sample_id}"
                )
            records[sample_id] = PcrTruth(
                sample_id=sample_id,
                pcr_status=status,
                pcr_species=species,
                pcr_species_source_text=source_species,
                pcr_assay_or_source=(row.get("pcr_assay_or_source") or "").strip(),
                pcr_notes=(row.get("pcr_notes") or "").strip(),
                include_in_primary_comparison=include,
            )
    manifest_ids = [sample.sample_id for sample in samples]
    missing = [sample_id for sample_id in manifest_ids if sample_id not in records]
    extras = sorted(set(records).difference(manifest_ids))
    if missing or extras:
        raise ValueError(
            "PCR truth sample coverage differs from the FASTQ manifest; "
            f"missing={missing}, extra={extras}"
        )
    return tuple(records[sample_id] for sample_id in manifest_ids)


def pcr_truth_rows(*, records: Sequence[PcrTruth]) -> list[dict[str, str]]:
    """Serialise validated PCR records to a deterministic TSV-ready form.

    Args:
        records: Validated PCR truth records.

    Returns:
        Ordered dictionaries containing source and canonical species labels.
    """
    return [
        {
            "sample_id": record.sample_id,
            "pcr_status": record.pcr_status,
            "pcr_species": record.pcr_species_source_text,
            "pcr_species_canonical": "; ".join(record.pcr_species),
            "pcr_assay_or_source": record.pcr_assay_or_source,
            "pcr_notes": record.pcr_notes,
            "include_in_primary_comparison": str(record.include_in_primary_comparison).lower(),
        }
        for record in records
    ]


def build_pcr_concordance(
    *,
    truth_records: Sequence[PcrTruth],
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
    methods: Sequence[str] = METHODS,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Compare independent PCR truth with each classifier separately.

    Classifier failure and missing output remain unavailable observations and
    are never converted into biological non-detections.

    Args:
        truth_records: Validated PCR records.
        status_rows: Per-sample classifier terminal states.
        evidence_rows: Normalised classifier taxon evidence.
        methods: Ordered classifier names.

    Returns:
        Per-sample concordance rows and per-method exact-count summaries.
    """
    statuses = {
        (str(row.get("sample_id", "")), str(row.get("method", ""))): str(
            row.get("status", "missing")
        )
        for row in status_rows
    }
    detected: dict[tuple[str, str], dict[str, str]] = {}
    evidence_counts: dict[tuple[str, str, str], float] = {}
    for row in evidence_rows:
        if not _as_bool(row.get("detected")):
            continue
        sample_id = str(row.get("sample_id", ""))
        method = str(row.get("method", ""))
        species = canonical_species_name(value=str(row.get("taxon_name", "")))
        if not species.casefold().startswith("plasmodium "):
            continue
        species_key = species.casefold()
        detected.setdefault((sample_id, method), {}).setdefault(
            species_key,
            species,
        )
        key = (sample_id, method, species_key)
        evidence_counts[key] = evidence_counts.get(key, 0.0) + _number(row.get("evidence_count"))

    concordance: list[dict[str, str]] = []
    for truth in truth_records:
        expected = {species.casefold(): species for species in truth.pcr_species}
        for method in methods:
            method_status = statuses.get((truth.sample_id, method), "missing")
            observed = detected.get((truth.sample_id, method), {})
            found_keys = expected.keys() & observed.keys()
            missed_keys = expected.keys() - observed.keys()
            additional_keys = observed.keys() - expected.keys()
            found = {expected[key] for key in found_keys}
            missed = {expected[key] for key in missed_keys}
            additional = {observed[key] for key in additional_keys}
            comparison_status = _comparison_status(
                truth=truth,
                classifier_status=method_status,
                found=found,
                missed=missed,
                additional=additional,
            )
            expected_evidence = sum(
                evidence_counts.get((truth.sample_id, method, species_key), 0.0)
                for species_key in expected
            )
            concordance.append(
                {
                    "sample_id": truth.sample_id,
                    "method": method,
                    "pcr_status": truth.pcr_status,
                    "pcr_species": truth.pcr_species_source_text,
                    "pcr_species_canonical": "; ".join(sorted(expected.values())),
                    "include_in_primary_comparison": str(
                        truth.include_in_primary_comparison
                    ).lower(),
                    "classifier_status": method_status,
                    "detected_plasmodium_species": "; ".join(sorted(observed.values())),
                    "detected_expected_species": "; ".join(sorted(found)),
                    "missed_expected_species": "; ".join(sorted(missed)),
                    "additional_plasmodium_species": "; ".join(sorted(additional)),
                    "expected_species_count": str(len(expected)),
                    "detected_expected_species_count": str(len(found)),
                    "expected_species_evidence_count": _display_number(expected_evidence),
                    "all_expected_species_detected": _boolean_result(
                        method_status == "success" and not missed
                    ),
                    "species_exact_match": _boolean_result(
                        method_status == "success" and observed.keys() == expected.keys()
                    ),
                    "comparison_status": comparison_status,
                }
            )
    return concordance, _summarise_concordance(rows=concordance, methods=methods)


def canonical_species_name(*, value: str) -> str:
    """Return a conservative canonical species label for exact comparison.

    Args:
        value: Source taxon label.

    Returns:
        Whitespace-normalised label with ``P.`` expanded to ``Plasmodium``.
    """
    cleaned = " ".join(value.replace("_", " ").split())
    cleaned = PLASMODIUM_ABBREVIATION.sub("Plasmodium ", cleaned)
    return cleaned.strip()


def _split_species(*, value: str) -> tuple[str, ...]:
    species: list[str] = []
    seen: set[str] = set()
    for item in value.split(";"):
        canonical = canonical_species_name(value=item)
        if not canonical:
            continue
        folded = canonical.casefold()
        if folded not in seen:
            seen.add(folded)
            species.append(canonical)
    return tuple(species)


def _comparison_status(
    *,
    truth: PcrTruth,
    classifier_status: str,
    found: set[str],
    missed: set[str],
    additional: set[str],
) -> str:
    if not truth.include_in_primary_comparison:
        return "excluded"
    if truth.pcr_status not in {"positive", "negative"}:
        return "not_evaluable"
    if classifier_status != "success":
        return "classifier_unavailable"
    if truth.pcr_status == "negative":
        return "concordant_negative" if not additional else "unexpected_detection"
    if not missed and not additional:
        return "exact_species_match"
    if not missed:
        return "all_expected_plus_additional"
    if found:
        return "partial_expected_detection"
    return "expected_species_not_detected"


def _summarise_concordance(
    *,
    rows: Sequence[Mapping[str, str]],
    methods: Sequence[str],
) -> list[dict[str, str]]:
    summaries: list[dict[str, str]] = []
    for method in methods:
        method_rows = [
            row
            for row in rows
            if row["method"] == method and row["include_in_primary_comparison"] == "true"
        ]
        available = [row for row in method_rows if row["classifier_status"] == "success"]
        positive = [row for row in available if row["pcr_status"] == "positive"]
        negative = [row for row in available if row["pcr_status"] == "negative"]
        complete = sum(row["all_expected_species_detected"] == "true" for row in positive)
        exact = sum(row["species_exact_match"] == "true" for row in available)
        concordant_negative = sum(
            row["comparison_status"] == "concordant_negative" for row in negative
        )
        summaries.append(
            {
                "method": method,
                "primary_sample_count": str(len(method_rows)),
                "available_sample_count": str(len(available)),
                "unavailable_sample_count": str(len(method_rows) - len(available)),
                "pcr_positive_available_count": str(len(positive)),
                "all_expected_species_detected_count": str(complete),
                "exact_species_match_count": str(exact),
                "pcr_negative_available_count": str(len(negative)),
                "concordant_negative_count": str(concordant_negative),
            }
        )
    return summaries


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().casefold() in {"true", "yes", "1", "positive"}


def _number(value: Any) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return 0.0
    if not math.isfinite(numeric):
        return 0.0
    return max(0.0, numeric)


def _display_number(value: float | int) -> str:
    """Format finite integer or floating-point evidence without false decimals."""
    numeric = float(value)
    if not math.isfinite(numeric):
        return "0"
    return str(int(numeric)) if numeric.is_integer() else f"{numeric:.6g}"


def _boolean_result(value: bool) -> str:
    return str(value).lower()


def expected_pcr_species(*, records: Iterable[PcrTruth]) -> tuple[str, ...]:
    """Return unique canonical species required by included positive samples.

    Args:
        records: PCR truth records.

    Returns:
        Case-insensitively de-duplicated species names in first-seen order.
    """
    species: list[str] = []
    seen: set[str] = set()
    for record in records:
        if not record.include_in_primary_comparison or record.pcr_status != "positive":
            continue
        for item in record.pcr_species:
            folded = item.casefold()
            if folded not in seen:
                seen.add(folded)
                species.append(item)
    return tuple(species)
