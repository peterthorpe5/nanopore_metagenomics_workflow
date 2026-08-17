"""Truth labelling helpers for public ZymoBIOMICS D6300 validation.

The public Zymo mock-community screen is a species-level external validation
stress test. KmerSutra rows may represent an official expected organism, an
official expected reference label, another reference from the same expected
species, a near-neighbour, or a true non-target/off-target row. This module
keeps strict expected species/reference labels separate from compatible
same-species alternatives so external validation metrics do not credit
same-species competitor rows as official expected-target calls.
"""

from __future__ import annotations

import csv
import logging
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

from kmersutra.ai_calibration import build_call_feature_record, normalise_float
from kmersutra.table_io import read_records_table, write_records_table

EXPECTED_ZYMO_D6300_SPECIES = {
    "Pseudomonas aeruginosa",
    "Escherichia coli",
    "Salmonella enterica",
    "Limosilactobacillus fermentum",
    "Enterococcus faecalis",
    "Staphylococcus aureus",
    "Listeria monocytogenes",
    "Bacillus subtilis",
    "Saccharomyces cerevisiae",
    "Cryptococcus neoformans",
}

REFERENCE_LABEL_COLUMNS = [
    "reference_label",
    "report_label",
    "label",
    "species_name",
]

SPECIES_COLUMNS = [
    "original_species_name",
    "species",
    "species_name",
    # ``clade`` is intentionally the final fallback. Benchmark rows may use a
    # shared experimental clade (for example ``ATCC_MSA1003``) alongside the
    # actual per-row species name. Treating that clade as the species corrupts
    # truth labels and turns expected detections into apparent off-targets.
    "clade",
]

ROLE_COLUMNS = [
    "role",
    "reference_role",
    "panel_role",
]

REPORTABLE_CALL_TERMS = {
    "present",
    "species_detected",
    "detected",
    "positive",
    "reportable",
    "mixed_sample_species_detected",
}

BELOW_THRESHOLD_CALL_TERMS = {
    "observed_below_threshold",
    "below_threshold",
    "low_evidence",
    "neighbour_lineage_evidence",
    "neighbor_lineage_evidence",
}

NEAR_NEIGHBOUR_ROLES = {
    "near_neighbour",
    "near_neighbor",
    "same_genus_neighbour",
    "same_genus_neighbor",
    "same_family_neighbour",
    "same_family_neighbor",
}

OUTGROUP_ROLES = {
    "outgroup",
    "distant_outgroup",
}

TARGET_ROLE = "target_species"
SAME_SPECIES_ROLE = "same_species_competitor"


def normalise_text(value: object) -> str:
    """Return stripped text for robust category matching.

    Parameters
    ----------
    value : object
        Input value.

    Returns
    -------
    str
        Stripped text, or an empty string for missing values.
    """
    return str(value or "").strip()


def normalise_role(value: object) -> str:
    """Return a lower-case role label.

    Parameters
    ----------
    value : object
        Raw role value.

    Returns
    -------
    str
        Normalised role.
    """
    return normalise_text(value).lower()


def species_from_reference_label(value: object) -> str:
    """Infer a binomial species name from a reference-label string.

    Parameters
    ----------
    value : object
        Reference label such as ``Escherichia_coli__NRRL_B-1109``.

    Returns
    -------
    str
        Inferred species name, or an empty string if it cannot be inferred.
    """
    label = normalise_text(value)
    if not label or "__" not in label:
        return ""
    prefix = label.split("__", 1)[0]
    parts = [part for part in prefix.split("_") if part]
    if len(parts) < 2:
        return ""
    return " ".join(parts[:2])




def normalise_species_name(value: object) -> str:
    """Return a human-readable species name from row text.

    Parameters
    ----------
    value : object
        Species-like value or reference label.

    Returns
    -------
    str
        Normalised species name.
    """
    text = normalise_text(value)
    if not text:
        return ""
    if "__" in text:
        return species_from_reference_label(text)
    if " " in text:
        return text
    parts = [part for part in text.split("_") if part]
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return text


def first_non_empty(*, record: dict[str, object], columns: Iterable[str]) -> str:
    """Return the first non-empty value from a record.

    Parameters
    ----------
    record : dict[str, object]
        Input row.
    columns : iterable of str
        Candidate columns.

    Returns
    -------
    str
        First non-empty value, or an empty string.
    """
    for column in columns:
        value = normalise_text(record.get(column, ""))
        if value:
            return value
    return ""


def infer_reference_label(*, record: dict[str, object]) -> str:
    """Infer a reference label from a KmerSutra output row.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.

    Returns
    -------
    str
        Reference label if one is evident, otherwise an empty string.
    """
    for column in REFERENCE_LABEL_COLUMNS:
        value = normalise_text(record.get(column, ""))
        if "__" in value:
            return value
    return ""


def infer_species_name(*, record: dict[str, object]) -> str:
    """Infer a species name from a KmerSutra output row.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.

    Returns
    -------
    str
        Species name or inferred species name.
    """
    reference_label = infer_reference_label(record=record)
    for column in SPECIES_COLUMNS:
        value = normalise_text(record.get(column, ""))
        if value:
            species_name = normalise_species_name(value)
            if species_name:
                return species_name
    return species_from_reference_label(reference_label)


def read_reference_label_map(
    *,
    reference_label_map: str | Path | None,
    logger: logging.Logger | None = None,
) -> tuple[dict[str, str], set[str], set[str]]:
    """Read reference-label roles and expected target labels.

    Parameters
    ----------
    reference_label_map : str or pathlib.Path or None
        Optional reference-label map table.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    tuple[dict[str, str], set[str], set[str]]
        Role mapping, expected reference labels and expected species names.
    """
    roles: dict[str, str] = {}
    expected_reference_labels: set[str] = set()
    expected_species_names: set[str] = set(EXPECTED_ZYMO_D6300_SPECIES)

    if reference_label_map is None:
        if logger:
            logger.warning("No Zymo reference-label map supplied")
        return roles, expected_reference_labels, expected_species_names

    path = Path(reference_label_map)
    if not path.exists():
        if logger:
            logger.warning("Reference-label map not found: %s", path)
        return roles, expected_reference_labels, expected_species_names

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reference_label = infer_reference_label(record=dict(row))
            if not reference_label:
                continue
            role = normalise_role(first_non_empty(record=dict(row), columns=ROLE_COLUMNS))
            species_name = infer_species_name(record=dict(row))
            if role:
                roles[reference_label] = role
            if role == TARGET_ROLE:
                expected_reference_labels.add(reference_label)
                if species_name:
                    expected_species_names.add(species_name)

    if logger:
        logger.info("Read %d reference-label roles", len(roles))
        logger.info("Expected reference labels: %d", len(expected_reference_labels))
        logger.info("Expected species names: %d", len(expected_species_names))
    return roles, expected_reference_labels, expected_species_names


def infer_row_role(
    *,
    record: dict[str, object],
    reference_label_roles: dict[str, str],
) -> str:
    """Infer the panel role for a row.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.
    reference_label_roles : dict[str, str]
        Mapping from reference label to role.

    Returns
    -------
    str
        Normalised role, or an empty string.
    """
    explicit_role = normalise_role(first_non_empty(record=record, columns=ROLE_COLUMNS))
    if explicit_role:
        return explicit_role
    reference_label = infer_reference_label(record=record)
    return reference_label_roles.get(reference_label, "")


def row_has_evidence(*, record: dict[str, object]) -> bool:
    """Return whether a row has any KmerSutra evidence signal.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.

    Returns
    -------
    bool
        True if hit or positive-sequence evidence is non-zero.
    """
    evidence_columns = [
        "n_hits",
        "n_unique_kmers",
        "n_positive_sequences",
        "n_exact_hits",
        "n_fuzzy_hits",
    ]
    return any(normalise_float(record.get(column, 0.0)) > 0.0 for column in evidence_columns)


def row_is_reportable(*, record: dict[str, object]) -> bool:
    """Return whether a row appears reportable at the KmerSutra call layer.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.

    Returns
    -------
    bool
        True for reportable positive-style calls.
    """
    call = normalise_role(
        first_non_empty(
            record=record,
            columns=["call", "call_status", "species_call", "report_call"],
        )
    )
    if call in {"not_detected", "absent", "none", ""}:
        return False
    if call in BELOW_THRESHOLD_CALL_TERMS:
        return False
    if call in REPORTABLE_CALL_TERMS:
        return True
    return "reportable" in call or call.endswith("_detected") or "present" in call


def row_is_below_threshold(*, record: dict[str, object]) -> bool:
    """Return whether a row is explicitly below threshold.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row.

    Returns
    -------
    bool
        True for below-threshold/evidence-only calls.
    """
    call = normalise_role(
        first_non_empty(
            record=record,
            columns=["call", "call_status", "species_call", "report_call"],
        )
    )
    return call in BELOW_THRESHOLD_CALL_TERMS or "below" in call or "lineage" in call


def classify_zymo_truth_category(
    *,
    record: dict[str, object],
    reference_label_roles: dict[str, str],
    expected_reference_labels: set[str],
    expected_species_names: set[str] | None = None,
) -> dict[str, object]:
    """Classify one Zymo row into fine and coarse truth categories.

    Parameters
    ----------
    record : dict[str, object]
        KmerSutra row from a Zymo species-detection table.
    reference_label_roles : dict[str, str]
        Mapping from reference labels to panel roles.
    expected_reference_labels : set[str]
        Official expected Zymo reference labels.
    expected_species_names : set[str] or None, optional
        Expected organism species names.

    Returns
    -------
    dict[str, object]
        Fine category, coarse ML label and supporting flags.
    """
    expected_species = expected_species_names or EXPECTED_ZYMO_D6300_SPECIES
    reference_label = infer_reference_label(record=record)
    species_name = infer_species_name(record=record)
    role = infer_row_role(
        record=record,
        reference_label_roles=reference_label_roles,
    )
    has_evidence = row_has_evidence(record=record)
    is_reportable = row_is_reportable(record=record)
    is_below = row_is_below_threshold(record=record)
    is_expected_reference = reference_label in expected_reference_labels
    is_expected_species = species_name in expected_species

    if is_expected_reference and is_reportable:
        fine_category = "expected_reference_label"
        coarse_label = "expected_target"
    elif is_expected_species and not reference_label and is_reportable:
        fine_category = "expected_species"
        coarse_label = "expected_target"
    elif is_expected_species and reference_label and has_evidence:
        fine_category = "same_species_compatible_reference"
        coarse_label = "same_species_compatible_reference"
    elif is_reportable and not is_expected_species:
        fine_category = "true_off_target_reportable"
        coarse_label = "reportable_off_target_species"
    elif has_evidence and (role in NEAR_NEIGHBOUR_ROLES or is_below):
        fine_category = "near_neighbour_evidence"
        coarse_label = "observed_below_threshold"
    elif has_evidence and role in OUTGROUP_ROLES:
        fine_category = "non_target_background_evidence"
        coarse_label = "observed_below_threshold"
    elif has_evidence:
        fine_category = "unresolved_non_target_evidence"
        coarse_label = "observed_below_threshold"
    else:
        fine_category = "not_detected"
        coarse_label = "not_detected"

    return {
        "zymo_truth_category": fine_category,
        "ml_report_label": coarse_label,
        "zymo_species_name": species_name,
        "zymo_reference_label": reference_label,
        "zymo_reference_role": role,
        "is_expected_species": str(bool(is_expected_species)),
        "is_expected_reference_label": str(bool(is_expected_reference)),
        "has_zymo_evidence": str(bool(has_evidence)),
        "is_zymo_reportable": str(bool(is_reportable)),
    }


def build_zymo_ai_feature_records(
    *,
    records: Iterable[dict[str, object]],
    reference_label_roles: dict[str, str] | None = None,
    expected_reference_labels: set[str] | None = None,
    expected_species_names: set[str] | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Build AI feature records with clean Zymo truth categories.

    Parameters
    ----------
    records : iterable of dict
        Zymo KmerSutra call rows.
    reference_label_roles : dict[str, str] or None, optional
        Reference-label role map.
    expected_reference_labels : set[str] or None, optional
        Expected official reference labels.
    expected_species_names : set[str] or None, optional
        Expected species names.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        AI-ready feature records with fine Zymo truth labels.
    """
    roles = reference_label_roles or {}
    expected_labels = expected_reference_labels or set()
    expected_species = expected_species_names or EXPECTED_ZYMO_D6300_SPECIES
    output: list[dict[str, object]] = []

    for record in records:
        truth = classify_zymo_truth_category(
            record=record,
            reference_label_roles=roles,
            expected_reference_labels=expected_labels,
            expected_species_names=expected_species,
        )
        feature_record = build_call_feature_record(record=record)
        feature_record.update(truth)
        feature_record["public_call"] = first_non_empty(
            record=record,
            columns=["call", "call_status", "species_call", "report_call"],
        )
        feature_record["public_species_name"] = record.get("species_name", "")
        feature_record["public_report_label"] = record.get("report_label", "")
        feature_record["public_reference_label"] = record.get("reference_label", "")
        output.append(feature_record)

    if logger:
        counts = Counter(str(row["zymo_truth_category"]) for row in output)
        for category, count in sorted(counts.items()):
            logger.info("Zymo truth category %s: %d", category, count)
    return output


def count_records_by_column(
    *,
    records: Iterable[dict[str, object]],
    column: str,
) -> list[dict[str, object]]:
    """Count records by a single column.

    Parameters
    ----------
    records : iterable of dict
        Input records.
    column : str
        Column name.

    Returns
    -------
    list[dict[str, object]]
        Count table records.
    """
    counts = Counter(str(row.get(column, "")) for row in records)
    return [
        {column: value, "n_records": count}
        for value, count in sorted(counts.items())
    ]


def write_zymo_ai_feature_table(
    *,
    calls_table: str | Path,
    output_table: str | Path,
    reference_label_map: str | Path | None = None,
    category_counts_table: str | Path | None = None,
    coarse_label_counts_table: str | Path | None = None,
    logger: logging.Logger | None = None,
) -> list[dict[str, object]]:
    """Read Zymo calls and write a clean AI validation feature table.

    Parameters
    ----------
    calls_table : str or pathlib.Path
        KmerSutra species-detection calls table.
    output_table : str or pathlib.Path
        Output AI feature table.
    reference_label_map : str or pathlib.Path or None, optional
        Reference-label map with target/near-neighbour roles.
    category_counts_table : str or pathlib.Path or None, optional
        Optional fine-category count table.
    coarse_label_counts_table : str or pathlib.Path or None, optional
        Optional coarse ML-label count table.
    logger : logging.Logger or None, optional
        Logger.

    Returns
    -------
    list[dict[str, object]]
        Written feature records.
    """
    roles, expected_labels, expected_species = read_reference_label_map(
        reference_label_map=reference_label_map,
        logger=logger,
    )
    records = read_records_table(input_path=calls_table, logger=logger)
    feature_records = build_zymo_ai_feature_records(
        records=records,
        reference_label_roles=roles,
        expected_reference_labels=expected_labels,
        expected_species_names=expected_species,
        logger=logger,
    )
    if not feature_records:
        raise ValueError("No Zymo AI feature records were generated")

    write_records_table(
        records=feature_records,
        output_path=output_table,
        fieldnames=list(feature_records[0].keys()),
        logger=logger,
    )
    if category_counts_table is not None:
        write_records_table(
            records=count_records_by_column(
                records=feature_records,
                column="zymo_truth_category",
            ),
            output_path=category_counts_table,
            fieldnames=["zymo_truth_category", "n_records"],
            logger=logger,
        )
    if coarse_label_counts_table is not None:
        write_records_table(
            records=count_records_by_column(
                records=feature_records,
                column="ml_report_label",
            ),
            output_path=coarse_label_counts_table,
            fieldnames=["ml_report_label", "n_records"],
            logger=logger,
        )
    return feature_records
