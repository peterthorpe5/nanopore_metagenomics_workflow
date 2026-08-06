"""Controlled-reference minimap2 parsing and taxon-level summarisation."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


TAXID_PATTERN = re.compile(r"kraken:taxid\|([^|,\s]+)")
SPECIES_TOKEN_PATTERN = re.compile(r"(?:^|_)([A-Z]\.[A-Za-z0-9._-]+)$")


@dataclass(frozen=True)
class ReferenceRecord:
    """Taxonomic metadata parsed from one reference FASTA header."""

    reference_name: str
    tax_id: str
    taxon_name: str


@dataclass(frozen=True)
class PafHit:
    """Fields required to select and summarise one minimap2 PAF hit."""

    query_name: str
    target_name: str
    matches: int
    block_length: int
    mapq: int

    @property
    def score(self) -> tuple[int, int, int]:
        """Return the deterministic best-hit ordering used by the benchmark."""
        return (self.mapq, self.matches, self.block_length)


def parse_reference_fasta(*, path: Path) -> dict[str, ReferenceRecord]:
    """Parse reference identifiers and available taxon metadata.

    Args:
        path: Classification reference FASTA.

    Returns:
        Mapping from exact sequence identifier to parsed metadata.

    Raises:
        ValueError: If the FASTA has no records or duplicate identifiers.
    """
    records: dict[str, ReferenceRecord] = {}
    with _open_text(path=path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Blank FASTA header at line {line_number}: {path}")
            reference_name = header.split()[0]
            if reference_name in records:
                raise ValueError(f"Duplicate FASTA reference identifier {reference_name!r}: {path}")
            tax_id = _tax_id(header=header)
            records[reference_name] = ReferenceRecord(
                reference_name=reference_name,
                tax_id=tax_id,
                taxon_name=_taxon_name(
                    reference_name=reference_name,
                    header=header,
                    tax_id=tax_id,
                ),
            )
    if not records:
        raise ValueError(f"Classification reference FASTA contains no records: {path}")
    return records


def parse_paf(*, path: Path) -> Iterator[PafHit]:
    """Yield validated minimap2 PAF hits from a plain or gzip file.

    Args:
        path: Minimap2 PAF or PAF.GZ path.

    Yields:
        Parsed minimal PAF records.

    Raises:
        ValueError: If a non-empty row is malformed.
    """
    with _open_text(path=path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 fields: {path}")
            try:
                yield PafHit(
                    query_name=fields[0],
                    target_name=fields[5],
                    matches=int(fields[9]),
                    block_length=int(fields[10]),
                    mapq=int(fields[11]),
                )
            except ValueError as error:
                raise ValueError(
                    f"PAF line {line_number} has invalid numeric fields: {path}"
                ) from error


def summarise_minimap_paf(
    *,
    paf_path: Path,
    reference_fasta: Path,
    taxon_report_path: Path,
    mapping_summary_path: Path,
    sample_id: str,
    minimum_mapq: int,
    minimum_alignment: int,
) -> None:
    """Summarise filtered best alignments by reported reference taxon.

    Minimap2 emits records for each query together. The implementation streams
    one query group at a time, avoiding a read-sized in-memory dictionary.

    Args:
        paf_path: Compressed or plain minimap2 PAF output.
        reference_fasta: FASTA used to build the minimap2 index.
        taxon_report_path: Output taxon-level TSV.
        mapping_summary_path: Output one-row mapping summary TSV.
        sample_id: Logical sample identifier.
        minimum_mapq: Minimum mapping quality retained.
        minimum_alignment: Minimum PAF alignment block length retained.

    Raises:
        ValueError: If thresholds, PAF rows or reference identifiers are invalid.
    """
    if minimum_mapq < 0:
        raise ValueError("minimum_mapq must be non-negative")
    if minimum_alignment <= 0:
        raise ValueError("minimum_alignment must be positive")
    references = parse_reference_fasta(path=reference_fasta)
    counters: dict[str, dict[str, object]] = {}
    retained_alignment_count = 0
    mapped_read_count = 0

    for hits in _query_groups(
        hits=(
            hit
            for hit in parse_paf(path=paf_path)
            if hit.mapq >= minimum_mapq and hit.block_length >= minimum_alignment
        )
    ):
        mapped_read_count += 1
        retained_alignment_count += len(hits)
        for hit in hits:
            record = _reference_for_hit(hit=hit, references=references, paf_path=paf_path)
            counter = _counter(counters=counters, record=record)
            counter["alignment_count"] = int(counter["alignment_count"]) + 1
            cast_references = counter["reference_names"]
            assert isinstance(cast_references, set)
            cast_references.add(record.reference_name)

        best_score = max(hit.score for hit in hits)
        best_hits = [hit for hit in hits if hit.score == best_score]
        best_taxa = {
            _reference_for_hit(
                hit=hit,
                references=references,
                paf_path=paf_path,
            ).tax_id
            for hit in best_hits
        }
        for tax_id in best_taxa:
            counter = counters[tax_id]
            field = "best_read_count" if len(best_taxa) == 1 else "ambiguous_best_read_count"
            counter[field] = int(counter[field]) + 1

    rows = []
    for tax_id, values in sorted(counters.items()):
        reference_names = values["reference_names"]
        assert isinstance(reference_names, set)
        rows.append(
            {
                "sample_id": sample_id,
                "method": "minimap2",
                "tax_id": tax_id,
                "taxon_name": values["taxon_name"],
                "best_read_count": values["best_read_count"],
                "ambiguous_best_read_count": values["ambiguous_best_read_count"],
                "alignment_count": values["alignment_count"],
                "reference_count": len(reference_names),
                "min_mapq": minimum_mapq,
                "min_alignment": minimum_alignment,
            }
        )
    _write_tsv(
        path=taxon_report_path,
        rows=rows,
        fieldnames=(
            "sample_id",
            "method",
            "tax_id",
            "taxon_name",
            "best_read_count",
            "ambiguous_best_read_count",
            "alignment_count",
            "reference_count",
            "min_mapq",
            "min_alignment",
        ),
    )
    _write_tsv(
        path=mapping_summary_path,
        rows=[
            {
                "sample_id": sample_id,
                "method": "minimap2",
                "mapped_read_count": mapped_read_count,
                "retained_alignment_count": retained_alignment_count,
                "reported_taxon_count": len(rows),
                "min_mapq": minimum_mapq,
                "min_alignment": minimum_alignment,
            }
        ],
        fieldnames=(
            "sample_id",
            "method",
            "mapped_read_count",
            "retained_alignment_count",
            "reported_taxon_count",
            "min_mapq",
            "min_alignment",
        ),
    )


def _query_groups(*, hits: Iterable[PafHit]) -> Iterator[list[PafHit]]:
    """Yield consecutive PAF records grouped by query identifier."""
    group: list[PafHit] = []
    current_query: str | None = None
    for hit in hits:
        if current_query is not None and hit.query_name != current_query:
            yield group
            group = []
        group.append(hit)
        current_query = hit.query_name
    if group:
        yield group


def _reference_for_hit(
    *,
    hit: PafHit,
    references: dict[str, ReferenceRecord],
    paf_path: Path,
) -> ReferenceRecord:
    """Return reference metadata or fail for a mismatched FASTA/index pair."""
    try:
        return references[hit.target_name]
    except KeyError as error:
        raise ValueError(
            f"PAF target {hit.target_name!r} is absent from the configured FASTA; "
            f"the minimap2 index and reference may not match: {paf_path}"
        ) from error


def _counter(
    *,
    counters: dict[str, dict[str, object]],
    record: ReferenceRecord,
) -> dict[str, object]:
    """Return a mutable counter initialised for one reported taxon."""
    return counters.setdefault(
        record.tax_id,
        {
            "taxon_name": record.taxon_name,
            "best_read_count": 0,
            "ambiguous_best_read_count": 0,
            "alignment_count": 0,
            "reference_names": set(),
        },
    )


def _tax_id(*, header: str) -> str:
    """Extract a taxon identifier, retaining unknown references separately."""
    match = TAXID_PATTERN.search(header)
    if match:
        return match.group(1)
    return f"unknown:{header.split()[0]}"


def _taxon_name(*, reference_name: str, header: str, tax_id: str) -> str:
    """Derive a conservative display label without inventing taxonomy."""
    species_match = SPECIES_TOKEN_PATTERN.search(reference_name)
    if species_match:
        return species_match.group(1).replace("_", " ")
    description = header[len(reference_name) :].strip()
    if description and "kraken:taxid" not in description:
        return description
    if not tax_id.startswith("unknown:"):
        return f"taxid_{tax_id}"
    return reference_name


def _open_text(*, path: Path) -> TextIO:
    """Open a plain or gzip-compressed text file."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _write_tsv(
    *, path: Path, rows: Iterable[dict[str, object]], fieldnames: tuple[str, ...]
) -> None:
    """Write a deterministic tab-separated table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
