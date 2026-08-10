"""Controlled-reference minimap2 parsing and taxon-level summarisation."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

from nanopore_realdata.reference import species_name_from_header


TAXID_PATTERN = re.compile(r"kraken:taxid\|([^|,\s]+)")


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
            inferred_species = species_name_from_header(header=header)
            tax_id = _tax_id(header=header, inferred_species=inferred_species)
            records[reference_name] = ReferenceRecord(
                reference_name=reference_name,
                tax_id=tax_id,
                taxon_name=_taxon_name(
                    reference_name=reference_name,
                    header=header,
                    tax_id=tax_id,
                    inferred_species=inferred_species,
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
    input_read_count: int,
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
        input_read_count: Validated number of reads supplied to minimap2.

    Raises:
        ValueError: If thresholds, PAF rows or reference identifiers are invalid.
    """
    if minimum_mapq < 0:
        raise ValueError("minimum_mapq must be non-negative")
    if minimum_alignment <= 0:
        raise ValueError("minimum_alignment must be positive")
    if input_read_count < 0:
        raise ValueError("input_read_count must be non-negative")
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
        if mapped_read_count > input_read_count:
            raise ValueError(
                "minimap2 mapped-read groups exceed the validated input-read count; "
                "the index may be multipart or the PAF may contain duplicated query blocks"
            )
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
                "input_read_count": input_read_count,
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
            "input_read_count",
            "min_mapq",
            "min_alignment",
        ),
    )


def _query_groups(*, hits: Iterable[PafHit]) -> Iterator[list[PafHit]]:
    """Yield query groups while rejecting repeated non-consecutive queries.

    Minimap2 emits the same read block again for each part of a multipart
    index. Tracking completed query names turns that otherwise silent count
    inflation into a hard validation failure.
    """
    group: list[PafHit] = []
    current_query: str | None = None
    completed_queries: set[str] = set()
    for hit in hits:
        if current_query is not None and hit.query_name != current_query:
            yield group
            completed_queries.add(current_query)
            group = []
        if hit.query_name in completed_queries:
            raise ValueError(
                "PAF query appears in more than one non-consecutive block: "
                f"{hit.query_name}; multipart index output is not permitted"
            )
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


def _tax_id(*, header: str, inferred_species: str | None) -> str:
    """Extract a taxon identifier, retaining unknown references separately."""
    match = TAXID_PATTERN.search(header)
    if match:
        return match.group(1)
    if inferred_species:
        return "name:" + inferred_species.replace(" ", "_")
    return f"unknown:{header.split()[0]}"


def _taxon_name(
    *,
    reference_name: str,
    header: str,
    tax_id: str,
    inferred_species: str | None,
) -> str:
    """Derive a conservative display label without inventing taxonomy."""
    if inferred_species:
        return inferred_species
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
