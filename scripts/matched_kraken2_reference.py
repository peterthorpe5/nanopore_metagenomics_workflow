#!/usr/bin/env python3
"""Prepare and validate a held-out, source-matched Kraken2 database."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO


LOGGER = logging.getLogger("matched_kraken2_reference")
REQUIRED_GENOME_COLUMNS = {
    "genome_fasta",
    "species_name",
    "taxid",
    "assembly_accession",
    "role",
    "source",
}
REQUIRED_TRUTH_COLUMNS = {
    "species_name",
    "accepted_species_names",
    "expected",
    "ncbi_taxid",
    "truth_assembly_accessions",
    "truth_sequence_accessions",
}
ACCESSION_PATTERN = re.compile(r"[A-Za-z]{1,6}_[A-Za-z0-9]+(?:\.\d+)?")
SAFE_HEADER_PATTERN = re.compile(r"[^A-Za-z0-9_.:-]+")
KRAKEN_DATABASE_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")


@dataclass(frozen=True)
class GenomeRecord:
    """One source assembly used in the matched database."""

    genome_fasta: Path
    species_name: str
    taxid: int
    assembly_accession: str
    role: str
    source: str


@dataclass(frozen=True)
class TruthRecord:
    """One expected organism and its accepted identifiers."""

    species_name: str
    accepted_species_names: tuple[str, ...]
    ncbi_taxid: int
    truth_assembly_accessions: tuple[str, ...]
    truth_sequence_accessions: tuple[str, ...]

    @property
    def all_species_names(self) -> tuple[str, ...]:
        """Return the published name and accepted synonyms."""
        return (self.species_name, *self.accepted_species_names)


def _utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _split_values(value: str) -> tuple[str, ...]:
    """Split semicolon-delimited manifest values."""
    return tuple(item.strip() for item in value.split(";") if item.strip())


def _normalise_accession(value: str) -> str:
    """Normalise an accession for leakage comparison."""
    return value.strip().upper().split(".", maxsplit=1)[0]


def _parse_positive_int(value: str, *, field: str, row_number: int) -> int:
    """Parse a positive integer with a useful TSV error."""
    text = value.strip()
    if not text.isdigit() or int(text) <= 0:
        raise ValueError(f"Invalid {field} on row {row_number}: {value!r}")
    return int(text)


def _resolve_input_path(value: str, *, base: Path, row_number: int) -> Path:
    """Resolve a genome path and require a non-empty file."""
    expanded = os.path.expandvars(os.path.expanduser(value.strip()))
    if not expanded:
        raise ValueError(f"Blank genome_fasta on row {row_number}")
    path = Path(expanded)
    if not path.is_absolute():
        path = base / path
    resolved = path.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"Genome FASTA is missing or empty on row {row_number}: {resolved}")
    return resolved


def load_genome_config(path: Path) -> tuple[GenomeRecord, ...]:
    """Load and strictly validate the KmerSutra genome configuration.

    Args:
        path: Tab-separated KmerSutra genome configuration.

    Returns:
        Genome records in source order.

    Raises:
        ValueError: If required fields, identifiers or uniqueness checks fail.
        FileNotFoundError: If the configuration or a FASTA is unavailable.
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Genome configuration does not exist: {resolved}")
    records: list[GenomeRecord] = []
    seen_paths: set[Path] = set()
    seen_assemblies: set[str] = set()
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = REQUIRED_GENOME_COLUMNS - fields
        if missing:
            raise ValueError(
                "Genome configuration is missing columns: " + ", ".join(sorted(missing))
            )
        for row_number, row in enumerate(reader, start=2):
            species_name = (row.get("species_name") or "").strip()
            assembly = (row.get("assembly_accession") or "").strip()
            if not species_name:
                raise ValueError(f"Blank species_name on row {row_number}")
            if not assembly:
                raise ValueError(f"Blank assembly_accession on row {row_number}")
            genome_path = _resolve_input_path(
                row.get("genome_fasta") or "",
                base=resolved.parent,
                row_number=row_number,
            )
            if genome_path in seen_paths:
                raise ValueError(f"Genome FASTA is listed more than once: {genome_path}")
            normalised_assembly = _normalise_accession(assembly)
            if normalised_assembly in seen_assemblies:
                raise ValueError(f"Assembly accession is listed more than once: {assembly}")
            seen_paths.add(genome_path)
            seen_assemblies.add(normalised_assembly)
            records.append(
                GenomeRecord(
                    genome_fasta=genome_path,
                    species_name=species_name,
                    taxid=_parse_positive_int(
                        row.get("taxid") or "",
                        field="taxid",
                        row_number=row_number,
                    ),
                    assembly_accession=assembly,
                    role=(row.get("role") or "").strip(),
                    source=(row.get("source") or "").strip(),
                )
            )
    if not records:
        raise ValueError("Genome configuration contains no data rows")
    return tuple(records)


def load_truth_manifest(path: Path) -> tuple[TruthRecord, ...]:
    """Load expected organisms from the locked ATCC truth manifest."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Truth manifest does not exist: {resolved}")
    records: list[TruthRecord] = []
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = REQUIRED_TRUTH_COLUMNS - fields
        if missing:
            raise ValueError(f"Truth manifest is missing columns: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            expected = (row.get("expected") or "").strip().casefold()
            if expected not in {"true", "false"}:
                raise ValueError(f"Invalid expected value on row {row_number}: {expected!r}")
            if expected == "false":
                continue
            species_name = (row.get("species_name") or "").strip()
            if not species_name:
                raise ValueError(f"Blank species_name on truth row {row_number}")
            records.append(
                TruthRecord(
                    species_name=species_name,
                    accepted_species_names=_split_values(
                        row.get("accepted_species_names") or ""
                    ),
                    ncbi_taxid=_parse_positive_int(
                        row.get("ncbi_taxid") or "",
                        field="ncbi_taxid",
                        row_number=row_number,
                    ),
                    truth_assembly_accessions=_split_values(
                        row.get("truth_assembly_accessions") or ""
                    ),
                    truth_sequence_accessions=_split_values(
                        row.get("truth_sequence_accessions") or ""
                    ),
                )
            )
    if not records:
        raise ValueError("Truth manifest contains no expected organisms")
    return tuple(records)


def validate_reference_gate(path: Path) -> None:
    """Require a completed KmerSutra reference gate with no failing checks."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"Reference gate summary is missing or empty: {resolved}")
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    tokens = {cell.strip().casefold() for row in rows for cell in row if cell.strip()}
    if "fail" in tokens or "failed" in tokens:
        raise ValueError(f"Reference gate contains a failing status: {resolved}")
    if "pass" not in tokens and "passed" not in tokens:
        raise ValueError(f"Reference gate does not contain an explicit PASS status: {resolved}")


def _open_fasta(path: Path) -> TextIO:
    """Open a plain or gzip-compressed FASTA as text."""
    if path.name.casefold().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _header_accessions(header: str) -> set[str]:
    """Extract normalised accession-like identifiers from a FASTA header."""
    return {_normalise_accession(match.group(0)) for match in ACCESSION_PATTERN.finditer(header)}


def _safe_identifier(value: str) -> str:
    """Make a bounded, deterministic FASTA header component."""
    cleaned = SAFE_HEADER_PATTERN.sub("_", value.strip()).strip("_")
    return (cleaned or "record")[:120]


def _sha256(path: Path) -> str:
    """Calculate a file SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    """Write JSON atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _validate_expected_counts(
    genomes: tuple[GenomeRecord, ...],
    truth: tuple[TruthRecord, ...],
    *,
    expected_genome_count: int,
    expected_source_species_count: int,
    expected_truth_species_count: int,
) -> tuple[str, ...]:
    """Validate the locked benchmark dimensions and truth coverage."""
    if len(genomes) != expected_genome_count:
        raise ValueError(f"Expected {expected_genome_count} genomes, observed {len(genomes)}")
    source_names = {record.species_name.casefold() for record in genomes}
    if len(source_names) != expected_source_species_count:
        raise ValueError(
            f"Expected {expected_source_species_count} source species, observed {len(source_names)}"
        )
    if len(truth) != expected_truth_species_count:
        raise ValueError(
            f"Expected {expected_truth_species_count} truth species, observed {len(truth)}"
        )
    missing = tuple(
        record.species_name
        for record in truth
        if not any(name.casefold() in source_names for name in record.all_species_names)
    )
    if missing:
        raise ValueError(f"Truth species absent from the source genomes: {', '.join(missing)}")
    return missing


def prepare_library(
    *,
    genome_config: Path,
    truth_manifest: Path,
    gate_summary: Path,
    output_fasta: Path,
    output_manifest: Path,
    output_summary: Path,
    expected_genome_count: int,
    expected_source_species_count: int,
    expected_truth_species_count: int,
) -> dict[str, object]:
    """Create one taxid-annotated FASTA after locked benchmark checks."""
    validate_reference_gate(gate_summary)
    genomes = load_genome_config(genome_config)
    truth = load_truth_manifest(truth_manifest)
    _validate_expected_counts(
        genomes,
        truth,
        expected_genome_count=expected_genome_count,
        expected_source_species_count=expected_source_species_count,
        expected_truth_species_count=expected_truth_species_count,
    )
    held_out_assemblies = {
        _normalise_accession(value)
        for record in truth
        for value in record.truth_assembly_accessions
    }
    held_out_sequences = {
        _normalise_accession(value)
        for record in truth
        for value in record.truth_sequence_accessions
    }
    leaked_assemblies = sorted(
        {
            record.assembly_accession
            for record in genomes
            if _normalise_accession(record.assembly_accession) in held_out_assemblies
        }
    )
    if leaked_assemblies:
        raise ValueError(
            "Held-out truth assemblies occur in the source configuration: "
            + ", ".join(leaked_assemblies)
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    fasta_temporary = output_fasta.with_name(f".{output_fasta.name}.{os.getpid()}.tmp")
    manifest_temporary = output_manifest.with_name(f".{output_manifest.name}.{os.getpid()}.tmp")
    record_count = 0
    base_count = 0
    try:
        with fasta_temporary.open("w", encoding="utf-8", newline="\n") as fasta_out, \
                manifest_temporary.open("w", encoding="utf-8", newline="") as manifest_out:
            fieldnames = [
                "genome_fasta",
                "species_name",
                "taxid",
                "assembly_accession",
                "role",
                "source",
                "source_size_bytes",
                "sequence_record_count",
                "sequence_base_count",
            ]
            writer = csv.DictWriter(manifest_out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for genome_index, genome in enumerate(genomes, start=1):
                genome_records = 0
                genome_bases = 0
                current_header_seen = False
                with _open_fasta(genome.genome_fasta) as source_handle:
                    for line_number, raw_line in enumerate(source_handle, start=1):
                        line = raw_line.strip()
                        if not line:
                            continue
                        if line.startswith(">"):
                            header = line[1:].strip()
                            if not header:
                                raise ValueError(
                                    f"Blank FASTA header in {genome.genome_fasta}:{line_number}"
                                )
                            leaked = _header_accessions(header) & held_out_sequences
                            if leaked:
                                raise ValueError(
                                    "Held-out truth sequence accession found in source FASTA "
                                    f"{genome.genome_fasta}:{line_number}: "
                                    f"{', '.join(sorted(leaked))}"
                                )
                            genome_records += 1
                            record_count += 1
                            current_header_seen = True
                            original_id = header.split(maxsplit=1)[0]
                            assembly_identifier = _safe_identifier(
                                genome.assembly_accession
                            )
                            record_identifier = _safe_identifier(original_id)
                            fasta_out.write(
                                f">kraken:taxid|{genome.taxid}|{assembly_identifier}"
                                f"|g{genome_index}|r{genome_records}|{record_identifier}\n"
                            )
                            continue
                        if not current_header_seen:
                            raise ValueError(
                                f"Sequence occurs before a FASTA header in "
                                f"{genome.genome_fasta}:{line_number}"
                            )
                        sequence = "".join(line.split()).upper()
                        if not sequence:
                            continue
                        genome_bases += len(sequence)
                        base_count += len(sequence)
                        fasta_out.write(sequence + "\n")
                if genome_records == 0 or genome_bases == 0:
                    raise ValueError(f"FASTA contains no sequence records: {genome.genome_fasta}")
                writer.writerow(
                    {
                        "genome_fasta": str(genome.genome_fasta),
                        "species_name": genome.species_name,
                        "taxid": genome.taxid,
                        "assembly_accession": genome.assembly_accession,
                        "role": genome.role,
                        "source": genome.source,
                        "source_size_bytes": genome.genome_fasta.stat().st_size,
                        "sequence_record_count": genome_records,
                        "sequence_base_count": genome_bases,
                    }
                )
        fasta_temporary.replace(output_fasta)
        manifest_temporary.replace(output_manifest)
    except Exception:
        fasta_temporary.unlink(missing_ok=True)
        manifest_temporary.unlink(missing_ok=True)
        raise

    summary: dict[str, object] = {
        "status": "prepared",
        "created_utc": _utc_now(),
        "genome_count": len(genomes),
        "source_species_count": len({record.species_name.casefold() for record in genomes}),
        "truth_species_count": len(truth),
        "truth_species_represented": len(truth),
        "held_out_assembly_leakage_count": 0,
        "held_out_sequence_leakage_count": 0,
        "sequence_record_count": record_count,
        "sequence_base_count": base_count,
        "genome_config": str(genome_config.expanduser().resolve()),
        "truth_manifest": str(truth_manifest.expanduser().resolve()),
        "gate_summary": str(gate_summary.expanduser().resolve()),
        "combined_fasta": str(output_fasta.resolve()),
        "combined_fasta_sha256": _sha256(output_fasta),
        "assembly_manifest": str(output_manifest.resolve()),
        "assembly_manifest_sha256": _sha256(output_manifest),
    }
    _atomic_json(output_summary, summary)
    return summary


def _parse_kraken_inspect(path: Path) -> tuple[set[int], int]:
    """Return taxids and row count from standard Kraken2 inspect output."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"Kraken2 inspect output is missing or empty: {resolved}")
    taxids: set[int] = set()
    row_count = 0
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"Malformed kraken2-inspect row {line_number}: {line!r}")
            taxid_text = fields[4].strip()
            if not taxid_text.isdigit():
                raise ValueError(
                    f"Invalid taxid in kraken2-inspect row {line_number}: {taxid_text!r}"
                )
            taxids.add(int(taxid_text))
            row_count += 1
    if not row_count:
        raise ValueError("Kraken2 inspect output contains no data rows")
    return taxids, row_count


def validate_database(
    *,
    database: Path,
    inspect_report: Path,
    truth_manifest: Path,
    preparation_summary: Path,
    output_summary: Path | None,
) -> dict[str, object]:
    """Validate database files and the representation of every truth taxid."""
    resolved_database = database.expanduser().resolve()
    if not resolved_database.is_dir():
        raise FileNotFoundError(
            f"Kraken2 database directory does not exist: {resolved_database}"
        )
    missing_database_files = [
        name
        for name in KRAKEN_DATABASE_FILES
        if not (resolved_database / name).is_file()
        or (resolved_database / name).stat().st_size == 0
    ]
    if missing_database_files:
        raise FileNotFoundError(
            "Kraken2 database is incomplete; missing: " + ", ".join(missing_database_files)
        )
    truth = load_truth_manifest(truth_manifest)
    inspect_taxids, inspect_rows = _parse_kraken_inspect(inspect_report)
    missing_truth = [record for record in truth if record.ncbi_taxid not in inspect_taxids]
    if missing_truth:
        details = ", ".join(
            f"{record.species_name} ({record.ncbi_taxid})" for record in missing_truth
        )
        raise ValueError(f"Expected truth taxids absent from matched Kraken2 database: {details}")
    with preparation_summary.expanduser().resolve().open("r", encoding="utf-8") as handle:
        preparation = json.load(handle)
    required_preparation = {
        "genome_count": 475,
        "source_species_count": 396,
        "truth_species_count": 20,
        "held_out_assembly_leakage_count": 0,
        "held_out_sequence_leakage_count": 0,
    }
    for key, expected in required_preparation.items():
        if preparation.get(key) != expected:
            raise ValueError(
                f"Preparation summary has {key}={preparation.get(key)!r}; expected {expected!r}"
            )
    summary: dict[str, object] = {
        "status": "complete",
        "validated_utc": _utc_now(),
        "database": str(resolved_database),
        "database_files": {
            name: (resolved_database / name).stat().st_size for name in KRAKEN_DATABASE_FILES
        },
        "kraken2_inspect": str(inspect_report.expanduser().resolve()),
        "kraken2_inspect_sha256": _sha256(inspect_report.expanduser().resolve()),
        "kraken2_inspect_row_count": inspect_rows,
        "truth_species_count": len(truth),
        "truth_taxids_represented": len(truth),
        "missing_truth_taxids": [],
        "preparation": preparation,
    }
    if output_summary is not None:
        _atomic_json(output_summary, summary)
    return summary


def build_parser() -> argparse.ArgumentParser:
    """Build the named-option-only command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--action",
        required=True,
        choices=("prepare-library", "validate-database"),
    )
    parser.add_argument("--genome-config", type=Path)
    parser.add_argument("--truth-manifest", type=Path, required=True)
    parser.add_argument("--gate-summary", type=Path)
    parser.add_argument("--output-fasta", type=Path)
    parser.add_argument("--output-manifest", type=Path)
    parser.add_argument("--output-summary", type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--inspect-report", type=Path)
    parser.add_argument("--preparation-summary", type=Path)
    parser.add_argument("--expected-genome-count", type=int, default=475)
    parser.add_argument("--expected-source-species-count", type=int, default=396)
    parser.add_argument("--expected-truth-species-count", type=int, default=20)
    parser.add_argument("--verbose", action="store_true")
    return parser


def _require_arguments(args: argparse.Namespace, names: Iterable[str]) -> None:
    """Require action-specific named arguments."""
    missing = [name for name in names if getattr(args, name) is None]
    if missing:
        options = ", ".join("--" + name.replace("_", "-") for name in missing)
        raise ValueError(f"Missing arguments for {args.action}: {options}")


def main(argv: list[str] | None = None) -> int:
    """Run the selected preparation or validation action."""
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
    )
    try:
        if args.action == "prepare-library":
            _require_arguments(
                args,
                (
                    "genome_config",
                    "gate_summary",
                    "output_fasta",
                    "output_manifest",
                    "output_summary",
                ),
            )
            summary = prepare_library(
                genome_config=args.genome_config,
                truth_manifest=args.truth_manifest,
                gate_summary=args.gate_summary,
                output_fasta=args.output_fasta,
                output_manifest=args.output_manifest,
                output_summary=args.output_summary,
                expected_genome_count=args.expected_genome_count,
                expected_source_species_count=args.expected_source_species_count,
                expected_truth_species_count=args.expected_truth_species_count,
            )
        else:
            _require_arguments(
                args,
                ("database", "inspect_report", "preparation_summary"),
            )
            summary = validate_database(
                database=args.database,
                inspect_report=args.inspect_report,
                truth_manifest=args.truth_manifest,
                preparation_summary=args.preparation_summary,
                output_summary=args.output_summary,
            )
    except (FileNotFoundError, OSError, ValueError, json.JSONDecodeError) as error:
        LOGGER.error("%s", error)
        return 2
    LOGGER.info("Action %s completed: %s", args.action, json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
