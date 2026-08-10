"""Controlled minimap2 reference construction and index-log validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import re
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence, TextIO

from nanopore_realdata.pcr import canonical_species_name
from nanopore_realdata.runtime import sha256_file


SAFE_TOKEN = re.compile(r"[^A-Za-z0-9._-]+")
TAXON_NAME_TOKEN = re.compile(r"(?:^|\s)taxon_name=([^\s]+)")
PLASMODIUM_NAME_PATTERNS = (
    re.compile(
        r"(?<![A-Za-z])Plasmodium[\s_.|=:-]+([a-z][A-Za-z0-9-]*)",
        flags=re.IGNORECASE,
    ),
    re.compile(
        r"(?<![A-Za-z])Plas[\s_.|=:-]+([a-z][A-Za-z0-9-]*)",
        flags=re.IGNORECASE,
    ),
    re.compile(
        r"(?<![A-Za-z])P[._-]+([a-z][A-Za-z0-9-]*)",
        flags=re.IGNORECASE,
    ),
)
NON_SPECIES_EPITHETS = frozenset({"sp", "spp", "species", "strain", "isolate", "genome"})
INDEX_PART_PATTERN = re.compile(r"\[M::mm_idx_stat\].*\btotal length:\s*([0-9]+)")
IUPAC_DNA = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
VALID_GENOME_ROLES = frozenset(
    {
        "target_species",
        "target_clade_member",
        "target",
        "near_neighbour",
        "near_neighbor",
        "outgroup",
        "apicomplexan_outgroup",
        "distant_outgroup",
        "host",
        "host_or_background",
        "background_host",
        "host_background",
        "background",
        "background_pathogen",
        "environmental_background",
        "non_target",
        "downloaded",
        "exclude",
    }
)


@dataclass(frozen=True)
class GenomeSource:
    """One genome FASTA selected for the controlled reference."""

    genome_fasta: Path
    species_name: str
    taxid: str
    assembly_accession: str
    role: str
    clade: str


@dataclass(frozen=True)
class ReferenceBuildStats:
    """Summary statistics for a generated controlled FASTA."""

    genome_count: int
    reference_record_count: int
    total_bases: int
    reference_sha256: str
    manifest_sha256: str


def load_genome_sources(*, config_path: Path) -> tuple[GenomeSource, ...]:
    """Load the KmerSutra genome configuration used as reference provenance.

    Args:
        config_path: KmerSutra ``kmersutra_genome_config.tsv`` path.

    Returns:
        Validated source genomes in table order.

    Raises:
        FileNotFoundError: If the table or a selected FASTA is unavailable.
        ValueError: If required metadata are missing or duplicated.
    """
    resolved = config_path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"Genome configuration is missing or empty: {resolved}")
    required = {
        "genome_fasta",
        "species_name",
        "taxid",
        "assembly_accession",
        "role",
        "clade",
    }
    sources: list[GenomeSource] = []
    seen_paths: set[Path] = set()
    taxid_species: dict[str, str] = {}
    species_taxid: dict[str, str] = {}
    with resolved.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValueError(f"Genome configuration is missing columns: {missing}")
        for row_number, row in enumerate(reader, start=2):
            role = (row.get("role") or "unspecified").strip()
            if role not in VALID_GENOME_ROLES:
                raise ValueError(f"Unsupported genome role on row {row_number}: {role!r}")
            if role == "exclude":
                continue
            fasta_text = (row.get("genome_fasta") or "").strip()
            if not fasta_text:
                raise ValueError(f"Blank genome_fasta on row {row_number}")
            fasta = Path(fasta_text).expanduser()
            if not fasta.is_absolute():
                fasta = resolved.parent / fasta
            fasta = fasta.resolve()
            if not fasta.is_file() or fasta.stat().st_size == 0:
                raise FileNotFoundError(
                    f"Genome FASTA is missing or empty on row {row_number}: {fasta}"
                )
            if fasta in seen_paths:
                raise ValueError(f"Genome FASTA is listed more than once: {fasta}")
            seen_paths.add(fasta)
            species = canonical_species_name(value=(row.get("species_name") or ""))
            taxid = (row.get("taxid") or "").strip()
            if not species or not taxid:
                raise ValueError(f"species_name and taxid are required on genome row {row_number}")
            if not taxid.isdecimal():
                raise ValueError(
                    f"taxid must contain digits only on genome row {row_number}: {taxid!r}"
                )
            prior_species = taxid_species.setdefault(taxid, species)
            if prior_species.casefold() != species.casefold():
                raise ValueError(
                    f"taxid {taxid} maps to conflicting species names: "
                    f"{prior_species!r} and {species!r}"
                )
            species_key = species.casefold()
            prior_taxid = species_taxid.setdefault(species_key, taxid)
            if prior_taxid != taxid:
                raise ValueError(
                    f"Species {species!r} maps to conflicting taxids: {prior_taxid} and {taxid}"
                )
            sources.append(
                GenomeSource(
                    genome_fasta=fasta,
                    species_name=species,
                    taxid=taxid,
                    assembly_accession=(row.get("assembly_accession") or fasta.stem).strip(),
                    role=role,
                    clade=(row.get("clade") or "unspecified").strip(),
                )
            )
    if not sources:
        raise ValueError(f"Genome configuration contains no data rows: {resolved}")
    return tuple(sources)


def validate_required_species(
    *,
    sources: Sequence[GenomeSource],
    required_species: Sequence[str],
) -> None:
    """Require every PCR-expected species in the controlled source panel.

    Args:
        sources: Validated genome sources.
        required_species: Canonical PCR-expected species names.

    Raises:
        ValueError: If any expected species has no source genome.
    """
    available = {source.species_name.casefold() for source in sources}
    missing = [species for species in required_species if species.casefold() not in available]
    if missing:
        raise ValueError(
            "Controlled minimap2 sources lack PCR-expected species: " + "; ".join(missing)
        )


def validate_required_reference_species(
    *,
    reference_fasta: Path,
    required_species: Sequence[str],
) -> tuple[str, ...]:
    """Require PCR-expected species labels in a prebuilt reference FASTA.

    The reduced masked Plasmodium reference predates the controlled-reference
    header format. Its headers can therefore use full names (``Plasmodium
    inui``), filename-derived labels (``Plas_inui``), or abbreviations
    (``P.inui``). This check deliberately uses only explicit header evidence;
    an accession-only reference cannot support auditable species-level PCR
    comparison and is rejected for the minimap2 branch.

    Args:
        reference_fasta: Existing plain or gzip-compressed FASTA.
        required_species: Canonical PCR-expected species names.

    Returns:
        Unique species labels found in FASTA-header order.

    Raises:
        FileNotFoundError: If the FASTA is missing or empty.
        ValueError: If the FASTA has no headers or lacks a required species.
    """
    resolved = reference_fasta.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FileNotFoundError(f"minimap2 reference is missing or empty: {resolved}")
    inventory: list[str] = []
    seen: set[str] = set()
    header_count = 0
    with _open_text(path=resolved) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header_count += 1
            species = species_name_from_header(header=line[1:].strip())
            if species is None:
                continue
            folded = species.casefold()
            if folded not in seen:
                seen.add(folded)
                inventory.append(species)
    if header_count == 0:
        raise ValueError(f"Classification reference FASTA contains no records: {resolved}")
    missing = [species for species in required_species if species.casefold() not in seen]
    if missing:
        shown = "; ".join(inventory[:20]) if inventory else "none with auditable species labels"
        raise ValueError(
            "Prebuilt minimap2 reference lacks PCR-expected species labels: "
            + "; ".join(missing)
            + f". Parsed reference species: {shown}"
        )
    return tuple(inventory)


def build_controlled_reference(
    *,
    sources: Sequence[GenomeSource],
    output_fasta: Path,
    output_manifest: Path,
    maximum_reference_bases: int,
) -> ReferenceBuildStats:
    """Build a taxon-labelled FASTA from an explicit genome configuration.

    Args:
        sources: Ordered source genome records.
        output_fasta: Destination uncompressed FASTA.
        output_manifest: Destination sequence-level TSV manifest.
        maximum_reference_bases: Hard upper bound preventing runaway indexes.

    Returns:
        Reference counts and durable checksums.

    Raises:
        ValueError: If FASTA content is malformed or exceeds the hard bound.
    """
    if not sources:
        raise ValueError("At least one genome source is required")
    if maximum_reference_bases <= 0:
        raise ValueError("maximum_reference_bases must be positive")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    fasta_temporary = output_fasta.with_name(
        f".{output_fasta.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    manifest_temporary = output_manifest.with_name(
        f".{output_manifest.name}.partial.{os.getpid()}.{uuid.uuid4().hex}"
    )
    total_bases = 0
    record_count = 0
    try:
        with (
            fasta_temporary.open("w", encoding="utf-8", newline="\n") as fasta_handle,
            manifest_temporary.open("w", encoding="utf-8", newline="") as manifest_handle,
        ):
            fields = (
                "reference_id",
                "species_name",
                "taxid",
                "assembly_accession",
                "role",
                "clade",
                "source_fasta",
                "source_fasta_sha256",
                "source_record_id",
                "sequence_length",
            )
            writer = csv.DictWriter(
                manifest_handle,
                fieldnames=list(fields),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            seen_reference_ids: set[str] = set()
            for genome_index, source in enumerate(sources, start=1):
                source_digest = sha256_file(path=source.genome_fasta)
                for sequence_index, (source_id, sequence) in enumerate(
                    _read_fasta(path=source.genome_fasta),
                    start=1,
                ):
                    total_bases += len(sequence)
                    if total_bases > maximum_reference_bases:
                        raise ValueError(
                            "Controlled minimap2 reference exceeds the configured hard "
                            f"limit of {maximum_reference_bases} bases"
                        )
                    record_count += 1
                    reference_id = _reference_id(
                        source=source,
                        source_id=source_id,
                        genome_index=genome_index,
                        sequence_index=sequence_index,
                    )
                    if reference_id in seen_reference_ids:
                        raise ValueError(f"Generated duplicate reference ID: {reference_id}")
                    seen_reference_ids.add(reference_id)
                    taxon_token = source.species_name.replace(" ", "_")
                    fasta_handle.write(
                        f">{reference_id} kraken:taxid|{source.taxid}| "
                        f"taxon_name={taxon_token} role={_token(source.role)} "
                        f"clade={_token(source.clade)}\n"
                    )
                    for offset in range(0, len(sequence), 80):
                        fasta_handle.write(sequence[offset : offset + 80] + "\n")
                    writer.writerow(
                        {
                            "reference_id": reference_id,
                            "species_name": source.species_name,
                            "taxid": source.taxid,
                            "assembly_accession": source.assembly_accession,
                            "role": source.role,
                            "clade": source.clade,
                            "source_fasta": str(source.genome_fasta),
                            "source_fasta_sha256": source_digest,
                            "source_record_id": source_id,
                            "sequence_length": len(sequence),
                        }
                    )
        if record_count == 0 or total_bases == 0:
            raise ValueError("Controlled minimap2 reference contains no sequence")
        os.replace(fasta_temporary, output_fasta)
        os.replace(manifest_temporary, output_manifest)
    finally:
        fasta_temporary.unlink(missing_ok=True)
        manifest_temporary.unlink(missing_ok=True)
    return ReferenceBuildStats(
        genome_count=len(sources),
        reference_record_count=record_count,
        total_bases=total_bases,
        reference_sha256=sha256_file(path=output_fasta),
        manifest_sha256=sha256_file(path=output_manifest),
    )


def validate_single_part_index_log(
    *,
    log_path: Path,
    maximum_reference_bases: int,
) -> dict[str, int]:
    """Validate that minimap2 built exactly one bounded index part.

    Args:
        log_path: Persistent minimap2 index build log.
        maximum_reference_bases: Hard maximum expected reference length.

    Returns:
        Parsed part count and indexed base total.

    Raises:
        ValueError: If the log is absent, multipart or inconsistent.
    """
    if not log_path.is_file() or log_path.stat().st_size == 0:
        raise ValueError(f"minimap2 index log is missing or empty: {log_path}")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    lengths = [int(value) for value in INDEX_PART_PATTERN.findall(text)]
    if len(lengths) != 1:
        raise ValueError(
            "minimap2 index must contain exactly one part; "
            f"the build log recorded {len(lengths)} parts"
        )
    if lengths[0] <= 0 or lengths[0] > maximum_reference_bases:
        raise ValueError(f"minimap2 indexed base count is outside the allowed range: {lengths[0]}")
    return {"index_part_count": 1, "indexed_bases": lengths[0]}


def taxon_name_from_header(*, header: str) -> str | None:
    """Read the explicit ``taxon_name`` token emitted by this module.

    Args:
        header: FASTA header without the leading ``>``.

    Returns:
        Canonical taxon name, or ``None`` when the token is absent.
    """
    match = TAXON_NAME_TOKEN.search(header)
    if match is None:
        return None
    return canonical_species_name(value=match.group(1))


def species_name_from_header(*, header: str) -> str | None:
    """Infer an auditable species label from established FASTA header forms.

    Args:
        header: FASTA header without the leading ``>``.

    Returns:
        Explicit ``taxon_name`` metadata or a conservatively parsed
        Plasmodium binomial; otherwise ``None``.
    """
    explicit = taxon_name_from_header(header=header)
    if explicit:
        return explicit
    for pattern in PLASMODIUM_NAME_PATTERNS:
        match = pattern.search(header)
        if match is None:
            continue
        epithet = match.group(1).strip("-_.").casefold()
        if epithet and epithet not in NON_SPECIES_EPITHETS:
            return f"Plasmodium {epithet}"
    return None


def _read_fasta(*, path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    sequence_parts: list[str] = []
    with _open_text(path=path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header is not None:
                    yield (
                        header,
                        _validated_sequence(
                            parts=sequence_parts,
                            path=path,
                            header=header,
                        ),
                    )
                header = stripped[1:].split()[0]
                if not header:
                    raise ValueError(f"Blank FASTA header at line {line_number}: {path}")
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence precedes the first FASTA header: {path}")
                sequence_parts.append(stripped)
    if header is not None:
        yield header, _validated_sequence(parts=sequence_parts, path=path, header=header)


def _validated_sequence(*, parts: Sequence[str], path: Path, header: str) -> str:
    sequence = "".join(parts)
    if not sequence:
        raise ValueError(f"Empty FASTA sequence for {header!r}: {path}")
    invalid = sorted(set(sequence).difference(IUPAC_DNA))
    if invalid:
        raise ValueError(f"Invalid FASTA character(s) for {header!r}: {''.join(invalid)}")
    return sequence.upper().replace("U", "T")


def _reference_id(
    *,
    source: GenomeSource,
    source_id: str,
    genome_index: int,
    sequence_index: int,
) -> str:
    assembly = _token(source.assembly_accession or f"genome_{genome_index}")
    original = _token(source_id)
    digest = hashlib.sha1(source_id.encode("utf-8")).hexdigest()[:10]
    return f"{assembly}.{genome_index:04d}.{sequence_index:06d}.{original[:48]}.{digest}"


def _token(value: str) -> str:
    cleaned = SAFE_TOKEN.sub("_", value.strip()).strip("._-")
    return cleaned or "unknown"


def _open_text(*, path: Path) -> TextIO:
    if str(path).casefold().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")
