#!/usr/bin/env python3
"""Run and summarise the ATCC HiFi classifier operating-point sweep."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


LOGGER = logging.getLogger("atcc_hifi_classifier_sweep")
SUPPORTED_METHODS = ("kraken2", "metabuli")
READ_SUPPORT_THRESHOLDS = (1, 2, 10, 100)
SPECIES_RANKS = {"s", "species"}


@dataclass(frozen=True)
class OperatingPoint:
    """One predeclared classifier operating point."""

    method: str
    setting_id: str
    confidence: float | None
    min_score: float | None
    min_sp_score: float | None
    precise: int | None
    reuse_existing: bool
    setting_context: str
    description: str


@dataclass(frozen=True)
class ReportRow:
    """A normalised direct-assignment row from a classifier report."""

    percentage: float
    clade_reads: int
    direct_reads: int
    rank: str
    taxid: int
    name: str


def configure_logging(verbose: bool) -> None:
    """Configure deterministic stderr logging.

    Args:
        verbose: Emit debug-level messages when true.
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
    )


def utc_now() -> str:
    """Return the current UTC time in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat()


def parse_optional_float(value: str, field_name: str) -> float | None:
    """Parse an optional unit-interval floating-point value.

    Args:
        value: Raw TSV value.
        field_name: Field name used in validation errors.

    Returns:
        Parsed float or ``None`` for a blank field.

    Raises:
        ValueError: If the value is not numeric or is outside zero to one.
    """
    if not value.strip():
        return None
    try:
        parsed = float(value)
    except ValueError as error:
        raise ValueError(f"{field_name} must be numeric: {value!r}") from error
    if not 0.0 <= parsed <= 1.0:
        raise ValueError(f"{field_name} must be between 0 and 1: {parsed}")
    return parsed


def parse_boolean(value: str, field_name: str) -> bool:
    """Parse a strict true/false TSV value.

    Args:
        value: Raw TSV value.
        field_name: Field name used in validation errors.

    Returns:
        Parsed boolean.

    Raises:
        ValueError: If the field is not ``true`` or ``false``.
    """
    normalised = value.strip().lower()
    if normalised == "true":
        return True
    if normalised == "false":
        return False
    raise ValueError(f"{field_name} must be true or false: {value!r}")


def parse_optional_integer(value: str, field_name: str) -> int | None:
    """Parse an optional integer.

    Args:
        value: Raw TSV value.
        field_name: Field name used in validation errors.

    Returns:
        Parsed integer or ``None`` for a blank field.

    Raises:
        ValueError: If the field is not an integer.
    """
    if not value.strip():
        return None
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"{field_name} must be an integer: {value!r}") from error


def load_operating_points(path: Path) -> list[OperatingPoint]:
    """Load and validate a predeclared operating-point manifest.

    Args:
        path: Tab-separated manifest.

    Returns:
        Ordered operating points.

    Raises:
        ValueError: If the manifest is malformed or internally inconsistent.
    """
    required = {
        "method",
        "setting_id",
        "confidence",
        "min_score",
        "min_sp_score",
        "precise",
        "reuse_existing",
        "setting_context",
        "description",
    }
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Operating-point manifest is missing or empty: {path}")
    points: list[OperatingPoint] = []
    seen: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Operating-point manifest lacks columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            method = row["method"].strip().lower()
            setting_id = row["setting_id"].strip()
            if method not in SUPPORTED_METHODS:
                raise ValueError(f"Unsupported method on row {row_number}: {method}")
            if not setting_id or any(
                character not in "abcdefghijklmnopqrstuvwxyz0123456789_"
                for character in setting_id
            ):
                raise ValueError(f"Unsafe setting_id on row {row_number}: {setting_id!r}")
            key = (method, setting_id)
            if key in seen:
                raise ValueError(f"Duplicate operating point: {method}/{setting_id}")
            seen.add(key)
            point = OperatingPoint(
                method=method,
                setting_id=setting_id,
                confidence=parse_optional_float(row["confidence"], "confidence"),
                min_score=parse_optional_float(row["min_score"], "min_score"),
                min_sp_score=parse_optional_float(row["min_sp_score"], "min_sp_score"),
                precise=parse_optional_integer(row["precise"], "precise"),
                reuse_existing=parse_boolean(row["reuse_existing"], "reuse_existing"),
                setting_context=row["setting_context"].strip(),
                description=row["description"].strip(),
            )
            validate_operating_point(point)
            points.append(point)
    if not points:
        raise ValueError("Operating-point manifest contains no data rows")
    return points


def validate_operating_point(point: OperatingPoint) -> None:
    """Validate method-specific parameters for one operating point.

    Args:
        point: Operating point to validate.

    Raises:
        ValueError: If required values are missing or incompatible values are set.
    """
    if point.method == "kraken2":
        if point.confidence is None:
            raise ValueError(f"Kraken2 point lacks confidence: {point.setting_id}")
        if (
            point.min_score is not None
            or point.min_sp_score is not None
            or point.precise is not None
        ):
            raise ValueError(f"Kraken2 point contains Metabuli settings: {point.setting_id}")
    elif point.method == "metabuli":
        if point.confidence is not None:
            raise ValueError(f"Metabuli point contains Kraken2 confidence: {point.setting_id}")
        if point.precise not in {0, 1, 2}:
            raise ValueError(f"Metabuli precise preset must be 0, 1 or 2: {point.setting_id}")
        if point.precise == 0 and (point.min_score is None or point.min_sp_score is None):
            raise ValueError(f"Metabuli point lacks explicit score settings: {point.setting_id}")
        if point.precise in {1, 2} and (
            point.min_score is not None or point.min_sp_score is not None
        ):
            raise ValueError(
                f"Metabuli preset point must not override score settings: {point.setting_id}"
            )


def pending_points(points: Iterable[OperatingPoint], method: str) -> list[OperatingPoint]:
    """Return non-reused points for one method in manifest order.

    Args:
        points: All operating points.
        method: Classifier method.

    Returns:
        Points that require a new classification.
    """
    if method not in SUPPORTED_METHODS:
        raise ValueError(f"Unsupported method: {method}")
    return [point for point in points if point.method == method and not point.reuse_existing]


def require_file(path: Path, label: str) -> None:
    """Require a non-empty regular file.

    Args:
        path: File to validate.
        label: Human-readable label for errors.

    Raises:
        FileNotFoundError: If the file is missing or empty.
    """
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is missing or empty: {path}")


def require_directory(path: Path, label: str) -> None:
    """Require an existing directory.

    Args:
        path: Directory to validate.
        label: Human-readable label for errors.

    Raises:
        NotADirectoryError: If the directory is absent.
    """
    if not path.is_dir():
        raise NotADirectoryError(f"{label} is missing: {path}")


def build_classifier_command(
    point: OperatingPoint,
    *,
    input_fastq: Path,
    database: Path,
    staging_directory: Path,
    threads: int,
    max_ram_gb: int,
) -> list[str]:
    """Build a classifier command without invoking a shell.

    Args:
        point: Selected operating point.
        input_fastq: Classification-ready FASTQ file.
        database: Classifier database directory.
        staging_directory: Temporary output directory.
        threads: CPU thread count.
        max_ram_gb: Metabuli RAM limit in GiB.

    Returns:
        Command argument list.
    """
    if threads < 1:
        raise ValueError("threads must be at least one")
    if max_ram_gb < 1:
        raise ValueError("max_ram_gb must be at least one")
    if point.method == "kraken2":
        assert point.confidence is not None
        return [
            "kraken2",
            "--db",
            str(database),
            "--threads",
            str(threads),
            "--confidence",
            format(point.confidence, ".6g"),
            "--gzip-compressed",
            "--report",
            str(staging_directory / "report.tsv"),
            "--output",
            str(staging_directory / "classifications.tsv"),
            str(input_fastq),
        ]
    assert point.precise is not None
    command = [
        "metabuli",
        "classify",
        "--seq-mode",
        "3",
        str(input_fastq),
        str(database),
        str(staging_directory),
        point.setting_id,
        "--threads",
        str(threads),
        "--max-ram",
        str(max_ram_gb),
        "--precise",
        str(point.precise),
    ]
    if point.precise == 0:
        assert point.min_score is not None
        assert point.min_sp_score is not None
        command.extend(
            [
                "--min-score",
                format(point.min_score, ".6g"),
                "--min-sp-score",
                format(point.min_sp_score, ".6g"),
            ]
        )
    return command


def tool_version(executable: str) -> str:
    """Return a best-effort classifier version string.

    Args:
        executable: Executable name.

    Returns:
        First non-empty version line, or ``unavailable``.
    """
    attempts = ([executable, "--version"], [executable, "version"])
    for command in attempts:
        try:
            result = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        output = "\n".join((result.stdout, result.stderr)).strip()
        if output:
            return output.splitlines()[0].strip()
    return "unavailable"


def sha256_file(path: Path) -> str:
    """Calculate a SHA-256 digest for a file.

    Args:
        path: File to hash.

    Returns:
        Lower-case hexadecimal digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_file(path: Path) -> Path:
    """Compress a file deterministically and remove the source.

    Args:
        path: Uncompressed input file.

    Returns:
        Gzip output path.
    """
    output = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as source, output.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    path.unlink()
    return output


def normalise_metabuli_outputs(staging_directory: Path, setting_id: str) -> Path:
    """Normalise Metabuli output names and compress per-read classifications.

    Args:
        staging_directory: Metabuli output directory.
        setting_id: Job identifier supplied to Metabuli.

    Returns:
        Compressed classification path.

    Raises:
        RuntimeError: If required Metabuli outputs are absent or ambiguous.
    """
    report_candidates = list(staging_directory.glob(f"{setting_id}*_report.tsv"))
    classification_candidates = list(
        staging_directory.glob(f"{setting_id}*_classifications.tsv")
    )
    if len(report_candidates) != 1:
        raise RuntimeError(
            f"Expected one Metabuli report for {setting_id}; found {report_candidates}"
        )
    if len(classification_candidates) != 1:
        raise RuntimeError(
            "Expected one Metabuli classifications file for "
            f"{setting_id}; found {classification_candidates}"
        )
    report_candidates[0].replace(staging_directory / "report.tsv")
    classification_candidates[0].replace(staging_directory / "classifications.tsv")
    krona_candidates = list(staging_directory.glob(f"{setting_id}*_krona.html"))
    if len(krona_candidates) == 1:
        krona_candidates[0].replace(staging_directory / "krona.html")
    return gzip_file(staging_directory / "classifications.tsv")


def write_json(path: Path, payload: object) -> None:
    """Write stable, human-readable JSON.

    Args:
        path: Destination path.
        payload: JSON-serialisable content.
    """
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_operating_point(
    point: OperatingPoint,
    *,
    input_fastq: Path,
    database: Path,
    output_root: Path,
    threads: int,
    max_ram_gb: int,
) -> Path:
    """Run one classifier point and publish outputs atomically.

    Args:
        point: Operating point to run.
        input_fastq: Classification-ready reads.
        database: Classifier database.
        output_root: Root of the sweep result tree.
        threads: CPU thread count.
        max_ram_gb: Metabuli RAM limit.

    Returns:
        Completed setting directory.
    """
    if point.reuse_existing:
        raise ValueError(f"Cannot rerun a reused operating point: {point.setting_id}")
    require_file(input_fastq, "Input FASTQ")
    require_directory(database, f"{point.method} database")
    executable = shutil.which(point.method)
    if executable is None:
        raise FileNotFoundError(f"Classifier executable is unavailable: {point.method}")

    method_root = output_root / point.method
    method_root.mkdir(parents=True, exist_ok=True)
    final_directory = method_root / point.setting_id
    complete_path = final_directory / "complete.json"
    if complete_path.is_file() and complete_path.stat().st_size > 0:
        LOGGER.info("Operating point is already complete: %s", final_directory)
        return final_directory
    if final_directory.exists():
        raise FileExistsError(
            "Incomplete result directory already exists; inspect before rerunning: "
            f"{final_directory}"
        )

    staging_directory = Path(
        tempfile.mkdtemp(prefix=f".{point.setting_id}.", dir=method_root)
    )
    started_at = utc_now()
    command = build_classifier_command(
        point,
        input_fastq=input_fastq,
        database=database,
        staging_directory=staging_directory,
        threads=threads,
        max_ram_gb=max_ram_gb,
    )
    LOGGER.info("Running %s setting %s", point.method, point.setting_id)
    LOGGER.debug("Command: %s", command)
    try:
        with (staging_directory / "classifier.stdout.log").open(
            "w", encoding="utf-8"
        ) as stdout_handle, (staging_directory / "classifier.stderr.log").open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            subprocess.run(
                command,
                check=True,
                stdout=stdout_handle,
                stderr=stderr_handle,
            )
        if point.method == "kraken2":
            require_file(staging_directory / "report.tsv", "Kraken2 report")
            require_file(
                staging_directory / "classifications.tsv",
                "Kraken2 classifications",
            )
            classification_path = gzip_file(staging_directory / "classifications.tsv")
        else:
            classification_path = normalise_metabuli_outputs(
                staging_directory,
                point.setting_id,
            )
        report_path = staging_directory / "report.tsv"
        require_file(report_path, "Classifier report")
        require_file(classification_path, "Compressed classifications")
        completed_at = utc_now()
        metadata = {
            "command": command,
            "completed_at_utc": completed_at,
            "database": str(database.resolve()),
            "input_fastq": str(input_fastq.resolve()),
            "operating_point": asdict(point),
            "started_at_utc": started_at,
            "status": "success",
            "threads": threads,
            "tool_version": tool_version(point.method),
        }
        write_json(staging_directory / "metadata.json", metadata)
        completion = {
            "completed_at_utc": completed_at,
            "method": point.method,
            "outputs": {
                "classifications": {
                    "path": "classifications.tsv.gz",
                    "size_bytes": classification_path.stat().st_size,
                },
                "metadata": {
                    "path": "metadata.json",
                    "size_bytes": (staging_directory / "metadata.json").stat().st_size,
                },
                "report": {
                    "path": "report.tsv",
                    "sha256": sha256_file(report_path),
                    "size_bytes": report_path.stat().st_size,
                },
            },
            "setting_id": point.setting_id,
            "status": "success",
        }
        write_json(staging_directory / "complete.json", completion)
        staging_directory.replace(final_directory)
    except Exception as error:
        failure = {
            "error": str(error),
            "failed_at_utc": utc_now(),
            "method": point.method,
            "setting_id": point.setting_id,
            "status": "failure",
        }
        write_json(staging_directory / "failure.json", failure)
        failed_directory = method_root / f"{point.setting_id}.failed.{os.getpid()}"
        staging_directory.replace(failed_directory)
        LOGGER.error("Classifier point failed; diagnostics retained in %s", failed_directory)
        raise
    LOGGER.info("Completed operating point: %s", final_directory)
    return final_directory


def parse_classifier_report(path: Path) -> list[ReportRow]:
    """Parse Kraken2- or Metabuli-style six-column report output.

    Args:
        path: Report TSV.

    Returns:
        Parsed report rows.

    Raises:
        ValueError: If a non-comment data row has fewer than six fields.
    """
    require_file(path, "Classifier report")
    rows: list[ReportRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", maxsplit=5)
            if len(fields) != 6:
                raise ValueError(f"Malformed classifier report row {line_number}: {path}")
            rows.append(
                ReportRow(
                    percentage=float(fields[0].strip()),
                    clade_reads=int(fields[1].strip()),
                    direct_reads=int(fields[2].strip()),
                    rank=fields[3].strip(),
                    taxid=int(fields[4].strip()),
                    name=fields[5].strip(),
                )
            )
    if not rows:
        raise ValueError(f"Classifier report contains no rows: {path}")
    return rows


def load_truth_species(path: Path) -> dict[int, str]:
    """Load expected taxids and display names from the truth manifest.

    Args:
        path: Detailed ATCC truth manifest.

    Returns:
        Mapping from expected NCBI taxid to canonical display name.
    """
    require_file(path, "Truth-species manifest")
    truth: dict[int, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_name", "accepted_species_names", "expected", "ncbi_taxid"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Truth manifest lacks columns: {sorted(missing)}")
        for row in reader:
            if row["expected"].strip().lower() != "true":
                continue
            taxid = int(row["ncbi_taxid"])
            accepted = row["accepted_species_names"].strip()
            display_name = (
                accepted.split(";")[0].strip()
                if accepted
                else row["species_name"].strip()
            )
            if taxid in truth:
                raise ValueError(f"Duplicate expected taxid in truth manifest: {taxid}")
            truth[taxid] = display_name
    if not truth:
        raise ValueError("Truth manifest contains no expected species")
    return truth


def report_path_for_point(
    point: OperatingPoint,
    *,
    sweep_root: Path,
    existing_kraken2_report: Path,
    existing_metabuli_report: Path,
) -> Path:
    """Resolve the report for one completed or reused operating point.

    Args:
        point: Operating point.
        sweep_root: Root holding newly run settings.
        existing_kraken2_report: Existing confidence-zero Kraken2 report.
        existing_metabuli_report: Existing score-0.008 Metabuli report.

    Returns:
        Report path.
    """
    if point.reuse_existing:
        return existing_kraken2_report if point.method == "kraken2" else existing_metabuli_report
    return sweep_root / point.method / point.setting_id / "report.tsv"


def summarise_report(
    rows: Sequence[ReportRow],
    *,
    expected_taxids: set[int],
    focus_taxid: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Summarise expected recovery and additional species by read support.

    Args:
        rows: Parsed classifier report.
        expected_taxids: Expected species taxids.
        focus_taxid: Additional taxid to report separately.

    Returns:
        Summary metrics and species-level direct-count rows.
    """
    root_rows = [row for row in rows if row.taxid == 1]
    unclassified_rows = [row for row in rows if row.taxid == 0]
    if len(root_rows) != 1 or len(unclassified_rows) != 1:
        raise ValueError("Report must contain exactly one root and one unclassified row")
    classified_reads = root_rows[0].clade_reads
    unclassified_reads = unclassified_rows[0].clade_reads
    total_reads = classified_reads + unclassified_reads
    species_rows = [
        row for row in rows if row.rank.lower() in SPECIES_RANKS and row.direct_reads > 0
    ]
    summary: dict[str, object] = {
        "additional_species_ge_1": 0,
        "classified_percent": round(100.0 * classified_reads / total_reads, 6),
        "classified_reads": classified_reads,
        "expected_species_total": len(expected_taxids),
        "focus_taxid_direct_reads": sum(
            row.direct_reads for row in species_rows if row.taxid == focus_taxid
        ),
        "species_labels_ge_1": len(species_rows),
        "total_reads": total_reads,
        "unclassified_reads": unclassified_reads,
    }
    species_output: list[dict[str, object]] = []
    for row in sorted(species_rows, key=lambda item: (-item.direct_reads, item.taxid)):
        species_output.append(
            {
                "direct_reads": row.direct_reads,
                "expected": row.taxid in expected_taxids,
                "species_name": row.name,
                "taxid": row.taxid,
            }
        )
    for threshold in READ_SUPPORT_THRESHOLDS:
        supported = [row for row in species_rows if row.direct_reads >= threshold]
        expected = {row.taxid for row in supported if row.taxid in expected_taxids}
        additional = {row.taxid for row in supported if row.taxid not in expected_taxids}
        summary[f"expected_species_ge_{threshold}"] = len(expected)
        summary[f"additional_species_ge_{threshold}"] = len(additional)
    return summary, species_output


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    """Write dictionaries as a tab-separated table.

    Args:
        path: Destination TSV.
        fieldnames: Ordered output columns.
        rows: Output row dictionaries.
    """
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def summarise_sweep(
    *,
    manifest: Path,
    truth_manifest: Path,
    sweep_root: Path,
    existing_kraken2_report: Path,
    existing_metabuli_report: Path,
    output_directory: Path,
    focus_taxid: int,
) -> Path:
    """Create final cross-method operating-point summary tables.

    Args:
        manifest: Predeclared operating-point manifest.
        truth_manifest: Expected-species manifest.
        sweep_root: Root containing new classifier outputs.
        existing_kraken2_report: Reused confidence-zero report.
        existing_metabuli_report: Reused score-0.008 report.
        output_directory: Final summary directory.
        focus_taxid: Taxid reported separately.

    Returns:
        Completed summary directory.
    """
    points = load_operating_points(manifest)
    truth = load_truth_species(truth_manifest)
    expected_taxids = set(truth)
    if len(expected_taxids) != 20:
        raise ValueError(f"Expected exactly 20 truth species; found {len(expected_taxids)}")
    if output_directory.exists():
        raise FileExistsError(f"Summary output already exists: {output_directory}")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    staging_directory = Path(
        tempfile.mkdtemp(prefix=f".{output_directory.name}.", dir=output_directory.parent)
    )
    summary_rows: list[dict[str, object]] = []
    species_rows: list[dict[str, object]] = []
    try:
        for point in points:
            report_path = report_path_for_point(
                point,
                sweep_root=sweep_root,
                existing_kraken2_report=existing_kraken2_report,
                existing_metabuli_report=existing_metabuli_report,
            )
            rows = parse_classifier_report(report_path)
            metrics, point_species = summarise_report(
                rows,
                expected_taxids=expected_taxids,
                focus_taxid=focus_taxid,
            )
            summary_rows.append(
                {
                    "method": point.method,
                    "setting_id": point.setting_id,
                    "setting_context": point.setting_context,
                    "confidence": "" if point.confidence is None else point.confidence,
                    "min_score": "" if point.min_score is None else point.min_score,
                    "min_sp_score": "" if point.min_sp_score is None else point.min_sp_score,
                    "precise": "" if point.precise is None else point.precise,
                    "reused_existing_run": point.reuse_existing,
                    **metrics,
                    "report_path": str(report_path),
                }
            )
            for species in point_species:
                species_rows.append(
                    {
                        "method": point.method,
                        "setting_id": point.setting_id,
                        **species,
                    }
                )
        summary_fields = [
            "method",
            "setting_id",
            "setting_context",
            "confidence",
            "min_score",
            "min_sp_score",
            "precise",
            "reused_existing_run",
            "total_reads",
            "classified_reads",
            "unclassified_reads",
            "classified_percent",
            "expected_species_total",
            "species_labels_ge_1",
            "expected_species_ge_1",
            "additional_species_ge_1",
            "expected_species_ge_2",
            "additional_species_ge_2",
            "expected_species_ge_10",
            "additional_species_ge_10",
            "expected_species_ge_100",
            "additional_species_ge_100",
            "focus_taxid_direct_reads",
            "report_path",
        ]
        write_tsv(
            staging_directory / "classifier_operating_point_summary.tsv",
            summary_fields,
            summary_rows,
        )
        write_tsv(
            staging_directory / "classifier_species_direct_read_counts.tsv",
            ["method", "setting_id", "taxid", "species_name", "expected", "direct_reads"],
            species_rows,
        )
        shutil.copy2(manifest, staging_directory / "operating_points.tsv")
        shutil.copy2(truth_manifest, staging_directory / "truth_species.tsv")
        readme = (
            "ATCC MSA-1003 PacBio HiFi classifier operating-point sweep\n\n"
            "Kraken2 reports are compared across confidence thresholds. Metabuli reports "
            "are compared across minimum-score and species-rank minimum-score settings. "
            "Counts use direct positive assignments at species rank. Species not present "
            "in the locked 20-species truth manifest are labelled additional. These are "
            "known-composition benchmark labels, not estimates of clinical specificity.\n"
        )
        (staging_directory / "README.txt").write_text(readme, encoding="utf-8")
        complete = {
            "completed_at_utc": utc_now(),
            "focus_taxid": focus_taxid,
            "method_count": len({point.method for point in points}),
            "operating_point_count": len(points),
            "status": "complete",
            "truth_species_count": len(expected_taxids),
        }
        write_json(staging_directory / "complete.json", complete)
        staging_directory.replace(output_directory)
    except Exception:
        LOGGER.exception("Sweep summary failed; staging directory retained: %s", staging_directory)
        raise
    LOGGER.info("Wrote operating-point summary: %s", output_directory)
    return output_directory


def build_parser() -> argparse.ArgumentParser:
    """Build the named-argument command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--action",
        required=True,
        choices=("validate-manifest", "task-count", "run", "summarise"),
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--method", choices=SUPPORTED_METHODS)
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--input-fastq", type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--max-ram-gb", type=int, default=120)
    parser.add_argument("--truth-manifest", type=Path)
    parser.add_argument("--sweep-root", type=Path)
    parser.add_argument("--existing-kraken2-report", type=Path)
    parser.add_argument("--existing-metabuli-report", type=Path)
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--focus-taxid", type=int, default=99158)
    parser.add_argument("--verbose", action="store_true")
    return parser


def required_argument(value: object, name: str) -> object:
    """Require an action-specific parsed argument.

    Args:
        value: Parsed argument value.
        name: CLI option name.

    Returns:
        The validated value.

    Raises:
        ValueError: If the option was omitted.
    """
    if value is None:
        raise ValueError(f"{name} is required for this action")
    return value


def main(arguments: Sequence[str] | None = None) -> int:
    """Run the requested validation, classifier task or summary action.

    Args:
        arguments: Optional command-line arguments for tests.

    Returns:
        Process exit code.
    """
    parser = build_parser()
    args = parser.parse_args(arguments)
    configure_logging(args.verbose)
    try:
        points = load_operating_points(args.manifest)
        if args.action == "validate-manifest":
            LOGGER.info("Validated %d operating points", len(points))
            return 0
        if args.action == "task-count":
            method = str(required_argument(args.method, "--method"))
            print(len(pending_points(points, method)))
            return 0
        if args.action == "run":
            method = str(required_argument(args.method, "--method"))
            task_index = int(required_argument(args.task_index, "--task-index"))
            selected = pending_points(points, method)
            if not 0 <= task_index < len(selected):
                raise IndexError(
                    f"Task index {task_index} is outside 0..{len(selected) - 1} for {method}"
                )
            run_operating_point(
                selected[task_index],
                input_fastq=Path(required_argument(args.input_fastq, "--input-fastq")),
                database=Path(required_argument(args.database, "--database")),
                output_root=Path(required_argument(args.output_root, "--output-root")),
                threads=args.threads,
                max_ram_gb=args.max_ram_gb,
            )
            return 0
        summarise_sweep(
            manifest=args.manifest,
            truth_manifest=Path(
                required_argument(args.truth_manifest, "--truth-manifest")
            ),
            sweep_root=Path(required_argument(args.sweep_root, "--sweep-root")),
            existing_kraken2_report=Path(
                required_argument(
                    args.existing_kraken2_report,
                    "--existing-kraken2-report",
                )
            ),
            existing_metabuli_report=Path(
                required_argument(
                    args.existing_metabuli_report,
                    "--existing-metabuli-report",
                )
            ),
            output_directory=Path(
                required_argument(args.output_directory, "--output-directory")
            ),
            focus_taxid=args.focus_taxid,
        )
        return 0
    except (OSError, RuntimeError, ValueError, IndexError, subprocess.SubprocessError) as error:
        LOGGER.error("%s", error)
        return 2


if __name__ == "__main__":
    sys.exit(main())
