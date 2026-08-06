"""Safe command construction for real Nanopore classification tools."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def minimap2_index_command(*, reference: Path, output_index: Path) -> list[str]:
    """Build the minimap2 host-index command.

    Args:
        reference: Host reference FASTA.
        output_index: Destination MMI path.

    Returns:
        Argument vector suitable for ``subprocess`` without a shell.
    """
    return ["minimap2", "-d", str(output_index), str(reference)]


def minimap2_host_command(
    *,
    host_index: Path,
    fastq_paths: Sequence[Path],
    threads: int,
) -> list[str]:
    """Build a primary-only ONT host-alignment command.

    Args:
        host_index: Host reference or MMI index.
        fastq_paths: One or more FASTQ parts for a logical sample.
        threads: Worker threads.

    Returns:
        Minimap2 argument vector.
    """
    _require_positive(threads=threads)
    if not fastq_paths:
        raise ValueError("At least one FASTQ is required for host depletion")
    return [
        "minimap2",
        "-ax",
        "map-ont",
        "--secondary=no",
        "-t",
        str(threads),
        str(host_index),
        *(str(path) for path in fastq_paths),
    ]


def minimap2_classification_command(
    *,
    reference_index: Path,
    input_fastq: Path,
    threads: int,
) -> list[str]:
    """Build the controlled-reference ONT classification command.

    Args:
        reference_index: Minimap2 MMI index built from the classification FASTA.
        input_fastq: Host-removed analysis reads.
        threads: Worker threads.

    Returns:
        Command producing PAF with CIGAR tags on standard output.
    """
    _require_positive(threads=threads)
    return [
        "minimap2",
        "-x",
        "map-ont",
        "--secondary=yes",
        "-c",
        "-t",
        str(threads),
        str(reference_index),
        str(input_fastq),
    ]


def samtools_bam_command(*, output_bam: Path, threads: int) -> list[str]:
    """Build the SAM-to-BAM command used after host alignment."""
    _require_positive(threads=threads)
    return ["samtools", "view", "-b", "-@", str(threads), "-o", str(output_bam), "-"]


def samtools_count_command(
    *,
    bam_path: Path,
    threads: int,
    unmapped_only: bool,
) -> list[str]:
    """Build a primary-record count command for host-depletion metrics."""
    _require_positive(threads=threads)
    command = ["samtools", "view", "-c", "-@", str(threads), "-F", "2304"]
    if unmapped_only:
        command.extend(["-f", "4"])
    command.append(str(bam_path))
    return command


def samtools_non_host_fastq_command(*, bam_path: Path, threads: int) -> list[str]:
    """Build a command that emits primary unmapped reads as FASTQ."""
    _require_positive(threads=threads)
    return [
        "samtools",
        "fastq",
        "-@",
        str(threads),
        "-f",
        "4",
        "-F",
        "2304",
        str(bam_path),
    ]


def pigz_command(*, threads: int) -> list[str]:
    """Build a deterministic gzip-compression command."""
    _require_positive(threads=threads)
    return ["pigz", "-c", "-p", str(threads)]


def kraken2_command(
    *,
    input_fastq: Path,
    database: Path,
    classifications: Path,
    report: Path,
    threads: int,
    confidence: float,
) -> list[str]:
    """Build the Kraken2 read-classification command."""
    _require_positive(threads=threads)
    if not 0.0 <= confidence <= 1.0:
        raise ValueError("Kraken2 confidence must be between zero and one")
    command = [
        "kraken2",
        "--db",
        str(database),
        "--threads",
        str(threads),
        "--confidence",
        str(confidence),
        "--report",
        str(report),
        "--output",
        str(classifications),
    ]
    if str(input_fastq).lower().endswith(".gz"):
        command.append("--gzip-compressed")
    command.append(str(input_fastq))
    return command


def metabuli_command(
    *,
    input_fastq: Path,
    database: Path,
    output_directory: Path,
    output_prefix: str,
    threads: int,
    maximum_ram_gb: int,
    minimum_score: float,
) -> list[str]:
    """Build the Metabuli ONT read-classification command."""
    _require_positive(threads=threads)
    if maximum_ram_gb <= 0:
        raise ValueError("Metabuli maximum RAM must be positive")
    if not 0.0 <= minimum_score <= 1.0:
        raise ValueError("Metabuli minimum score must be between zero and one")
    if not output_prefix or "/" in output_prefix:
        raise ValueError("Metabuli output prefix must be a plain filename prefix")
    return [
        "metabuli",
        "classify",
        "--seq-mode",
        "3",
        "--threads",
        str(threads),
        "--max-ram",
        str(maximum_ram_gb),
        "--min-score",
        str(minimum_score),
        str(input_fastq),
        str(database),
        str(output_directory),
        output_prefix,
    ]


def kmersutra_command(
    *,
    input_fastq: Path,
    panel: Path,
    sample_id: str,
    output_directory: Path,
    threads: int,
    screen_preset: str,
    call_preset: str,
    same_genus_fraction: float,
    write_parquet: bool,
) -> list[str]:
    """Build the KmerSutra exact-evidence screening command."""
    _require_positive(threads=threads)
    if not 0.0 <= same_genus_fraction <= 1.0:
        raise ValueError("Same-genus reportability fraction must be between zero and one")
    command = [
        "kmersutra-screen",
        "--input",
        str(input_fastq),
        "--input_format",
        "fastq",
        "--panel",
        str(panel),
        "--sample_id",
        sample_id,
        "--out_dir",
        str(output_directory),
        "--screen_mode",
        "flat",
        "--screen_preset",
        screen_preset,
        "--call_preset",
        call_preset,
        "--same_genus_reportable_min_fraction",
        str(same_genus_fraction),
        "--consolidate_species_calls",
        "--threads",
        str(threads),
        "--chunk_size",
        "10000",
        "--decompressor",
        "auto",
        "--use_panel_cache",
        "--no_read_level_hits",
        "--profile",
        "--verbose",
    ]
    if write_parquet:
        command.append("--write_parquet_outputs")
    return command


def required_executables(*, action: str) -> tuple[str, ...]:
    """Return external executables required for a workflow action."""
    mapping = {
        "build-host-index": ("minimap2",),
        "accept-host-removed": ("rsync",),
        "host-deplete": ("minimap2", "samtools", "pigz", "rsync"),
        "classify-kraken2": ("kraken2", "pigz", "rsync"),
        "classify-metabuli": ("metabuli", "pigz", "rsync"),
        "classify-minimap2": ("minimap2", "pigz", "rsync"),
        "classify-kmersutra": ("kmersutra-screen", "rsync"),
    }
    try:
        return mapping[action]
    except KeyError as error:
        raise ValueError(f"Unknown workflow action: {action}") from error


def _require_positive(*, threads: int) -> None:
    if threads <= 0:
        raise ValueError("Thread count must be positive")
