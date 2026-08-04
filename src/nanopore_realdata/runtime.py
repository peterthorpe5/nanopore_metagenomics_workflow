"""Runtime, staging and provenance helpers for the real-data workflow."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Iterator, Mapping, Sequence


LOGGER = logging.getLogger(__name__)
NETWORK_FILESYSTEMS = {"cifs", "gpfs", "lustre", "nfs", "nfs4", "smb3"}


def utc_now() -> str:
    """Return the current UTC time in an ISO-8601 representation."""
    return datetime.now(timezone.utc).isoformat()


def write_json_atomic(*, path: Path, payload: Mapping[str, Any]) -> None:
    """Write JSON atomically in the destination directory.

    Args:
        path: Final JSON path.
        payload: JSON-compatible mapping.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def sha256_file(*, path: Path, block_size: int = 1024 * 1024) -> str:
    """Calculate a file SHA-256 digest without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def metadata_fingerprint(*, paths: Sequence[Path], checksum_files: bool) -> str:
    """Fingerprint files or shallow database directories.

    Database directories are deliberately inspected only one level deep. This
    avoids the damaging recursive database audits observed on shared storage.

    Args:
        paths: Input files or directories.
        checksum_files: Hash regular files instead of using size and mtime.

    Returns:
        SHA-256 digest of the normalised resource inventory.
    """
    rows: list[dict[str, Any]] = []
    for path in sorted((item.resolve() for item in paths), key=str):
        if not path.exists():
            raise FileNotFoundError(f"Cannot fingerprint missing path: {path}")
        stat = path.stat()
        row: dict[str, Any] = {
            "path": str(path),
            "type": "directory" if path.is_dir() else "file",
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
        if path.is_file() and checksum_files:
            row["sha256"] = sha256_file(path=path)
        rows.append(row)
        if path.is_dir():
            for child in sorted(path.iterdir(), key=lambda item: item.name):
                child_stat = child.stat()
                rows.append(
                    {
                        "path": str(child.relative_to(path)),
                        "type": "directory" if child.is_dir() else "file",
                        "size": child_stat.st_size,
                        "mtime_ns": child_stat.st_mtime_ns,
                    }
                )
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def task_signature(
    *,
    task: Mapping[str, Any],
    inputs: Sequence[Path],
    checksum_files: bool,
) -> str:
    """Bind a task to its settings and input fingerprints."""
    payload = {
        "task": task,
        "inputs": metadata_fingerprint(paths=inputs, checksum_files=checksum_files),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def completion_is_valid(
    *,
    completion_path: Path,
    signature: str,
    outputs: Sequence[Path],
) -> bool:
    """Validate a completion token and every declared output."""
    if not completion_path.is_file():
        return False
    try:
        payload = json.loads(completion_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return False
    if payload.get("status") != "success" or payload.get("signature") != signature:
        return False
    return all(path.is_file() and path.stat().st_size > 0 for path in outputs)


def write_completion(
    *,
    completion_path: Path,
    signature: str,
    outputs: Sequence[Path],
    extra: Mapping[str, Any] | None = None,
) -> None:
    """Write a success token after validating declared outputs."""
    missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Cannot mark task complete; outputs are missing or empty: {missing}")
    payload: dict[str, Any] = {
        "status": "success",
        "completed_at_utc": utc_now(),
        "signature": signature,
        "outputs": [
            {"path": str(path), "size_bytes": path.stat().st_size}
            for path in outputs
        ],
    }
    if extra:
        payload.update(extra)
    write_json_atomic(path=completion_path, payload=payload)


def run_command(
    *,
    command: Sequence[str],
    log_path: Path,
    stdout_path: Path | None = None,
    timeout_seconds: int | None = None,
) -> None:
    """Run one external command without shell interpretation.

    Args:
        command: Complete argument vector.
        log_path: Persistent or local standard-error log.
        stdout_path: Optional file receiving standard output.
        timeout_seconds: Optional wall-clock timeout.

    Raises:
        RuntimeError: If the command exits unsuccessfully.
    """
    if not command:
        raise ValueError("Command must not be empty")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Running command: %s", " ".join(command))
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"[{utc_now()}] COMMAND: {' '.join(command)}\n")
        output_handle: BinaryIO | None = None
        try:
            if stdout_path is not None:
                # External tools such as pigz write compressed bytes. Opening the
                # descriptor in binary mode prevents a text wrapper from altering
                # or rejecting those bytes.
                output_handle = stdout_path.open("wb")
            try:
                completed = subprocess.run(
                    list(command),
                    check=False,
                    stdout=output_handle if output_handle is not None else log_handle,
                    stderr=log_handle,
                    text=True,
                    timeout=timeout_seconds,
                )
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(
                    f"Command exceeded its {timeout_seconds}-second time limit: {command[0]}"
                ) from error
        finally:
            if output_handle is not None:
                output_handle.close()
        if completed.returncode != 0:
            raise RuntimeError(
                f"Command failed with exit code {completed.returncode}: {command[0]}; "
                f"see {log_path}"
            )


def run_pipeline(
    *,
    commands: Sequence[Sequence[str]],
    log_path: Path,
    stdout_path: Path | None = None,
) -> None:
    """Run a shell-free subprocess pipeline and check every component."""
    if not commands or any(not command for command in commands):
        raise ValueError("Pipeline commands must not be empty")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(
            f"[{utc_now()}] PIPELINE: "
            + " | ".join(" ".join(command) for command in commands)
            + "\n"
        )
        processes: list[subprocess.Popen[str]] = []
        output_handle: BinaryIO | None = None
        previous_stdout = None
        try:
            if stdout_path is not None:
                output_handle = stdout_path.open("wb")
            for index, command in enumerate(commands):
                is_last = index == len(commands) - 1
                process = subprocess.Popen(
                    list(command),
                    stdin=previous_stdout,
                    stdout=(
                        output_handle
                        if is_last and output_handle is not None
                        else log_handle
                        if is_last
                        else subprocess.PIPE
                    ),
                    stderr=log_handle,
                    text=True,
                )
                if previous_stdout is not None:
                    previous_stdout.close()
                previous_stdout = process.stdout if not is_last else None
                processes.append(process)
            return_codes = [process.wait() for process in processes]
        finally:
            if previous_stdout is not None:
                previous_stdout.close()
            if output_handle is not None:
                output_handle.close()
        failures = [
            f"{commands[index][0]}={return_code}"
            for index, return_code in enumerate(return_codes)
            if return_code != 0
        ]
        if failures:
            raise RuntimeError(f"Pipeline failed ({', '.join(failures)}); see {log_path}")


def capture_output(*, command: Sequence[str]) -> str:
    """Run a small command and return stripped standard output."""
    completed = subprocess.run(
        list(command),
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Version or count command failed: {' '.join(command)}")
    return completed.stdout.strip()


def filesystem_type(*, path: Path) -> str | None:
    """Return the filesystem type when ``findmnt`` is available."""
    executable = shutil.which("findmnt")
    if executable is None:
        return None
    completed = subprocess.run(
        [executable, "--noheadings", "--output", "FSTYPE", "--target", str(path)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    if completed.returncode != 0:
        return None
    value = completed.stdout.strip().splitlines()
    return value[0].strip().lower() if value else None


def validate_scratch(*, scratch_root: Path, minimum_gb: int) -> dict[str, Any]:
    """Verify scratch is writable, large enough and not network-hosted."""
    resolved = scratch_root.resolve()
    if not resolved.is_dir():
        raise ValueError(f"Scratch root is not a directory: {resolved}")
    filesystem = filesystem_type(path=resolved)
    if filesystem in NETWORK_FILESYSTEMS:
        raise ValueError(
            f"Scratch root uses network filesystem {filesystem!r}: {resolved}"
        )
    usage = shutil.disk_usage(resolved)
    required = minimum_gb * 1024**3
    if usage.free < required:
        raise ValueError(
            f"Scratch root has {usage.free / 1024**3:.1f} GiB free; "
            f"at least {minimum_gb} GiB is required"
        )
    with tempfile.NamedTemporaryFile(dir=resolved, prefix="kmersutra_write_probe."):
        pass
    return {
        "path": str(resolved),
        "filesystem_type": filesystem or "unknown",
        "free_bytes": usage.free,
        "minimum_required_bytes": required,
    }


@contextmanager
def scratch_workspace(*, scratch_root: Path, label: str) -> Iterator[Path]:
    """Create and always remove a bounded job-local workspace."""
    safe_label = "".join(character if character.isalnum() else "_" for character in label)
    workspace = Path(
        tempfile.mkdtemp(prefix=f"nanopore_realdata.{safe_label}.", dir=scratch_root)
    )
    LOGGER.info("Created node-local workspace: %s", workspace)
    try:
        yield workspace
    finally:
        shutil.rmtree(workspace, ignore_errors=True)
        LOGGER.info("Removed node-local workspace: %s", workspace)


def stage_resource(*, source: Path, destination_root: Path, log_path: Path) -> Path:
    """Stage a file or directory to an empty node-local resource directory."""
    destination_root.mkdir(parents=True, exist_ok=False)
    destination = destination_root / source.name
    if source.is_dir():
        destination.mkdir()
        command = [
            "rsync",
            "--archive",
            "--copy-links",
            f"{source}/",
            f"{destination}/",
        ]
    elif source.is_file():
        command = ["rsync", "--archive", "--copy-links", str(source), str(destination)]
    else:
        raise ValueError(f"Resource is neither a file nor directory: {source}")
    run_command(command=command, log_path=log_path)
    _validate_staged_inventory(source=source, destination=destination)
    return destination


def publish_directory(*, source: Path, destination: Path, log_path: Path) -> None:
    """Rsync a completed result directory and replace its destination atomically."""
    if not source.is_dir():
        raise ValueError(f"Publication source is not a directory: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = destination.parent / f".{destination.name}.partial.{uuid.uuid4().hex}"
    backup = destination.parent / f".{destination.name}.previous.{uuid.uuid4().hex}"
    staging.mkdir()
    try:
        run_command(
            command=["rsync", "--archive", f"{source}/", f"{staging}/"],
            log_path=log_path,
        )
        if not any(staging.iterdir()):
            raise RuntimeError(f"Refusing to publish an empty directory: {source}")
        if destination.exists():
            destination.rename(backup)
        staging.rename(destination)
        shutil.rmtree(backup, ignore_errors=True)
    except Exception:
        if backup.exists() and not destination.exists():
            backup.rename(destination)
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _validate_staged_inventory(*, source: Path, destination: Path) -> None:
    if source.is_file():
        if not destination.is_file() or source.stat().st_size != destination.stat().st_size:
            raise RuntimeError(f"Staged file validation failed: {source}")
        return
    source_rows = _relative_size_inventory(root=source)
    destination_rows = _relative_size_inventory(root=destination)
    if source_rows != destination_rows:
        raise RuntimeError(f"Staged directory inventory differs from source: {source}")


def _relative_size_inventory(*, root: Path) -> list[tuple[str, int]]:
    return sorted(
        (str(path.relative_to(root)), path.stat().st_size)
        for path in root.rglob("*")
        if path.is_file()
    )
