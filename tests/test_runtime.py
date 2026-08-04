"""Tests for staging, subprocess, atomic publication and provenance helpers."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from nanopore_realdata.runtime import (
    completion_is_valid,
    metadata_fingerprint,
    publish_directory,
    run_command,
    run_pipeline,
    scratch_workspace,
    stage_resource,
    task_signature,
    validate_scratch,
    write_completion,
    write_json_atomic,
)


class TestRuntime(unittest.TestCase):
    """Exercise defensive runtime paths without scientific databases."""

    def test_atomic_json_and_completion_validation(self) -> None:
        """Completion requires matching signatures and non-empty outputs."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output.tsv"
            output.write_text("header\n", encoding="utf-8")
            completion = root / "complete.json"
            write_completion(
                completion_path=completion,
                signature="abc",
                outputs=[output],
                extra={"sample_id": "sample_1"},
            )
            self.assertTrue(
                completion_is_valid(
                    completion_path=completion,
                    signature="abc",
                    outputs=[output],
                )
            )
            self.assertFalse(
                completion_is_valid(
                    completion_path=completion,
                    signature="different",
                    outputs=[output],
                )
            )
            output.write_text("", encoding="utf-8")
            self.assertFalse(
                completion_is_valid(
                    completion_path=completion,
                    signature="abc",
                    outputs=[output],
                )
            )

    def test_completion_rejects_empty_output_and_invalid_json(self) -> None:
        """Tokens cannot hide empty outputs or malformed prior state."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "empty"
            output.touch()
            with self.assertRaisesRegex(RuntimeError, "missing or empty"):
                write_completion(
                    completion_path=root / "complete.json",
                    signature="abc",
                    outputs=[output],
                )
            completion = root / "bad.json"
            completion.write_text("not json", encoding="utf-8")
            self.assertFalse(
                completion_is_valid(
                    completion_path=completion,
                    signature="abc",
                    outputs=[output],
                )
            )

    def test_fingerprints_change_with_file_content_and_task(self) -> None:
        """Checksummed files and task settings should bind restart state."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.txt"
            path.write_text("first", encoding="utf-8")
            first = metadata_fingerprint(paths=[path], checksum_files=True)
            signature = task_signature(
                task={"threads": 1},
                inputs=[path],
                checksum_files=True,
            )
            path.write_text("second", encoding="utf-8")
            second = metadata_fingerprint(paths=[path], checksum_files=True)
            changed = task_signature(
                task={"threads": 2},
                inputs=[path],
                checksum_files=True,
            )
        self.assertNotEqual(first, second)
        self.assertNotEqual(signature, changed)

    def test_metadata_fingerprint_is_shallow_for_directories(self) -> None:
        """Database fingerprinting should use direct metadata, not file content."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = root / "db"
            database.mkdir()
            (database / "core.dat").write_text("data", encoding="utf-8")
            nested = database / "library"
            nested.mkdir()
            fingerprint = metadata_fingerprint(
                paths=[database],
                checksum_files=False,
            )
        self.assertEqual(len(fingerprint), 64)

    def test_run_command_writes_stdout_and_reports_failure(self) -> None:
        """Shell-free commands should expose success and failure clearly."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "command.log"
            output = root / "output.txt"
            run_command(
                command=[sys.executable, "-c", "print('hello')"],
                log_path=log,
                stdout_path=output,
            )
            self.assertEqual(output.read_text(encoding="utf-8"), "hello\n")
            with self.assertRaisesRegex(RuntimeError, "exit code 3"):
                run_command(
                    command=[sys.executable, "-c", "raise SystemExit(3)"],
                    log_path=log,
                )

    def test_run_command_timeout_is_a_controlled_failure(self) -> None:
        """A timed-out KmerSutra process should become a normal runtime error."""
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(RuntimeError, "time limit"):
                run_command(
                    command=[sys.executable, "-c", "import time; time.sleep(2)"],
                    log_path=Path(temporary) / "timeout.log",
                    timeout_seconds=1,
                )

    def test_run_pipeline_checks_every_component(self) -> None:
        """A pipeline should publish output and detect upstream failures."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "upper.txt"
            log = root / "pipeline.log"
            run_pipeline(
                commands=[["printf", "hello"], ["tr", "a-z", "A-Z"]],
                log_path=log,
                stdout_path=output,
            )
            self.assertEqual(output.read_text(encoding="utf-8"), "HELLO")
            with self.assertRaisesRegex(RuntimeError, "Pipeline failed"):
                run_pipeline(
                    commands=[["false"], ["cat"]],
                    log_path=log,
                )

    def test_scratch_validation_rejects_network_and_low_capacity(self) -> None:
        """A GPFS scratch path or insufficient capacity should fail safely."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch("nanopore_realdata.runtime.filesystem_type", return_value="gpfs"):
                with self.assertRaisesRegex(ValueError, "network filesystem"):
                    validate_scratch(scratch_root=root, minimum_gb=1)
            fake_usage = shutil._ntuple_diskusage(total=10, used=9, free=1)
            with (
                patch("nanopore_realdata.runtime.filesystem_type", return_value="xfs"),
                patch("nanopore_realdata.runtime.shutil.disk_usage", return_value=fake_usage),
            ):
                with self.assertRaisesRegex(ValueError, "at least"):
                    validate_scratch(scratch_root=root, minimum_gb=1)

    def test_scratch_workspace_is_removed_after_failure(self) -> None:
        """Node-local data must not accumulate when a tool raises."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            observed: Path | None = None
            with self.assertRaisesRegex(RuntimeError, "test failure"):
                with scratch_workspace(scratch_root=root, label="test") as workspace:
                    observed = workspace
                    (workspace / "large.tmp").write_text("temporary", encoding="utf-8")
                    raise RuntimeError("test failure")
            self.assertIsNotNone(observed)
            self.assertFalse(observed.exists())  # type: ignore[union-attr]

    def test_stage_resource_validates_file_and_directory_inventories(self) -> None:
        """Staging should copy and validate both files and database directories."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_file = root / "source.txt"
            source_file.write_text("content", encoding="utf-8")
            database = root / "database"
            database.mkdir()
            (database / "core").write_text("db", encoding="utf-8")

            def fake_rsync(*, command, log_path, stdout_path=None, timeout_seconds=None):
                del log_path, stdout_path, timeout_seconds
                source = Path(command[-2].rstrip("/"))
                destination = Path(command[-1].rstrip("/"))
                if source.is_dir():
                    shutil.copytree(source, destination, dirs_exist_ok=True)
                else:
                    shutil.copy2(source, destination)

            with patch("nanopore_realdata.runtime.run_command", side_effect=fake_rsync):
                staged_file = stage_resource(
                    source=source_file,
                    destination_root=root / "stage_file",
                    log_path=root / "stage.log",
                )
                staged_database = stage_resource(
                    source=database,
                    destination_root=root / "stage_database",
                    log_path=root / "stage.log",
                )
            self.assertEqual(staged_file.read_text(encoding="utf-8"), "content")
            self.assertTrue((staged_database / "core").is_file())

    def test_publish_directory_replaces_only_after_staging(self) -> None:
        """A prior result should survive until a complete replacement is ready."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "result.tsv").write_text("new", encoding="utf-8")
            destination = root / "destination"
            destination.mkdir()
            (destination / "result.tsv").write_text("old", encoding="utf-8")

            def fake_rsync(*, command, log_path, stdout_path=None, timeout_seconds=None):
                del log_path, stdout_path, timeout_seconds
                shutil.copytree(
                    Path(command[-2].rstrip("/")),
                    Path(command[-1].rstrip("/")),
                    dirs_exist_ok=True,
                )

            with patch("nanopore_realdata.runtime.run_command", side_effect=fake_rsync):
                publish_directory(
                    source=source,
                    destination=destination,
                    log_path=root / "publish.log",
                )
            self.assertEqual(
                (destination / "result.tsv").read_text(encoding="utf-8"),
                "new",
            )

    def test_write_json_atomic_produces_parseable_payload(self) -> None:
        """Atomic JSON helper should leave no partial sibling."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "metadata.json"
            write_json_atomic(path=path, payload={"status": "success"})
            payload = json.loads(path.read_text(encoding="utf-8"))
            partials = list(path.parent.glob(".metadata.json.partial.*"))
        self.assertEqual(payload["status"], "success")
        self.assertFalse(partials)


if __name__ == "__main__":
    unittest.main()
