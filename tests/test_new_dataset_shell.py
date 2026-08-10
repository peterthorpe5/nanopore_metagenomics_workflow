"""Integration tests for the defensive new-dataset shell helper."""

from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY_ROOT / "scripts" / "start_new_dataset.sh"


class TestNewDatasetShell(unittest.TestCase):
    """Verify template initialisation without invoking external classifiers."""

    def test_help_describes_named_actions(self) -> None:
        """The shell interface should be self-documenting and non-positional."""
        completed = subprocess.run(
            ["bash", str(SCRIPT), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0)
        self.assertIn("--action initialise|validate|dry-run|run", completed.stdout)

    def test_initialise_creates_run_templates_and_refuses_overwrite(self) -> None:
        """Initialisation creates all templates once and preserves user edits."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = Path(temporary) / "new_dataset.yaml"
            command = [
                "bash",
                str(SCRIPT),
                "--action",
                "initialise",
                "--config",
                str(config_path),
            ]
            first = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            sample_path = config_path.with_suffix(".samples.tsv")
            pcr_path = config_path.with_suffix(".pcr_truth.tsv")
            self.assertEqual(first.returncode, 0)
            self.assertTrue(config_path.is_file())
            self.assertTrue(sample_path.is_file())
            self.assertTrue(pcr_path.is_file())
            config_path.write_text("user edit\n", encoding="utf-8")
            second = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(second.returncode, 3)
            self.assertEqual(config_path.read_text(encoding="utf-8"), "user edit\n")

    def test_failed_workflow_attempts_report_but_preserves_failure_exit_code(self) -> None:
        """A scheduler-style failure should trigger salvage without being hidden."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary_directory = root / "bin"
            binary_directory.mkdir()
            marker = root / "report.called"
            fake = binary_directory / "nanopore-realdata"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                'if [[ " $* " == *" --action report "* ]]; then\n'
                '  : > "${REPORT_MARKER}"\n'
                "  exit 0\n"
                "fi\n"
                "exit 9\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            config = root / "config.yaml"
            config.write_text("synthetic\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    "bash",
                    str(REPOSITORY_ROOT / "scripts" / "run_workflow.sh"),
                    "--config",
                    str(config),
                ],
                check=False,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{binary_directory}{os.pathsep}{os.environ.get('PATH', '')}",
                    "REPORT_MARKER": str(marker),
                },
            )
            report_was_attempted = marker.is_file()
        self.assertEqual(completed.returncode, 9)
        self.assertTrue(report_was_attempted)
        self.assertIn("attempting a truthful partial report", completed.stderr)


if __name__ == "__main__":
    unittest.main()
