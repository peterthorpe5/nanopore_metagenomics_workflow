"""Tests for detached, failure-aware Slurm submission planning."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from helpers import build_test_project
from nanopore_realdata.config import load_workflow_config
from nanopore_realdata.slurm import (
    JobPlan,
    JobResources,
    _dependency_ids,
    _load_journal,
    _refuse_active_jobs,
    _sbatch_command,
    _submit,
    build_submission_plan,
    planned_commands,
    submit_workflow,
)


class TestSlurmPlan(unittest.TestCase):
    """Keep method arrays independent and final reporting reachable."""

    def test_plan_is_per_sample_and_failure_aware(self) -> None:
        """No classifier should wait for a different classifier to succeed."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(
                root=Path(temporary),
                input_read_state="host_removed",
            )
            workflow = load_workflow_config(config_path=config_path)
            plans = {plan.key: plan for plan in build_submission_plan(workflow=workflow)}
        for method in ("kraken2", "metabuli", "minimap2", "kmersutra"):
            plan = plans[f"classify_{method}"]
            self.assertEqual(plan.array, "0-0%1")
            self.assertEqual(plan.dependency_mode, "afterany")
            self.assertNotIn("classify_", " ".join(plan.dependency_keys))
        self.assertIn(
            "build_minimap_index",
            plans["classify_minimap2"].dependency_keys,
        )
        self.assertNotIn(
            "build_minimap_index",
            plans["classify_kmersutra"].dependency_keys,
        )
        self.assertEqual(plans["classify_kmersutra"].resources.qos, "test_long")
        self.assertEqual(plans["classify_kraken2"].resources.qos, "")
        self.assertEqual(plans["aggregate"].dependency_mode, "afterany")

    def test_sbatch_command_passes_array_index_to_named_cli(self) -> None:
        """Each array task must resolve exactly one manifest sample."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config_path = build_test_project(
                root=root,
                input_read_state="host_removed",
            )
            workflow = load_workflow_config(config_path=config_path)
            plan = next(
                item
                for item in build_submission_plan(workflow=workflow)
                if item.key == "classify_kmersutra"
            )
            command = _sbatch_command(
                workflow=workflow,
                plan=plan,
                dependency_ids=("101",),
                stage_script=root / "run_stage.sh",
                repository_root=root,
                log_root=root / "logs",
            )
        self.assertIn("--array=0-0%1", command)
        self.assertIn("--qos=test_long", command)
        self.assertIn("--dependency=afterany:101", command)
        self.assertEqual(
            command[-3:],
            ["--method", "kmersutra", "--sample-index-from-slurm"],
        )

    def test_plan_variants_reject_raw_and_route_generated_reference(self) -> None:
        """Backend limits and the controlled-reference stage must be explicit."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(root=Path(temporary))
            workflow = load_workflow_config(config_path=config_path)
            with self.assertRaisesRegex(ValueError, "host_removed"):
                build_submission_plan(workflow=workflow)
            with self.assertRaisesRegex(ValueError, "At least one sample"):
                build_submission_plan(
                    workflow=replace(
                        workflow,
                        input_read_state="host_removed",
                        samples=(),
                    )
                )
            plans = build_submission_plan(
                workflow=replace(
                    workflow,
                    input_read_state="host_removed",
                    minimap_genome_config=Path("genomes.tsv"),
                )
            )
            keys = [plan.key for plan in plans]
            index = next(plan for plan in plans if plan.key == "build_minimap_index")
            disabled = build_submission_plan(
                workflow=replace(
                    workflow,
                    input_read_state="host_removed",
                    kmersutra_enabled=False,
                )
            )
        self.assertIn("build_minimap_reference", keys)
        self.assertEqual(index.dependency_keys, ("build_minimap_reference",))
        self.assertNotIn("classify_kmersutra", {plan.key for plan in disabled})


class TestSlurmSubmission(unittest.TestCase):
    """Exercise journalled submission, resume and low-level failures."""

    def test_submission_journals_every_job_and_refuses_ambiguous_repeat(self) -> None:
        """A complete submission is durable and cannot be repeated silently."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(
                root=Path(temporary),
                input_read_state="host_removed",
            )
            identifiers = (str(number) for number in range(100, 200))
            with patch(
                "nanopore_realdata.slurm._submit",
                side_effect=lambda **_: next(identifiers),
            ) as submitted:
                journal_path = submit_workflow(
                    config_path=config_path,
                    resume_submission=False,
                    new_attempt=False,
                )
            payload = json.loads(journal_path.read_text(encoding="utf-8"))
            plan = planned_commands(config_path=config_path)
            with self.assertRaisesRegex(RuntimeError, "already exists"):
                submit_workflow(
                    config_path=config_path,
                    resume_submission=False,
                    new_attempt=False,
                )
            with self.assertRaisesRegex(RuntimeError, "interrupted submission"):
                submit_workflow(
                    config_path=config_path,
                    resume_submission=True,
                    new_attempt=False,
                )
        self.assertEqual(payload["status"], "submitted")
        self.assertEqual(len(payload["jobs"]), 8)
        self.assertEqual(submitted.call_count, 8)
        self.assertEqual(len(plan), 8)
        self.assertEqual(plan[-1]["key"], "aggregate")

    def test_interrupted_submission_resumes_only_missing_jobs(self) -> None:
        """Resume must preserve accepted IDs and continue from the journal."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(
                root=Path(temporary),
                input_read_state="host_removed",
            )
            with patch(
                "nanopore_realdata.slurm._submit",
                side_effect=("101", RuntimeError("sbatch unavailable")),
            ):
                with self.assertRaisesRegex(RuntimeError, "sbatch unavailable"):
                    submit_workflow(
                        config_path=config_path,
                        resume_submission=False,
                        new_attempt=False,
                    )
            identifiers = (str(number) for number in range(102, 120))
            with patch(
                "nanopore_realdata.slurm._submit",
                side_effect=lambda **_: next(identifiers),
            ) as resumed:
                journal_path = submit_workflow(
                    config_path=config_path,
                    resume_submission=True,
                    new_attempt=False,
                )
            payload = json.loads(journal_path.read_text(encoding="utf-8"))
        self.assertEqual(payload["jobs"]["preflight"]["job_id"], "101")
        self.assertEqual(payload["status"], "submitted")
        self.assertEqual(resumed.call_count, 7)

    def test_new_attempt_archives_after_active_job_check(self) -> None:
        """A finished attempt may be archived and resubmitted deliberately."""
        with tempfile.TemporaryDirectory() as temporary:
            config_path = build_test_project(
                root=Path(temporary),
                input_read_state="host_removed",
            )
            identifiers = (str(number) for number in range(200, 300))
            with patch(
                "nanopore_realdata.slurm._submit",
                side_effect=lambda **_: next(identifiers),
            ):
                journal_path = submit_workflow(
                    config_path=config_path,
                    resume_submission=False,
                    new_attempt=False,
                )
            with (
                patch("nanopore_realdata.slurm._refuse_active_jobs") as checked,
                patch(
                    "nanopore_realdata.slurm._submit",
                    side_effect=lambda **_: next(identifiers),
                ),
            ):
                submit_workflow(
                    config_path=config_path,
                    resume_submission=False,
                    new_attempt=True,
                )
            history = list((journal_path.parent / "submission_history").glob("*.json"))
        checked.assert_called_once()
        self.assertEqual(len(history), 1)

    def test_journal_and_dependency_validation_reject_corruption(self) -> None:
        """Malformed state must stop before an unsafe sbatch call."""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "journal.json"
            path.write_text("not json", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Invalid"):
                _load_journal(path=path)
            path.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Malformed"):
                _load_journal(path=path)
        plan = JobPlan(
            key="dependent",
            name="dependent",
            action="aggregate",
            resources=JobResources(1, 100, 10),
            dependency_keys=("missing",),
        )
        with self.assertRaisesRegex(RuntimeError, "not been submitted"):
            _dependency_ids(plan=plan, jobs={})
        with self.assertRaisesRegex(RuntimeError, "invalid job ID"):
            _dependency_ids(plan=plan, jobs={"missing": {"job_id": "bad"}})

    def test_sbatch_response_and_active_job_checks(self) -> None:
        """Scheduler responses and retry safety checks are strict."""
        with patch(
            "nanopore_realdata.slurm.subprocess.run",
            return_value=SimpleNamespace(returncode=0, stdout="12345;cluster\n"),
        ):
            self.assertEqual(_submit(command=["sbatch"]), "12345")
        with patch(
            "nanopore_realdata.slurm.subprocess.run",
            return_value=SimpleNamespace(returncode=1, stdout="rejected"),
        ):
            with self.assertRaisesRegex(RuntimeError, "rejected"):
                _submit(command=["sbatch"])
        with patch(
            "nanopore_realdata.slurm.subprocess.run",
            return_value=SimpleNamespace(returncode=0, stdout="not-a-job"),
        ):
            with self.assertRaisesRegex(RuntimeError, "Unexpected"):
                _submit(command=["sbatch"])

        journal = {"jobs": {"preflight": {"job_id": "12345"}}}
        with patch.dict(os.environ, {"USER": ""}, clear=False):
            with self.assertRaisesRegex(RuntimeError, "USER"):
                _refuse_active_jobs(journal=journal)
        with (
            patch.dict(os.environ, {"USER": "tester"}, clear=False),
            patch(
                "nanopore_realdata.slurm.subprocess.run",
                return_value=SimpleNamespace(returncode=0, stdout="12345\n"),
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "still queued"):
                _refuse_active_jobs(journal=journal)
        with (
            patch.dict(os.environ, {"USER": "tester"}, clear=False),
            patch(
                "nanopore_realdata.slurm.subprocess.run",
                return_value=SimpleNamespace(returncode=0, stdout=""),
            ),
        ):
            _refuse_active_jobs(journal=journal)

        with self.assertRaisesRegex(RuntimeError, "jobs mapping"):
            _refuse_active_jobs(journal={"jobs": []})
        _refuse_active_jobs(journal={"jobs": {}})

    def test_mutually_exclusive_submission_flags_fail_before_loading(self) -> None:
        """Conflicting recovery intent is rejected deterministically."""
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            submit_workflow(
                config_path=Path("unused.yaml"),
                resume_submission=True,
                new_attempt=True,
            )


if __name__ == "__main__":
    unittest.main()
