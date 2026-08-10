"""Offline, failure-aware HTML reporting for classifier benchmark results."""

from __future__ import annotations

import html
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from nanopore_realdata.config import WorkflowConfig
from nanopore_realdata.runtime import utc_now, write_json_atomic


LOGGER = logging.getLogger(__name__)
METHODS = ("kraken2", "metabuli", "minimap2", "kmersutra")
METHOD_LABELS = {
    "kraken2": "Kraken2",
    "metabuli": "Metabuli",
    "minimap2": "minimap2",
    "kmersutra": "KmerSutra",
}
SUCCESS_STATUSES = {"success"}
INCOMPLETE_STATUSES = {
    "disabled",
    "failed",
    "invalid",
    "missing",
    "partial",
    "scheduler_failed",
    "skipped",
    "timeout",
    "unavailable",
}


def build_normalised_evidence(
    *,
    classifier_rows: Sequence[Mapping[str, Any]],
    minimap_rows: Sequence[Mapping[str, Any]],
    kmersutra_rows: Sequence[Mapping[str, Any]],
    focus_taxa: Sequence[str],
) -> list[dict[str, Any]]:
    """Convert method-specific records into a conservative comparison table.

    Args:
        classifier_rows: Harmonised Kraken2 and Metabuli report rows.
        minimap_rows: Controlled-reference minimap2 report rows.
        kmersutra_rows: KmerSutra species-call rows.
        focus_taxa: Case-insensitive taxon-name fragments to flag.

    Returns:
        Normalised evidence rows. No method-specific value is promoted to a
        cross-method consensus call.
    """
    evidence: list[dict[str, Any]] = []
    for row in classifier_rows:
        name = _first_text(row=row, keys=("taxon_name",))
        rank = _first_text(row=row, keys=("rank_code",))
        direct_reads = _non_negative_number(row.get("direct_reads"))
        clade_reads = _non_negative_number(row.get("clade_reads"))
        evidence.append(
            _evidence_row(
                row=row,
                taxon_name=name,
                tax_id=_first_text(row=row, keys=("tax_id",)),
                rank=rank,
                evidence_count=direct_reads,
                supporting_count=clade_reads,
                fraction=_non_negative_number(row.get("fraction_percent")),
                metric="direct reads",
                detected=direct_reads > 0,
                comparable=rank.upper().startswith("S"),
                focus_taxa=focus_taxa,
            )
        )
    for row in minimap_rows:
        name = _first_text(row=row, keys=("taxon_name", "reference_name"))
        best_reads = _non_negative_number(row.get("best_read_count"))
        evidence.append(
            _evidence_row(
                row=row,
                taxon_name=name,
                tax_id=_first_text(row=row, keys=("tax_id",)),
                rank="controlled_reference",
                evidence_count=best_reads,
                supporting_count=_non_negative_number(row.get("alignment_count")),
                fraction=0.0,
                metric="best-aligned reads",
                detected=best_reads > 0,
                comparable=True,
                focus_taxa=focus_taxa,
            )
        )
    for row in kmersutra_rows:
        name = _first_text(
            row=row,
            keys=(
                "species_name",
                "taxon_name",
                "organism_name",
                "species",
                "name",
            ),
        )
        count = _first_number(
            row=row,
            keys=(
                "n_unique_kmers",
                "n_exact_hits",
                "n_hits",
                "supporting_read_count",
                "total_supporting_reads",
                "read_count",
                "supporting_kmer_count",
                "observed_kmer_count",
                "matched_kmers",
            ),
        )
        call = _first_text(
            row=row,
            keys=(
                "detection_call",
                "final_call",
                "call",
                "decision",
                "status",
            ),
        )
        detected = _call_is_positive(value=call) if call else count > 0
        evidence.append(
            _evidence_row(
                row=row,
                taxon_name=name,
                tax_id=_first_text(row=row, keys=("tax_id", "taxonomy_id")),
                rank="species_call",
                evidence_count=count,
                supporting_count=count,
                fraction=0.0,
                metric="KmerSutra evidence",
                detected=detected and bool(name),
                comparable=True,
                focus_taxa=focus_taxa,
            )
        )
    return [row for row in evidence if row["taxon_name"]]


def generate_html_reports(
    *,
    workflow: WorkflowConfig,
    sample_rows: Sequence[Mapping[str, Any]],
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
    warning_rows: Sequence[Mapping[str, Any]],
    final_root: Path,
    pcr_concordance_rows: Sequence[Mapping[str, Any]] = (),
    pcr_method_summary_rows: Sequence[Mapping[str, Any]] = (),
) -> tuple[Path, ...]:
    """Write the complete offline HTML report suite.

    Args:
        workflow: Validated run configuration.
        sample_rows: Host/input and classifier status summary rows.
        status_rows: One terminal status row per sample and classifier.
        evidence_rows: Normalised method-specific evidence rows.
        warning_rows: Non-fatal parsing and completeness warnings.
        final_root: Final-results directory.
        pcr_concordance_rows: Per-sample, per-method PCR comparisons.
        pcr_method_summary_rows: Exact PCR comparison counts by method.

    Returns:
        Every generated HTML path, ordered from the final report onwards.
    """
    report_root = final_root / "reports"
    classifier_root = report_root / "classifiers"
    sample_root = report_root / "samples"
    classifier_root.mkdir(parents=True, exist_ok=True)
    sample_root.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []
    final_path = report_root / "index.html"
    final_path.write_text(
        _render_final_report(
            workflow=workflow,
            sample_rows=sample_rows,
            status_rows=status_rows,
            evidence_rows=evidence_rows,
            warning_rows=warning_rows,
        ),
        encoding="utf-8",
    )
    outputs.append(final_path)

    comparison_path = report_root / "comparison.html"
    comparison_path.write_text(
        _render_comparison_report(
            workflow=workflow,
            status_rows=status_rows,
            evidence_rows=evidence_rows,
        ),
        encoding="utf-8",
    )
    outputs.append(comparison_path)

    pcr_path = report_root / "pcr_comparison.html"
    pcr_path.write_text(
        _render_pcr_report(
            workflow=workflow,
            concordance_rows=pcr_concordance_rows,
            method_summary_rows=pcr_method_summary_rows,
        ),
        encoding="utf-8",
    )
    outputs.append(pcr_path)

    for method in METHODS:
        path = classifier_root / f"{method}.html"
        path.write_text(
            _render_classifier_report(
                workflow=workflow,
                method=method,
                status_rows=status_rows,
                evidence_rows=evidence_rows,
            ),
            encoding="utf-8",
        )
        outputs.append(path)

    for sample in workflow.samples:
        path = sample_root / f"{sample.sample_id}.html"
        path.write_text(
            _render_sample_report(
                workflow=workflow,
                sample_id=sample.sample_id,
                status_rows=status_rows,
                evidence_rows=evidence_rows,
            ),
            encoding="utf-8",
        )
        outputs.append(path)

    manifest_path = final_root / "report_manifest.json"
    write_json_atomic(
        path=manifest_path,
        payload={
            "generated_at_utc": utc_now(),
            "run_id": workflow.run_id,
            "entry_point": str(final_path),
            "html_files": [str(path) for path in outputs],
            "offline": True,
            "external_assets": False,
        },
    )
    LOGGER.info("Generated %d offline HTML reports beneath %s", len(outputs), report_root)
    return tuple(outputs)


def _render_final_report(
    *,
    workflow: WorkflowConfig,
    sample_rows: Sequence[Mapping[str, Any]],
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
    warning_rows: Sequence[Mapping[str, Any]],
) -> str:
    statuses = [str(row.get("status", "missing")) for row in status_rows]
    success_count = sum(status in SUCCESS_STATUSES for status in statuses)
    incomplete_count = len(statuses) - success_count
    detected = _detected_rows(evidence_rows=evidence_rows)
    focus = [row for row in detected if _as_bool(row.get("is_focus"))]
    method_cards = "".join(
        _method_card(method=method, status_rows=status_rows, evidence_rows=evidence_rows)
        for method in METHODS
    )
    status_matrix = _status_matrix(
        workflow=workflow,
        status_rows=status_rows,
        prefix="",
    )
    focus_table = _table(
        rows=_sort_evidence(rows=focus)[: workflow.report_max_table_rows],
        columns=(
            ("sample_id", "Sample"),
            ("method", "Classifier"),
            ("taxon_name", "Taxon"),
            ("evidence_count", "Primary evidence"),
            ("metric", "Metric"),
        ),
        table_id="focus-evidence",
        empty_message="No configured focus-taxon evidence was detected.",
    )
    warnings = _warning_panel(warning_rows=warning_rows)
    sample_links = "".join(
        '<a class="sample-chip" href="samples/{0}.html">{1}</a>'.format(
            html.escape(sample.sample_id, quote=True),
            html.escape(sample.sample_id),
        )
        for sample in workflow.samples
    )
    body = f"""
    {
        _hero(
            eyebrow="Final run report",
            title=workflow.run_id,
            subtitle=(
                "An offline, failure-aware view of real Nanopore evidence across four "
                "independent classifiers."
            ),
        )
    }
    <section class="metric-grid" aria-label="Run overview">
      {_metric("Samples", len(workflow.samples), "classification-ready inputs")}
      {_metric("Successful branches", success_count, f"of {len(statuses)} sample-method runs")}
      {_metric("Incomplete branches", incomplete_count, "failed, timed out or unavailable")}
      {_metric("Focus evidence rows", len(focus), "method-specific; not consensus")}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Classifier health</span>
      <h2>What completed?</h2></div></div>
      <div class="classifier-grid">{method_cards}</div>
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Run matrix</span>
      <h2>Status by sample and classifier</h2></div>
      <a class="action" href="comparison.html">Open comparison →</a></div>
      {status_matrix}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Configured focus</span>
      <h2>{html.escape(", ".join(workflow.report_focus_taxa))} evidence</h2></div></div>
      {_methodology_note()}
      {focus_table}
    </section>
    {warnings}
    <section>
      <div class="section-heading"><div><span class="kicker">Sample reports</span>
      <h2>Open one sample</h2></div></div>
      <div class="sample-cloud">{sample_links}</div>
    </section>
    <section class="provenance">
      <h2>Run provenance</h2>
      <dl>
        <div><dt>Input state</dt><dd>{html.escape(workflow.input_read_state)}</dd></div>
        <div><dt>Samples sheet</dt><dd>{html.escape(str(workflow.samples_path))}</dd></div>
        <div><dt>Output root</dt><dd>{html.escape(str(workflow.output_directory))}</dd></div>
        <div><dt>Scratch</dt><dd>{html.escape(str(workflow.scratch_root))}</dd></div>
      </dl>
    </section>
    """
    return _document(
        title=f"{workflow.run_id} · final report",
        body=body,
        prefix="",
        active="final",
        footer_note=f"{len(sample_rows)} sample summaries · generated offline",
    )


def _render_pcr_report(
    *,
    workflow: WorkflowConfig,
    concordance_rows: Sequence[Mapping[str, Any]],
    method_summary_rows: Sequence[Mapping[str, Any]],
) -> str:
    """Render independent PCR concordance without forcing classifier consensus."""
    primary_samples = {
        str(row.get("sample_id", ""))
        for row in concordance_rows
        if _as_bool(row.get("include_in_primary_comparison"))
    }
    available = sum(
        int(_first_number(row=row, keys=("available_sample_count",))) for row in method_summary_rows
    )
    expected = sum(
        int(_first_number(row=row, keys=("primary_sample_count",))) for row in method_summary_rows
    )
    body = (
        _hero(
            eyebrow="Independent reference",
            title="PCR concordance by classifier",
            subtitle=(
                "Each method is compared separately with the supplied PCR interpretation; "
                "missing or failed classifier runs remain unavailable, not non-detections."
            ),
        )
        + '<section class="metric-grid">'
        + _metric("Primary PCR samples", len(primary_samples), "Exact sample denominator")
        + _metric("Available method results", available, f"of {expected} expected")
        + _metric(
            "Excluded PCR records",
            len(
                {
                    str(row.get("sample_id", ""))
                    for row in concordance_rows
                    if not _as_bool(row.get("include_in_primary_comparison"))
                }
            ),
            "Retained transparently",
        )
        + "</section>"
        + '<section class="panel"><h2>Method-level exact counts</h2>'
        + _table(
            table_id="pcr-method-summary",
            rows=method_summary_rows,
            columns=(
                ("method", "Classifier"),
                ("primary_sample_count", "Primary samples"),
                ("available_sample_count", "Available"),
                ("unavailable_sample_count", "Unavailable"),
                ("all_expected_species_detected_count", "All expected detected"),
                ("exact_species_match_count", "Exact species match"),
            ),
            empty_message="No PCR summary is available.",
        )
        + "</section>"
        + '<section class="panel"><h2>Per-sample comparison</h2>'
        + _table(
            table_id="pcr-concordance",
            rows=concordance_rows,
            columns=(
                ("sample_id", "Sample"),
                ("method", "Classifier"),
                ("pcr_species", "PCR species"),
                ("classifier_status", "Classifier status"),
                ("detected_expected_species", "Expected detected"),
                ("missed_expected_species", "Expected missed"),
                ("additional_plasmodium_species", "Additional Plasmodium"),
                ("comparison_status", "PCR comparison"),
            ),
            empty_message="No PCR concordance rows are available.",
            status_columns={"classifier_status"},
        )
        + "</section>"
    )
    return _document(
        title=f"{workflow.run_id} · PCR concordance",
        body=body,
        prefix="",
        active="pcr",
        footer_note="Independent PCR comparison · exact counts and denominators",
    )


def _render_classifier_report(
    *,
    workflow: WorkflowConfig,
    method: str,
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
) -> str:
    label = METHOD_LABELS[method]
    method_statuses = [row for row in status_rows if row.get("method") == method]
    method_evidence = [row for row in evidence_rows if row.get("method") == method]
    detected = _detected_rows(evidence_rows=method_evidence)
    top = _top_taxa(rows=detected, limit=workflow.report_top_n)
    successful = sum(row.get("status") == "success" for row in method_statuses)
    focus = sum(_as_bool(row.get("is_focus")) for row in detected)
    status_table = _table(
        rows=method_statuses,
        columns=(
            ("sample_id", "Sample"),
            ("status", "Status"),
            ("message", "Message"),
            ("completed_at_utc", "Completed (UTC)"),
        ),
        table_id=f"{method}-status",
        status_columns={"status"},
    )
    evidence_table = _table(
        rows=_sort_evidence(rows=method_evidence)[: workflow.report_max_table_rows],
        columns=(
            ("sample_id", "Sample"),
            ("taxon_name", "Taxon"),
            ("rank", "Rank / evidence type"),
            ("evidence_count", "Primary evidence"),
            ("supporting_count", "Supporting evidence"),
            ("fraction", "Report %"),
            ("detected", "Detected"),
            ("is_focus", "Focus taxon"),
        ),
        table_id=f"{method}-evidence",
        empty_message=f"No validated {label} evidence tables were available.",
    )
    body = f"""
    {
        _hero(
            eyebrow="Classifier report",
            title=label,
            subtitle=f"Method-specific results for {workflow.run_id}; no cross-tool consensus is imposed.",
        )
    }
    <section class="metric-grid">
      {_metric("Successful samples", successful, f"of {len(workflow.samples)}")}
      {_metric("Detected rows", len(detected), "using this method's native evidence")}
      {_metric("Focus rows", focus, html.escape(", ".join(workflow.report_focus_taxa)))}
      {_metric("Evidence records", len(method_evidence), "validated and reportable")}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Overview</span>
      <h2>Top reported taxa</h2></div></div>
      {_bar_chart(items=top, empty_message="No positive evidence was available for ranking.")}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Completeness</span>
      <h2>Sample execution status</h2></div></div>
      {status_table}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Evidence explorer</span>
      <h2>All retained {html.escape(label)} evidence</h2></div></div>
      {_methodology_note()}
      {evidence_table}
    </section>
    """
    return _document(
        title=f"{workflow.run_id} · {label}",
        body=body,
        prefix="../",
        active=method,
        footer_note=f"{label} classifier report",
    )


def _render_comparison_report(
    *,
    workflow: WorkflowConfig,
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
) -> str:
    detected = _detected_rows(evidence_rows=evidence_rows)
    pairwise = _pairwise_overlap(workflow=workflow, evidence_rows=detected)
    matrix = _comparison_matrix(evidence_rows=detected, max_rows=workflow.report_max_table_rows)
    comparison_table = _table(
        rows=pairwise,
        columns=(
            ("classifier_a", "Classifier A"),
            ("classifier_b", "Classifier B"),
            ("samples_compared", "Samples compared"),
            ("mean_jaccard", "Mean Jaccard"),
            ("shared_taxa", "Shared taxa"),
            ("union_taxa", "Union taxa"),
        ),
        table_id="pairwise-overlap",
    )
    matrix_table = _table(
        rows=matrix,
        columns=(
            ("sample_id", "Sample"),
            ("taxon_name", "Taxon"),
            ("classifier_count", "Classifiers"),
            ("kraken2", "Kraken2"),
            ("metabuli", "Metabuli"),
            ("minimap2", "minimap2"),
            ("kmersutra", "KmerSutra"),
            ("focus", "Focus"),
        ),
        table_id="comparison-evidence",
    )
    body = f"""
    {
        _hero(
            eyebrow="Cross-classifier comparison",
            title="Agreement without forced consensus",
            subtitle=(
                "A descriptive comparison of reportable taxa. Differences in database scope, "
                "algorithms and evidence units remain visible."
            ),
        )
    }
    <section>
      <div class="section-heading"><div><span class="kicker">Completeness first</span>
      <h2>Status matrix</h2></div></div>
      {_status_matrix(workflow=workflow, status_rows=status_rows, prefix="")}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Descriptive overlap</span>
      <h2>Pairwise species-level agreement</h2></div></div>
      <div class="callout"><strong>Interpret carefully.</strong> Jaccard overlap is calculated only
      from positive, species-comparable rows in samples where both methods succeeded. It is not an
      accuracy score and does not make one classifier the truth set.</div>
      {comparison_table}
    </section>
    <section>
      <div class="section-heading"><div><span class="kicker">Evidence matrix</span>
      <h2>Where methods agree—and where they do not</h2></div></div>
      {_methodology_note()}
      {matrix_table}
    </section>
    """
    return _document(
        title=f"{workflow.run_id} · classifier comparison",
        body=body,
        prefix="",
        active="comparison",
        footer_note="Descriptive cross-classifier comparison",
    )


def _render_sample_report(
    *,
    workflow: WorkflowConfig,
    sample_id: str,
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
) -> str:
    statuses = [row for row in status_rows if row.get("sample_id") == sample_id]
    evidence = [row for row in evidence_rows if row.get("sample_id") == sample_id]
    focus = [row for row in evidence if _as_bool(row.get("is_focus"))]
    status_table = _table(
        rows=statuses,
        columns=(("method", "Classifier"), ("status", "Status"), ("message", "Message")),
        table_id="sample-status",
        status_columns={"status"},
    )
    evidence_table = _table(
        rows=_sort_evidence(rows=evidence)[: workflow.report_max_table_rows],
        columns=(
            ("method", "Classifier"),
            ("taxon_name", "Taxon"),
            ("rank", "Rank / evidence type"),
            ("evidence_count", "Primary evidence"),
            ("metric", "Metric"),
            ("detected", "Detected"),
            ("is_focus", "Focus"),
        ),
        table_id="sample-evidence",
    )
    body = f"""
    {
        _hero(
            eyebrow="Sample report",
            title=sample_id,
            subtitle=f"Every available method-specific result for {workflow.run_id}.",
        )
    }
    <section class="metric-grid">
      {
        _metric(
            "Classifier results",
            sum(row.get("status") == "success" for row in statuses),
            "successful",
        )
    }
      {_metric("Evidence rows", len(evidence), "all retained ranks")}
      {_metric("Positive rows", len(_detected_rows(evidence_rows=evidence)), "method-specific")}
      {_metric("Focus rows", len(focus), html.escape(", ".join(workflow.report_focus_taxa)))}
    </section>
    <section><div class="section-heading"><div><span class="kicker">Completeness</span>
      <h2>Classifier status</h2></div></div>{status_table}</section>
    <section><div class="section-heading"><div><span class="kicker">Evidence explorer</span>
      <h2>Reported evidence</h2></div></div>{_methodology_note()}{evidence_table}</section>
    """
    return _document(
        title=f"{workflow.run_id} · {sample_id}",
        body=body,
        prefix="../",
        active="samples",
        footer_note=f"Sample report · {html.escape(sample_id)}",
    )


def _document(*, title: str, body: str, prefix: str, active: str, footer_note: str) -> str:
    return f"""<!doctype html>
<html lang="en-GB">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="colour-scheme" content="dark light">
  <title>{html.escape(title)}</title>
  <style>{_stylesheet()}</style>
</head>
<body>
  <a class="skip-link" href="#main">Skip to report</a>
  {_navigation(prefix=prefix, active=active)}
  <main id="main" class="shell">{body}</main>
  <footer><div class="shell"><span>{footer_note}</span>
    <span>Research use only · no external web assets · safe to archive</span></div></footer>
  <script>{_javascript()}</script>
</body>
</html>
"""


def _navigation(*, prefix: str, active: str) -> str:
    links = [
        ("final", f"{prefix}index.html", "Final report"),
        ("comparison", f"{prefix}comparison.html", "Compare"),
        ("pcr", f"{prefix}pcr_comparison.html", "PCR"),
        *(
            (method, f"{prefix}classifiers/{method}.html", METHOD_LABELS[method])
            for method in METHODS
        ),
    ]
    rendered = "".join(
        '<a class="{0}" href="{1}">{2}</a>'.format(
            "active" if key == active else "",
            html.escape(url, quote=True),
            html.escape(label),
        )
        for key, url, label in links
    )
    return (
        '<header class="topbar"><div class="shell nav-wrap">'
        '<a class="brand" href="{0}index.html"><span class="brand-mark">N</span>'
        "<span>Nanopore<br><small>classifier observatory</small></span></a>"
        '<nav aria-label="Report navigation">{1}</nav></div></header>'
    ).format(html.escape(prefix, quote=True), rendered)


def _hero(*, eyebrow: str, title: str, subtitle: str) -> str:
    return f"""
    <section class="hero">
      <div><span class="eyebrow">{html.escape(eyebrow)}</span>
      <h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></div>
      <div class="hero-orbit" aria-hidden="true"><span></span><span></span><span></span></div>
    </section>
    <div class="safety-banner"><strong>Research classification output.</strong>
    Method-specific evidence is shown without treating agreement as clinical confirmation.</div>
    """


def _metric(label: str, value: Any, detail: str) -> str:
    return (
        '<article class="metric"><span>{0}</span><strong>{1}</strong><small>{2}</small></article>'
    ).format(html.escape(label), html.escape(str(value)), detail)


def _method_card(
    *,
    method: str,
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
) -> str:
    statuses = [row for row in status_rows if row.get("method") == method]
    successes = sum(row.get("status") == "success" for row in statuses)
    incomplete = len(statuses) - successes
    evidence = [row for row in evidence_rows if row.get("method") == method]
    overall = "success" if incomplete == 0 and statuses else "partial"
    if statuses and successes == 0:
        overall = str(statuses[0].get("status", "failed"))
    return f"""
    <a class="classifier-card" href="classifiers/{method}.html">
      <div><span class="status {html.escape(overall)}">{html.escape(overall)}</span>
      <h3>{html.escape(METHOD_LABELS[method])}</h3></div>
      <strong>{successes}/{len(statuses)}</strong>
      <small>{len(evidence)} evidence rows · {incomplete} incomplete</small>
    </a>
    """


def _status_matrix(
    *,
    workflow: WorkflowConfig,
    status_rows: Sequence[Mapping[str, Any]],
    prefix: str,
) -> str:
    indexed = {(str(row.get("sample_id")), str(row.get("method"))): row for row in status_rows}
    headers = "".join(f"<th>{html.escape(METHOD_LABELS[method])}</th>" for method in METHODS)
    rows = []
    for sample in workflow.samples:
        cells = []
        for method in METHODS:
            row = indexed.get((sample.sample_id, method), {})
            status = str(row.get("status", "missing"))
            message = str(row.get("message", ""))
            cells.append(
                '<td><span class="status {0}" title="{1}">{0}</span></td>'.format(
                    html.escape(status, quote=True),
                    html.escape(message, quote=True),
                )
            )
        rows.append(
            '<tr><th><a href="{0}samples/{1}.html">{2}</a></th>{3}</tr>'.format(
                html.escape(prefix, quote=True),
                html.escape(sample.sample_id, quote=True),
                html.escape(sample.sample_id),
                "".join(cells),
            )
        )
    return (
        '<div class="table-wrap"><table class="status-matrix"><thead><tr><th>Sample</th>'
        f"{headers}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )


def _table(
    *,
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[tuple[str, str]],
    table_id: str,
    empty_message: str = "No rows were available.",
    status_columns: set[str] | None = None,
) -> str:
    if not rows:
        return f'<div class="empty-state">{html.escape(empty_message)}</div>'
    status_columns = status_columns or set()
    headers = "".join(
        '<th scope="col"><button type="button" data-sort="{0}">{1}<span>↕</span></button></th>'.format(
            index,
            html.escape(label),
        )
        for index, (_, label) in enumerate(columns)
    )
    rendered_rows = []
    for row in rows:
        cells = []
        searchable = " ".join(str(row.get(key, "")) for key, _ in columns)
        for key, _ in columns:
            value = row.get(key, "")
            text = _display_value(value=value)
            if key in status_columns:
                status = str(value or "missing")
                cells.append(
                    f'<td><span class="status {html.escape(status, quote=True)}">'
                    f"{html.escape(status)}</span></td>"
                )
            elif isinstance(value, bool):
                cells.append(f'<td><span class="boolean">{str(value).lower()}</span></td>')
            else:
                cells.append(f"<td>{html.escape(text)}</td>")
        rendered_rows.append(
            f'<tr data-search="{html.escape(searchable.casefold(), quote=True)}">'
            + "".join(cells)
            + "</tr>"
        )
    return f"""
    <div class="table-tools">
      <label>Filter <input type="search" data-filter="{html.escape(table_id, quote=True)}"
      placeholder="sample, taxon, status…"></label>
      <button type="button" data-export="{html.escape(table_id, quote=True)}">Export visible TSV</button>
      <span data-count-for="{html.escape(table_id, quote=True)}">{len(rows)} rows</span>
    </div>
    <div class="table-wrap"><table id="{html.escape(table_id, quote=True)}" class="data-table">
      <thead><tr>{headers}</tr></thead><tbody>{"".join(rendered_rows)}</tbody>
    </table></div>
    """


def _bar_chart(*, items: Sequence[tuple[str, float]], empty_message: str) -> str:
    if not items:
        return f'<div class="empty-state">{html.escape(empty_message)}</div>'
    maximum = max(value for _, value in items) or 1.0
    bars = []
    for name, value in items:
        width = max(1.0, 100.0 * value / maximum)
        bars.append(
            '<div class="bar-row"><span title="{0}">{0}</span><div><i style="width:{1:.3f}%"></i>'
            "</div><strong>{2}</strong></div>".format(
                html.escape(name),
                width,
                html.escape(_display_value(value=value)),
            )
        )
    return f'<div class="bar-chart">{"".join(bars)}</div>'


def _methodology_note() -> str:
    return (
        '<div class="method-note"><span>i</span><p><strong>Evidence stays method-specific.</strong> '
        "Kraken2 and Metabuli counts, controlled-reference alignments and KmerSutra exact-k-mer "
        "calls are not interchangeable measurements. Tables expose the native evidence type.</p></div>"
    )


def _warning_panel(*, warning_rows: Sequence[Mapping[str, Any]]) -> str:
    if not warning_rows:
        return ""
    items = "".join(
        "<li><strong>{0}</strong> · {1}</li>".format(
            html.escape(str(row.get("source", "reporting"))),
            html.escape(str(row.get("message", "unspecified warning"))),
        )
        for row in warning_rows
    )
    return f"""
    <section class="warning-panel"><div><span class="kicker">Non-fatal warnings</span>
      <h2>What needs attention</h2></div><ul>{items}</ul></section>
    """


def _pairwise_overlap(
    *,
    workflow: WorkflowConfig,
    evidence_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in evidence_rows:
        if not _as_bool(row.get("comparable")):
            continue
        key = (str(row.get("sample_id", "")), str(row.get("method", "")))
        sets[key].add(str(row.get("taxon_name", "")).casefold())
    results = []
    for first_index, first in enumerate(METHODS):
        for second in METHODS[first_index + 1 :]:
            scores: list[float] = []
            shared_total = 0
            union_total = 0
            for sample in workflow.samples:
                left = sets[(sample.sample_id, first)]
                right = sets[(sample.sample_id, second)]
                union = left | right
                if not union:
                    continue
                intersection = left & right
                scores.append(len(intersection) / len(union))
                shared_total += len(intersection)
                union_total += len(union)
            results.append(
                {
                    "classifier_a": METHOD_LABELS[first],
                    "classifier_b": METHOD_LABELS[second],
                    "samples_compared": len(scores),
                    "mean_jaccard": f"{sum(scores) / len(scores):.3f}" if scores else "NA",
                    "shared_taxa": shared_total,
                    "union_taxa": union_total,
                }
            )
    return results


def _comparison_matrix(
    *,
    evidence_rows: Sequence[Mapping[str, Any]],
    max_rows: int,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for row in evidence_rows:
        if not _as_bool(row.get("comparable")):
            continue
        sample_id = str(row.get("sample_id", ""))
        taxon_name = str(row.get("taxon_name", ""))
        key = (sample_id, taxon_name.casefold())
        record = grouped.setdefault(
            key,
            {
                "sample_id": sample_id,
                "taxon_name": taxon_name,
                "kraken2": "—",
                "metabuli": "—",
                "minimap2": "—",
                "kmersutra": "—",
                "focus": "yes" if _as_bool(row.get("is_focus")) else "no",
            },
        )
        method = str(row.get("method", ""))
        if method in METHODS:
            record[method] = _display_value(value=row.get("evidence_count", 0))
        if _as_bool(row.get("is_focus")):
            record["focus"] = "yes"
    for record in grouped.values():
        record["classifier_count"] = sum(record[method] != "—" for method in METHODS)
    return sorted(
        grouped.values(),
        key=lambda row: (
            row["focus"] != "yes",
            -int(row["classifier_count"]),
            str(row["sample_id"]),
            str(row["taxon_name"]),
        ),
    )[:max_rows]


def _top_taxa(
    *,
    rows: Sequence[Mapping[str, Any]],
    limit: int,
) -> list[tuple[str, float]]:
    totals: dict[str, float] = defaultdict(float)
    labels: dict[str, str] = {}
    for row in rows:
        name = str(row.get("taxon_name", "")).strip()
        if not name:
            continue
        key = name.casefold()
        labels.setdefault(key, name)
        totals[key] += _non_negative_number(row.get("evidence_count"))
    return [
        (labels[key], value)
        for key, value in sorted(totals.items(), key=lambda item: (-item[1], item[0]))[:limit]
    ]


def _sort_evidence(*, rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            not _as_bool(row.get("is_focus")),
            not _as_bool(row.get("detected")),
            str(row.get("sample_id", "")),
            str(row.get("method", "")),
            -_non_negative_number(row.get("evidence_count")),
            str(row.get("taxon_name", "")),
        ),
    )


def _detected_rows(
    *,
    evidence_rows: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    return [row for row in evidence_rows if _as_bool(row.get("detected"))]


def _evidence_row(
    *,
    row: Mapping[str, Any],
    taxon_name: str,
    tax_id: str,
    rank: str,
    evidence_count: float,
    supporting_count: float,
    fraction: float,
    metric: str,
    detected: bool,
    comparable: bool,
    focus_taxa: Sequence[str],
) -> dict[str, Any]:
    return {
        "sample_id": str(row.get("sample_id", "")),
        "method": str(row.get("method", "")),
        "taxon_name": taxon_name,
        "tax_id": tax_id,
        "rank": rank,
        "evidence_count": evidence_count,
        "supporting_count": supporting_count,
        "fraction": fraction,
        "metric": metric,
        "detected": detected,
        "comparable": comparable,
        "is_focus": any(term.casefold() in taxon_name.casefold() for term in focus_taxa),
    }


def _first_text(*, row: Mapping[str, Any], keys: Sequence[str]) -> str:
    for key in keys:
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def _first_number(*, row: Mapping[str, Any], keys: Sequence[str]) -> float:
    for key in keys:
        if str(row.get(key, "")).strip():
            return _non_negative_number(row.get(key))
    return 0.0


def _non_negative_number(value: Any) -> float:
    try:
        numeric = float(str(value).strip())
    except (TypeError, ValueError):
        return 0.0
    return max(0.0, numeric)


def _call_is_positive(*, value: str) -> bool:
    normalised = value.strip().casefold().replace("-", "_").replace(" ", "_")
    positive = {
        "1",
        "detected",
        "mixed_species_present",
        "positive",
        "present",
        "present_high_confidence",
        "present_in_mixed_sample",
        "reportable",
        "true",
        "yes",
    }
    return normalised in positive


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().casefold() in {"1", "true", "yes", "detected", "positive"}


def _display_value(*, value: Any) -> str:
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.4g}"
    if value is None:
        return ""
    return str(value)


def _stylesheet() -> str:
    return """
:root{--bg:#07111f;--panel:#0d1b2d;--panel2:#10243a;--ink:#eaf2f8;--muted:#9fb2c4;
--cyan:#5ee7df;--blue:#45a3ff;--violet:#9b7bff;--amber:#ffca69;--red:#ff7185;
--green:#69e6a6;--line:rgba(168,196,220,.18);--shadow:0 24px 70px rgba(0,0,0,.25)}
*{box-sizing:border-box}html{scroll-behaviour:smooth}body{margin:0;background:
radial-gradient(circle at 84% 8%,rgba(69,163,255,.16),transparent 28rem),
radial-gradient(circle at 5% 38%,rgba(155,123,255,.10),transparent 24rem),var(--bg);
color:var(--ink);font:15px/1.55 Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:var(--cyan);text-decoration:none}a:hover{text-decoration:underline}.shell{width:min(1500px,94vw);margin:auto}
.skip-link{position:absolute;left:-10000px;top:auto}.skip-link:focus{left:1rem;top:1rem;z-index:99;
background:#fff;color:#000;padding:.7rem 1rem;border-radius:.5rem}.topbar{position:sticky;top:0;z-index:20;
background:rgba(7,17,31,.88);backdrop-filter:blur(18px);border-bottom:1px solid var(--line)}
.nav-wrap{display:flex;align-items:center;gap:2rem;min-height:72px}.brand{display:flex;align-items:center;gap:.7rem;
font-weight:750;color:var(--ink);letter-spacing:.01em;white-space:nowrap}.brand:hover{text-decoration:none}
.brand small{color:var(--muted);font-weight:500}.brand-mark{display:grid;place-items:center;width:36px;height:36px;
border-radius:12px;background:linear-gradient(135deg,var(--cyan),var(--blue));color:#06111d;font-size:1.15rem}
nav{display:flex;gap:.25rem;overflow:auto;padding:.5rem 0}nav a{color:var(--muted);padding:.55rem .72rem;
border-radius:.55rem;white-space:nowrap;font-size:.88rem}nav a.active,nav a:hover{background:rgba(94,231,223,.1);
color:var(--ink);text-decoration:none}.hero{min-height:340px;display:flex;justify-content:space-between;align-items:center;
gap:3rem;padding:5.8rem 0 3rem}.eyebrow,.kicker{color:var(--cyan);font-size:.73rem;font-weight:800;
letter-spacing:.17em;text-transform:uppercase}.hero h1{font-size:clamp(2.7rem,6vw,5.8rem);line-height:.98;
letter-spacing:-.055em;margin:.45rem 0 1.25rem;max-width:1050px}.hero p{max-width:760px;color:var(--muted);
font-size:1.15rem}.hero-orbit{position:relative;width:190px;height:190px;border:1px solid rgba(94,231,223,.25);
border-radius:50%;flex:0 0 auto}.hero-orbit:before,.hero-orbit:after{content:"";position:absolute;inset:23px;
border:1px solid rgba(69,163,255,.25);border-radius:50%}.hero-orbit:after{inset:53px;background:linear-gradient(135deg,
rgba(94,231,223,.65),rgba(69,163,255,.12));box-shadow:0 0 70px rgba(94,231,223,.25)}
.hero-orbit span{position:absolute;width:12px;height:12px;border-radius:50%;background:var(--cyan);
box-shadow:0 0 18px var(--cyan)}.hero-orbit span:nth-child(1){left:7px;top:86px}.hero-orbit span:nth-child(2){right:22px;top:23px;background:var(--violet)}
.hero-orbit span:nth-child(3){right:24px;bottom:20px;background:var(--blue)}.safety-banner,.callout,
.method-note{border:1px solid rgba(255,202,105,.3);background:rgba(255,202,105,.075);padding:1rem 1.2rem;
border-radius:14px;color:#f5dfb8}.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1rem;
margin:1.2rem 0 5rem}.metric{padding:1.3rem 1.4rem;border:1px solid var(--line);background:linear-gradient(145deg,
rgba(16,36,58,.92),rgba(10,24,40,.82));border-radius:16px;box-shadow:var(--shadow)}.metric span,.metric small{
display:block;color:var(--muted)}.metric strong{display:block;font-size:2rem;line-height:1.2;margin:.25rem 0}
section{margin:0 0 5.5rem}.section-heading{display:flex;align-items:end;justify-content:space-between;gap:1rem;
margin-bottom:1.4rem}.section-heading h2,.warning-panel h2,.provenance h2{font-size:clamp(1.7rem,3vw,2.6rem);
margin:.25rem 0;letter-spacing:-.035em}.action,.table-tools button{border:1px solid var(--line);padding:.55rem .85rem;
border-radius:.55rem;background:var(--panel);color:var(--cyan);cursor:pointer}.classifier-grid{display:grid;
grid-template-columns:repeat(4,minmax(0,1fr));gap:1rem}.classifier-card{display:grid;min-height:190px;padding:1.35rem;
border:1px solid var(--line);border-radius:16px;background:linear-gradient(145deg,var(--panel2),var(--panel));
color:var(--ink);transition:.2s transform,.2s border-color}.classifier-card:hover{transform:translateY(-3px);
border-color:rgba(94,231,223,.45);text-decoration:none}.classifier-card h3{font-size:1.5rem;margin:.7rem 0}
.classifier-card>strong{font-size:2.4rem;align-self:end}.classifier-card small{color:var(--muted)}.status{display:inline-flex;
align-items:center;width:max-content;border-radius:999px;padding:.18rem .55rem;font-weight:780;font-size:.72rem;
text-transform:uppercase;letter-spacing:.06em;background:rgba(159,178,196,.14);color:var(--muted)}
.status.success{background:rgba(105,230,166,.14);color:var(--green)}.status.partial,.status.timeout{
background:rgba(255,202,105,.14);color:var(--amber)}.status.failed,.status.invalid,.status.missing,
.status.scheduler_failed{
background:rgba(255,113,133,.14);color:var(--red)}.status.disabled,.status.skipped,.status.unavailable{
background:rgba(155,123,255,.14);color:#c7b5ff}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:14px;
background:rgba(13,27,45,.82)}table{width:100%;border-collapse:collapse;min-width:720px}th,td{text-align:left;
padding:.75rem .85rem;border-bottom:1px solid var(--line);vertical-align:top}thead th{position:sticky;top:0;
background:#10243a;z-index:2;font-size:.78rem;color:#bfd0df;text-transform:uppercase;letter-spacing:.05em}
tbody tr:hover{background:rgba(94,231,223,.045)}th button{all:unset;cursor:pointer;display:flex;gap:.5rem;align-items:center}
th button span{color:var(--muted)}.status-matrix td{text-align:center}.status-matrix th a{font-weight:700}
.table-tools{display:flex;align-items:center;gap:.8rem;margin:.8rem 0;flex-wrap:wrap}.table-tools label{color:var(--muted)}
.table-tools input{margin-left:.4rem;background:var(--panel);color:var(--ink);border:1px solid var(--line);
border-radius:.55rem;padding:.55rem .75rem;min-width:280px}.table-tools span{margin-left:auto;color:var(--muted)}
.method-note{display:flex;gap:.8rem;align-items:flex-start;border-color:rgba(94,231,223,.22);
background:rgba(94,231,223,.05);color:#c8e8e7;margin:0 0 1rem}.method-note>span{display:grid;place-items:center;
width:24px;height:24px;border-radius:50%;background:var(--cyan);color:#07111f;font-weight:900;flex:0 0 auto}
.method-note p{margin:0}.bar-chart{display:grid;gap:.65rem}.bar-row{display:grid;grid-template-columns:minmax(170px,1fr)
minmax(240px,4fr) 90px;gap:1rem;align-items:center}.bar-row>span{white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.bar-row>div{height:12px;background:rgba(159,178,196,.12);border-radius:99px;overflow:hidden}.bar-row i{display:block;
height:100%;border-radius:99px;background:linear-gradient(90deg,var(--cyan),var(--blue));box-shadow:0 0 18px rgba(94,231,223,.25)}
.bar-row strong{text-align:right;font-variant-numeric:tabular-nums}.empty-state{padding:2rem;border:1px dashed var(--line);
border-radius:14px;color:var(--muted);text-align:center}.warning-panel{border:1px solid rgba(255,113,133,.3);
background:rgba(255,113,133,.06);padding:1.5rem;border-radius:16px}.warning-panel ul{margin:.8rem 0 0;padding-left:1.4rem}
.sample-cloud{display:flex;flex-wrap:wrap;gap:.7rem}.sample-chip{padding:.65rem .9rem;border:1px solid var(--line);
border-radius:999px;background:var(--panel)}.provenance{padding:1.5rem;border:1px solid var(--line);border-radius:16px;
background:var(--panel)}.provenance dl{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem}
.provenance dl div{min-width:0}.provenance dt{color:var(--muted);font-size:.75rem;text-transform:uppercase;
letter-spacing:.08em}.provenance dd{margin:.25rem 0;overflow-wrap:anywhere}footer{border-top:1px solid var(--line);color:var(--muted);
padding:1.4rem 0;margin-top:6rem}footer .shell{display:flex;justify-content:space-between;gap:1rem}.boolean{color:var(--cyan)}
@media(max-width:980px){.metric-grid,.classifier-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.hero-orbit{display:none}}
@media(max-width:620px){.nav-wrap{display:block;padding:.65rem 0}.brand{margin-bottom:.4rem}.hero{padding:3.2rem 0 2rem;
min-height:auto}.metric-grid,.classifier-grid,.provenance dl{grid-template-columns:1fr}.bar-row{grid-template-columns:1fr}.bar-row strong{text-align:left}
.table-tools input{min-width:180px}.table-tools span{margin-left:0}footer .shell{display:block}}
@media print{.topbar,.table-tools,.hero-orbit{display:none}body{background:#fff;color:#111}.shell{width:100%}
.metric,.classifier-card,.table-wrap,.provenance{box-shadow:none;background:#fff;color:#111;border-color:#aaa}a{color:#111}}
"""


def _javascript() -> str:
    return """
(() => {
  const clean = value => (value || '').toString().trim();
  document.querySelectorAll('[data-filter]').forEach(input => {
    const table = document.getElementById(input.dataset.filter);
    if (!table) return;
    const rows = [...table.tBodies[0].rows];
    const counter = document.querySelector(`[data-count-for="${table.id}"]`);
    input.addEventListener('input', () => {
      const query = input.value.trim().toLocaleLowerCase('en-GB');
      let shown = 0;
      rows.forEach(row => {
        const visible = !query || (row.dataset.search || '').includes(query);
        row.hidden = !visible;
        if (visible) shown += 1;
      });
      if (counter) counter.textContent = `${shown} of ${rows.length} rows`;
    });
  });
  document.querySelectorAll('[data-sort]').forEach(button => {
    button.addEventListener('click', () => {
      const table = button.closest('table');
      const body = table.tBodies[0];
      const index = Number(button.dataset.sort);
      const direction = button.dataset.direction === 'asc' ? 'desc' : 'asc';
      button.dataset.direction = direction;
      [...body.rows].sort((left, right) => {
        const a = clean(left.cells[index].innerText);
        const b = clean(right.cells[index].innerText);
        const an = Number(a.replace(/,/g, ''));
        const bn = Number(b.replace(/,/g, ''));
        const compared = Number.isFinite(an) && Number.isFinite(bn)
          ? an - bn : a.localeCompare(b, 'en-GB', {numeric: true});
        return direction === 'asc' ? compared : -compared;
      }).forEach(row => body.appendChild(row));
    });
  });
  document.querySelectorAll('[data-export]').forEach(button => {
    button.addEventListener('click', () => {
      const table = document.getElementById(button.dataset.export);
      if (!table) return;
      const lines = [[...table.tHead.rows[0].cells].map(cell => clean(cell.innerText).replace('↕',''))];
      [...table.tBodies[0].rows].filter(row => !row.hidden).forEach(row => {
        lines.push([...row.cells].map(cell => clean(cell.innerText).replace(/[\t\r\n]+/g, ' ')));
      });
      const blob = new Blob([lines.map(row => row.join('\t')).join('\n') + '\n'],
        {type: 'text/tab-separated-values;charset=utf-8'});
      const link = document.createElement('a');
      link.href = URL.createObjectURL(blob);
      link.download = `${table.id}.visible.tsv`;
      link.click();
      URL.revokeObjectURL(link.href);
    });
  });
})();
"""


def serialisable_report_data(
    *,
    sample_rows: Sequence[Mapping[str, Any]],
    status_rows: Sequence[Mapping[str, Any]],
    evidence_rows: Sequence[Mapping[str, Any]],
    warning_rows: Sequence[Mapping[str, Any]],
    pcr_truth_rows: Sequence[Mapping[str, Any]] = (),
    pcr_concordance_rows: Sequence[Mapping[str, Any]] = (),
    pcr_method_summary_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Return JSON-safe report data for independent downstream visualisation."""
    payload = {
        "sample_summary": [dict(row) for row in sample_rows],
        "classifier_status": [dict(row) for row in status_rows],
        "normalised_evidence": [dict(row) for row in evidence_rows],
        "warnings": [dict(row) for row in warning_rows],
        "pcr_truth": [dict(row) for row in pcr_truth_rows],
        "pcr_concordance": [dict(row) for row in pcr_concordance_rows],
        "pcr_method_summary": [dict(row) for row in pcr_method_summary_rows],
    }
    # A defensive round-trip prevents unsupported custom values leaking into
    # the durable report-data contract.
    return json.loads(json.dumps(payload, ensure_ascii=False, default=str))
