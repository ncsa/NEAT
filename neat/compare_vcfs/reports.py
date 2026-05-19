"""
Report generation for `neat compare-vcfs`.

Three artifacts:
  - comparison_summary.json — machine-readable rollup
  - comparison_summary.txt  — human-readable rollup
  - FN_with_reasons.vcf     — hap.py's FN records, annotated with NEAT_REASON
"""
import json
import logging
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import pysam

from .. import __version__ as NEAT_VERSION

__all__ = [
    "REPORT_SCHEMA_VERSION",
    "build_comparison_summary",
    "compute_metrics",
    "render_summary_txt",
    "summarize_fn_reasons",
    "write_comparison_summary_json",
    "write_comparison_summary_txt",
    "write_fn_attribution_plot",
    "write_fn_with_reasons",
]

_LOG = logging.getLogger(__name__)

REPORT_SCHEMA_VERSION = "1"


def compute_metrics(counts: dict) -> dict:
    """
    Precision, recall, F1 from TP/FN/FP. Returns None for any metric that is
    undefined (e.g., empty truth set yields undefined recall).
    """
    tp = counts.get("TP", 0)
    fn = counts.get("FN", 0)
    fp = counts.get("FP", 0)
    precision = tp / (tp + fp) if (tp + fp) > 0 else None
    recall = tp / (tp + fn) if (tp + fn) > 0 else None
    if precision is None or recall is None or (precision + recall) == 0:
        f1 = None
    else:
        f1 = 2 * precision * recall / (precision + recall)
    return {"precision": precision, "recall": recall, "f1": f1}


def summarize_fn_reasons(fn_reasons) -> dict[str, int]:
    """Roll up reason tags across all FNs to a {reason: count} dict."""
    counts: dict[str, int] = {}
    for _, reasons in fn_reasons:
        for r in reasons:
            counts[r] = counts.get(r, 0) + 1
    return counts


def build_comparison_summary(
    *,
    golden_vcf: Path,
    called_vcf: Path,
    neat_run_dir: Path,
    simulation_summary_path: Path,
    happy_output_vcf: Path,
    happy_output_prefix: Path,
    counts: dict,
    fn_attribution: dict,
    fn_with_reasons_vcf: Path,
    comparison_summary_json: Path,
    comparison_summary_txt: Path,
) -> dict:
    """Assemble the comparison_summary dict from the run's artifacts."""
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "neat_version": NEAT_VERSION,
        "generated_at": _iso_utc(time.time()),
        "inputs": {
            "golden_vcf":          str(Path(golden_vcf).resolve()),
            "called_vcf":          str(Path(called_vcf).resolve()),
            "neat_run_dir":        str(Path(neat_run_dir).resolve()),
            "simulation_summary":  str(Path(simulation_summary_path).resolve()),
        },
        "happy": {
            "output_prefix": str(Path(happy_output_prefix).resolve()),
            "output_vcf":    str(Path(happy_output_vcf).resolve()),
        },
        "counts":         dict(counts),
        "metrics":        compute_metrics(counts),
        "fn_attribution": dict(fn_attribution),
        "outputs": {
            "fn_with_reasons_vcf":      str(Path(fn_with_reasons_vcf).resolve()),
            "comparison_summary_json":  str(Path(comparison_summary_json).resolve()),
            "comparison_summary_txt":   str(Path(comparison_summary_txt).resolve()),
        },
    }


def write_comparison_summary_json(summary: dict, path: Path) -> Path:
    """Atomically write the summary as pretty JSON."""
    path = Path(path)
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(summary, fh, indent=2)
    os.replace(tmp, path)
    return path


def write_comparison_summary_txt(summary: dict, path: Path) -> Path:
    """Atomically write a human-readable rollup."""
    path = Path(path)
    tmp = path.with_suffix(".txt.tmp")
    tmp.write_text(render_summary_txt(summary))
    os.replace(tmp, path)
    return path


def render_summary_txt(summary: dict) -> str:
    """Render the summary dict as a fixed-width text report."""
    counts = summary["counts"]
    metrics = summary["metrics"]
    fn_attr = summary["fn_attribution"]
    inputs = summary["inputs"]
    outputs = summary["outputs"]

    lines: list[str] = [
        "NEAT compare-vcfs report",
        "========================",
        "",
        f"Generated:    {summary['generated_at']}",
        f"NEAT version: {summary['neat_version']}",
        "",
        "Inputs",
        "------",
        f"  Truth (golden) VCF: {inputs['golden_vcf']}",
        f"  Called VCF:         {inputs['called_vcf']}",
        f"  NEAT run dir:       {inputs['neat_run_dir']}",
        "",
        "Classification (from hap.py)",
        "----------------------------",
        f"  True positives  (TP): {counts.get('TP', 0)}",
        f"  False negatives (FN): {counts.get('FN', 0)}",
        f"  False positives (FP): {counts.get('FP', 0)}",
        "",
        "Metrics",
        "-------",
        f"  Precision: {_fmt_metric(metrics['precision'])}  (TP / (TP + FP))",
        f"  Recall:    {_fmt_metric(metrics['recall'])}  (TP / (TP + FN))",
        f"  F1:        {_fmt_metric(metrics['f1'])}",
        "",
        "FN attribution",
        "--------------",
    ]
    if fn_attr:
        width = max(len(k) for k in fn_attr) + 2
        for reason in sorted(fn_attr):
            lines.append(f"  {reason:<{width}} {fn_attr[reason]}")
    else:
        lines.append("  (no false negatives)")

    lines += [
        "",
        "Outputs",
        "-------",
        f"  Annotated FN VCF:   {outputs['fn_with_reasons_vcf']}",
        f"  This report (JSON): {outputs['comparison_summary_json']}",
        f"  This report (text): {outputs['comparison_summary_txt']}",
        "",
    ]
    return "\n".join(lines)


def write_fn_with_reasons(
    source_happy_vcf: Path,
    fn_reasons: list,
    output_path: Path,
) -> Path:
    """
    Write only FN records to `output_path`, each annotated with a NEAT_REASON
    INFO tag. The source VCF is opened only to obtain a compatible header;
    records come from the already-classified `fn_reasons` list.
    """
    source_happy_vcf = Path(source_happy_vcf)
    output_path = Path(output_path)
    with pysam.VariantFile(str(source_happy_vcf)) as src:
        new_header = src.header.copy()
    new_header.info.add(
        "NEAT_REASON", ".", "String",
        "Comma-separated NEAT-aware false-negative attribution reasons",
    )
    with pysam.VariantFile(str(output_path), "w", header=new_header) as dst:
        for rec, reasons in fn_reasons:
            rec.translate(new_header)
            rec.info["NEAT_REASON"] = ",".join(reasons)
            dst.write(rec)
    return output_path


def write_fn_attribution_plot(fn_attribution: dict[str, int], path: Path) -> Path:
    """
    Render a horizontal bar chart of FN-reason counts to `path` (PNG).

    A true Venn diagram is awkward with this scheme (outside_simulated_contigs
    is mutually exclusive with the BED reasons; unknown can't co-occur with
    anything), so the artifact is a clearer bar chart.
    """
    # Import here so users running without --plot don't pay the matplotlib import cost
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(path)
    reasons = sorted(fn_attribution)
    counts = [fn_attribution[r] for r in reasons]

    fig, ax = plt.subplots(figsize=(8, max(2.0, 0.5 * len(reasons) + 1.0)))
    if reasons:
        ax.barh(reasons, counts, color="steelblue")
        for i, count in enumerate(counts):
            ax.text(count, i, f" {count}", va="center")
        ax.set_xlim(0, max(counts) * 1.15 if max(counts) else 1)
    else:
        ax.text(0.5, 0.5, "no false negatives", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_xlabel("FN count")
    ax.set_title("False-negative attribution")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def _iso_utc(epoch_seconds: float) -> str:
    return datetime.fromtimestamp(epoch_seconds, tz=timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _fmt_metric(v) -> str:
    return "N/A" if v is None else f"{v:.4f}"
