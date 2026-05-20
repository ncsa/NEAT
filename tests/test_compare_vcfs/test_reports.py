"""
Tests for neat/compare_vcfs/reports.py — metric computation, summary builders,
text rendering, and the annotated-FN VCF writer.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from neat.compare_vcfs.attribution import (
    REASON_OUTSIDE_MUTATION_BED,
    REASON_OUTSIDE_TARGET_BED,
    REASON_UNKNOWN,
)
from neat.compare_vcfs.reports import (
    REPORT_SCHEMA_VERSION,
    build_comparison_summary,
    compute_metrics,
    render_summary_txt,
    summarize_fn_reasons,
    write_comparison_summary_json,
    write_comparison_summary_txt,
    write_fn_attribution_plot,
    write_fn_with_reasons,
)


# ===========================================================================
# compute_metrics
# ===========================================================================

def test_compute_metrics_typical_counts():
    m = compute_metrics({"TP": 90, "FN": 10, "FP": 5})
    assert m["precision"] == pytest.approx(90 / 95)
    assert m["recall"] == pytest.approx(0.9)
    p, r = m["precision"], m["recall"]
    assert m["f1"] == pytest.approx(2 * p * r / (p + r))


def test_compute_metrics_perfect_call():
    m = compute_metrics({"TP": 100, "FN": 0, "FP": 0})
    assert m["precision"] == 1.0
    assert m["recall"] == 1.0
    assert m["f1"] == 1.0


def test_compute_metrics_no_truth_no_calls_returns_none_metrics():
    """Empty buckets → precision/recall/f1 all undefined."""
    m = compute_metrics({"TP": 0, "FN": 0, "FP": 0})
    assert m == {"precision": None, "recall": None, "f1": None}


def test_compute_metrics_no_calls_means_precision_undefined():
    m = compute_metrics({"TP": 0, "FN": 5, "FP": 0})
    assert m["precision"] is None
    assert m["recall"] == 0.0  # tp / (tp+fn) is defined (0 / 5)
    assert m["f1"] is None


def test_compute_metrics_no_truth_means_recall_undefined():
    m = compute_metrics({"TP": 0, "FN": 0, "FP": 5})
    assert m["precision"] == 0.0
    assert m["recall"] is None
    assert m["f1"] is None


# ===========================================================================
# summarize_fn_reasons
# ===========================================================================

def test_summarize_fn_reasons_counts_each_tag():
    fn_reasons = [
        (SimpleNamespace(), [REASON_UNKNOWN]),
        (SimpleNamespace(), [REASON_OUTSIDE_TARGET_BED]),
        (SimpleNamespace(), [REASON_OUTSIDE_TARGET_BED, REASON_OUTSIDE_MUTATION_BED]),
    ]
    assert summarize_fn_reasons(fn_reasons) == {
        REASON_UNKNOWN: 1,
        REASON_OUTSIDE_TARGET_BED: 2,
        REASON_OUTSIDE_MUTATION_BED: 1,
    }


def test_summarize_fn_reasons_empty_input():
    assert summarize_fn_reasons([]) == {}


# ===========================================================================
# build_comparison_summary
# ===========================================================================

_MISSING = object()


def _build_minimal_summary(tmp_path, counts=_MISSING, fn_attr=_MISSING):
    if counts is _MISSING:
        counts = {"TP": 1, "FN": 1, "FP": 1}
    if fn_attr is _MISSING:
        fn_attr = {REASON_UNKNOWN: 1}
    return build_comparison_summary(
        golden_vcf=tmp_path / "g.vcf",
        called_vcf=tmp_path / "c.vcf",
        neat_run_dir=tmp_path / "run",
        simulation_summary_path=tmp_path / "run" / "simulation_summary.json",
        happy_output_vcf=tmp_path / "out" / "happy.vcf.gz",
        happy_output_prefix=tmp_path / "out" / "happy",
        counts=counts,
        fn_attribution=fn_attr,
        fn_with_reasons_vcf=tmp_path / "out" / "FN_with_reasons.vcf",
        comparison_summary_json=tmp_path / "out" / "comparison_summary.json",
        comparison_summary_txt=tmp_path / "out" / "comparison_summary.txt",
    )


def test_build_comparison_summary_has_required_top_level_keys(tmp_path):
    s = _build_minimal_summary(tmp_path)
    assert set(s) == {
        "schema_version", "neat_version", "generated_at",
        "inputs", "happy", "counts", "metrics", "fn_attribution",
        "warnings", "outputs",
    }


def test_build_comparison_summary_schema_version_matches_constant(tmp_path):
    s = _build_minimal_summary(tmp_path)
    assert s["schema_version"] == REPORT_SCHEMA_VERSION


def test_build_comparison_summary_resolves_paths_absolute(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    s = build_comparison_summary(
        golden_vcf=Path("g.vcf"),
        called_vcf=Path("c.vcf"),
        neat_run_dir=Path("run"),
        simulation_summary_path=Path("run/simulation_summary.json"),
        happy_output_vcf=Path("out/happy.vcf.gz"),
        happy_output_prefix=Path("out/happy"),
        counts={"TP": 0, "FN": 0, "FP": 0},
        fn_attribution={},
        fn_with_reasons_vcf=Path("out/FN_with_reasons.vcf"),
        comparison_summary_json=Path("out/comparison_summary.json"),
        comparison_summary_txt=Path("out/comparison_summary.txt"),
    )
    for v in s["inputs"].values():
        assert Path(v).is_absolute()
    for v in s["outputs"].values():
        assert Path(v).is_absolute()


def test_build_comparison_summary_carries_counts_and_metrics(tmp_path):
    s = _build_minimal_summary(tmp_path, counts={"TP": 8, "FN": 2, "FP": 1})
    assert s["counts"] == {"TP": 8, "FN": 2, "FP": 1}
    assert s["metrics"]["precision"] == pytest.approx(8 / 9)
    assert s["metrics"]["recall"] == pytest.approx(0.8)


# ===========================================================================
# write_comparison_summary_json / write_comparison_summary_txt
# ===========================================================================

def test_write_comparison_summary_json_round_trips(tmp_path):
    s = _build_minimal_summary(tmp_path)
    path = write_comparison_summary_json(s, tmp_path / "comparison_summary.json")
    assert path.is_file()
    parsed = json.loads(path.read_text())
    assert parsed["schema_version"] == REPORT_SCHEMA_VERSION
    assert parsed["counts"] == {"TP": 1, "FN": 1, "FP": 1}


def test_write_comparison_summary_json_atomic(tmp_path):
    """The temp file should not linger after success."""
    s = _build_minimal_summary(tmp_path)
    write_comparison_summary_json(s, tmp_path / "report.json")
    assert not (tmp_path / "report.json.tmp").exists()


def test_write_comparison_summary_txt_contains_expected_sections(tmp_path):
    s = _build_minimal_summary(tmp_path)
    path = write_comparison_summary_txt(s, tmp_path / "report.txt")
    content = path.read_text()
    for section in ["NEAT compare-vcfs report", "Inputs", "Classification", "Metrics",
                    "FN attribution", "Outputs"]:
        assert section in content


def test_write_comparison_summary_txt_includes_counts_and_metrics(tmp_path):
    s = _build_minimal_summary(tmp_path, counts={"TP": 90, "FN": 10, "FP": 5})
    txt = write_comparison_summary_txt(s, tmp_path / "report.txt").read_text()
    assert "TP): 90" in txt
    assert "FN): 10" in txt
    assert "FP): 5" in txt
    assert "0.9474" in txt  # precision
    assert "0.9000" in txt  # recall


def test_render_summary_txt_handles_no_fns(tmp_path):
    s = _build_minimal_summary(tmp_path, counts={"TP": 10, "FN": 0, "FP": 0}, fn_attr={})
    txt = render_summary_txt(s)
    assert "(no false negatives)" in txt


def test_render_summary_txt_renders_NA_for_undefined_metrics(tmp_path):
    s = _build_minimal_summary(tmp_path, counts={"TP": 0, "FN": 0, "FP": 0}, fn_attr={})
    txt = render_summary_txt(s)
    assert "N/A" in txt


# ===========================================================================
# write_fn_with_reasons
# ===========================================================================

_HAPPY_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
"""


def _make_source_vcf(tmp_path: Path, lines: list[str]) -> Path:
    raw = tmp_path / "happy.vcf"
    raw.write_text(_HAPPY_HEADER + "\n".join(lines) + ("\n" if lines else ""))
    return Path(pysam.tabix_index(str(raw), preset="vcf", force=True))


def test_write_fn_with_reasons_adds_neat_reason_info(tmp_path):
    source = _make_source_vcf(tmp_path, [
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:FN\t.:.",
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:FN\t.:.",
    ])
    from neat.compare_vcfs.happy import parse_happy_output
    fn_records = parse_happy_output(source)["FN"]
    fn_reasons = [
        (fn_records[0], [REASON_UNKNOWN]),
        (fn_records[1], [REASON_OUTSIDE_MUTATION_BED, REASON_OUTSIDE_TARGET_BED]),
    ]
    out = write_fn_with_reasons(source, fn_reasons, tmp_path / "FN_with_reasons.vcf")
    assert out.is_file()

    # Re-open and verify the NEAT_REASON tag is present on each record
    with pysam.VariantFile(str(out)) as vf:
        assert "NEAT_REASON" in vf.header.info
        records = list(vf)
    assert len(records) == 2
    assert records[0].info["NEAT_REASON"] == (REASON_UNKNOWN,)
    assert tuple(records[1].info["NEAT_REASON"]) == (
        REASON_OUTSIDE_MUTATION_BED, REASON_OUTSIDE_TARGET_BED,
    )


def test_write_fn_with_reasons_writes_only_fn_records(tmp_path):
    """Source has TP/FN/FP; output should only contain the FN we passed in."""
    source = _make_source_vcf(tmp_path, [
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:TP\t1|0:TP",
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:FN\t.:.",
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:BD\t.:.\t1|0:FP",
    ])
    from neat.compare_vcfs.happy import parse_happy_output
    fn_records = parse_happy_output(source)["FN"]
    out = write_fn_with_reasons(
        source,
        [(fn_records[0], [REASON_UNKNOWN])],
        tmp_path / "FN_with_reasons.vcf",
    )
    with pysam.VariantFile(str(out)) as vf:
        records = list(vf)
    assert len(records) == 1
    assert records[0].pos == 200


def test_write_fn_with_reasons_handles_empty_fn_list(tmp_path):
    source = _make_source_vcf(tmp_path, [
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:TP\t1|0:TP",
    ])
    out = write_fn_with_reasons(source, [], tmp_path / "FN_with_reasons.vcf")
    with pysam.VariantFile(str(out)) as vf:
        assert "NEAT_REASON" in vf.header.info
        assert list(vf) == []


# ===========================================================================
# write_fn_attribution_plot
# ===========================================================================

def test_write_fn_attribution_plot_writes_png(tmp_path):
    path = write_fn_attribution_plot(
        {REASON_OUTSIDE_MUTATION_BED: 12, REASON_OUTSIDE_TARGET_BED: 5, REASON_UNKNOWN: 3},
        tmp_path / "fn_attribution.png",
    )
    assert path.is_file()
    # PNG magic bytes
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_write_fn_attribution_plot_handles_no_fns(tmp_path):
    """An empty attribution dict still produces a (placeholder) PNG."""
    path = write_fn_attribution_plot({}, tmp_path / "fn_attribution.png")
    assert path.is_file()
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
