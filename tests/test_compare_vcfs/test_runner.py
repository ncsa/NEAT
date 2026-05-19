"""
Tests for neat/compare_vcfs/runner.py — the scaffold for `neat compare-vcfs`.

Covers input validation, hap.py discovery, simulation_summary.json loading,
and the scaffold's NotImplementedError contract.
"""
import json
import os
import stat
from pathlib import Path

import pytest

from neat.compare_vcfs.runner import (
    HappyNotFoundError,
    SimulationSummaryError,
    compare_vcfs_runner,
    discover_happy,
    load_simulation_summary,
)
from neat.read_simulator.utils.simulation_summary import SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Synthetic-file builders
# ---------------------------------------------------------------------------

def _touch(path: Path, content: str = "stub") -> Path:
    """Write non-empty content; validate_input_path sys.exits on empty files."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _make_summary_dir(tmp_path: Path, schema_version: str = SCHEMA_VERSION,
                     extra: dict | None = None) -> Path:
    """Create a NEAT-run-dir with a valid (or version-bumped) simulation_summary.json."""
    run_dir = tmp_path / "neat_run"
    run_dir.mkdir()
    summary = {
        "schema_version": schema_version,
        "neat_version": "4.4.4",
        "run": {},
        "config": {},
        "outputs": {},
        "delivered": {
            "total_variants": 0,
            "contigs_simulated": ["chr1"],
        },
    }
    if extra:
        summary.update(extra)
    (run_dir / "simulation_summary.json").write_text(json.dumps(summary))
    return run_dir


def _make_executable(path: Path) -> Path:
    path.touch()
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


# ===========================================================================
# discover_happy
# ===========================================================================

def test_discover_happy_returns_path_from_PATH(tmp_path, monkeypatch):
    fake = _make_executable(tmp_path / "hap.py")
    monkeypatch.setenv("PATH", str(tmp_path), prepend=os.pathsep)
    assert discover_happy(None) == fake.resolve()


def test_discover_happy_raises_when_not_on_PATH(monkeypatch):
    monkeypatch.setenv("PATH", "")
    with pytest.raises(HappyNotFoundError, match="not found on \\$PATH"):
        discover_happy(None)


def test_discover_happy_returns_explicit_path_when_file_exists(tmp_path):
    fake = _make_executable(tmp_path / "my_happy")
    assert discover_happy(str(fake)) == fake.resolve()


def test_discover_happy_raises_when_explicit_path_missing(tmp_path):
    with pytest.raises(HappyNotFoundError, match="path does not exist"):
        discover_happy(str(tmp_path / "nope"))


def test_discover_happy_explicit_path_wins_over_PATH(tmp_path, monkeypatch):
    """An explicit --happy-bin should be used even if hap.py is also on PATH."""
    on_path = _make_executable(tmp_path / "hap.py")
    explicit = _make_executable(tmp_path / "alt_hap.py")
    monkeypatch.setenv("PATH", str(tmp_path), prepend=os.pathsep)
    assert discover_happy(str(explicit)) == explicit.resolve()
    assert discover_happy(str(explicit)) != on_path.resolve()


# ===========================================================================
# load_simulation_summary
# ===========================================================================

def test_load_simulation_summary_returns_dict_for_valid_file(tmp_path):
    run_dir = _make_summary_dir(tmp_path)
    data = load_simulation_summary(run_dir)
    assert data["schema_version"] == SCHEMA_VERSION
    assert data["neat_version"] == "4.4.4"


def test_load_simulation_summary_raises_when_missing(tmp_path):
    run_dir = tmp_path / "empty_run"
    run_dir.mkdir()
    with pytest.raises(SimulationSummaryError, match="not found"):
        load_simulation_summary(run_dir)


def test_load_simulation_summary_raises_on_malformed_json(tmp_path):
    run_dir = tmp_path / "broken"
    run_dir.mkdir()
    (run_dir / "simulation_summary.json").write_text("{not valid json")
    with pytest.raises(SimulationSummaryError, match="not valid JSON"):
        load_simulation_summary(run_dir)


def test_load_simulation_summary_raises_on_version_mismatch(tmp_path):
    run_dir = _make_summary_dir(tmp_path, schema_version="99")
    with pytest.raises(SimulationSummaryError, match="schema_version"):
        load_simulation_summary(run_dir)


def test_load_simulation_summary_raises_on_missing_schema_version(tmp_path):
    """A summary without schema_version is treated as a version mismatch (None ≠ '1')."""
    run_dir = tmp_path / "no_version"
    run_dir.mkdir()
    (run_dir / "simulation_summary.json").write_text(json.dumps({"foo": "bar"}))
    with pytest.raises(SimulationSummaryError, match="schema_version"):
        load_simulation_summary(run_dir)


# ===========================================================================
# compare_vcfs_runner — input validation
# ===========================================================================

def test_runner_exits_when_golden_missing(tmp_path):
    """Missing truth VCF: validate_input_path sys.exit(5)."""
    run_dir = _make_summary_dir(tmp_path)
    called = _touch(tmp_path / "called.vcf")
    with pytest.raises(SystemExit) as excinfo:
        compare_vcfs_runner(
            golden_vcf=str(tmp_path / "nope.vcf"),
            called_vcf=str(called),
            neat_run_dir=str(run_dir),
            output_dir=str(tmp_path / "out"),
        )
    assert excinfo.value.code == 5


def test_runner_exits_when_called_missing(tmp_path):
    run_dir = _make_summary_dir(tmp_path)
    golden = _touch(tmp_path / "golden.vcf")
    with pytest.raises(SystemExit) as excinfo:
        compare_vcfs_runner(
            golden_vcf=str(golden),
            called_vcf=str(tmp_path / "nope.vcf"),
            neat_run_dir=str(run_dir),
            output_dir=str(tmp_path / "out"),
        )
    assert excinfo.value.code == 5


def test_runner_raises_when_run_dir_missing(tmp_path):
    golden = _touch(tmp_path / "golden.vcf")
    called = _touch(tmp_path / "called.vcf")
    with pytest.raises(FileNotFoundError, match="--neat-run-dir"):
        compare_vcfs_runner(
            golden_vcf=str(golden), called_vcf=str(called),
            neat_run_dir=str(tmp_path / "nope_dir"),
            output_dir=str(tmp_path / "out"),
        )


def test_runner_raises_when_happy_missing(tmp_path, monkeypatch):
    run_dir = _make_summary_dir(tmp_path)
    golden = _touch(tmp_path / "golden.vcf")
    called = _touch(tmp_path / "called.vcf")
    monkeypatch.setenv("PATH", "")  # nothing on PATH
    with pytest.raises(HappyNotFoundError):
        compare_vcfs_runner(
            golden_vcf=str(golden), called_vcf=str(called),
            neat_run_dir=str(run_dir),
            output_dir=str(tmp_path / "out"),
        )


def test_runner_completes_end_to_end_writing_all_reports(tmp_path, monkeypatch):
    """End-to-end happy path: mocked hap.py writes a tiny output VCF; runner
    produces all three reports and does not raise."""
    run_dir = _make_summary_dir(tmp_path)
    golden = _touch(tmp_path / "golden.vcf")
    called = _touch(tmp_path / "called.vcf")
    _make_executable(tmp_path / "hap.py")
    monkeypatch.setenv("PATH", str(tmp_path), prepend=os.pathsep)
    out_dir = tmp_path / "out"

    from neat.compare_vcfs import runner as runner_mod

    happy_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:TP\t1|0:TP
chr1\t200\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:FN\t.:.
chr1\t300\t.\tA\tT\t.\t.\t.\tGT:BD\t.:.\t1|0:FP
"""

    def fake_run_happy(happy_bin, golden_vcf, called_vcf, output_prefix, **kw):
        import pysam
        raw = Path(str(output_prefix) + ".vcf")
        raw.write_text(happy_header)
        bgz = Path(pysam.tabix_index(str(raw), preset="vcf", force=True))
        return bgz

    monkeypatch.setattr(runner_mod, "run_happy", fake_run_happy)

    compare_vcfs_runner(
        golden_vcf=str(golden), called_vcf=str(called),
        neat_run_dir=str(run_dir),
        output_dir=str(out_dir),
    )

    assert (out_dir / "comparison_summary.json").is_file()
    assert (out_dir / "comparison_summary.txt").is_file()
    assert (out_dir / "FN_with_reasons.vcf").is_file()

    import json as _json
    report = _json.loads((out_dir / "comparison_summary.json").read_text())
    assert report["counts"] == {"TP": 1, "FN": 1, "FP": 1}
    assert report["metrics"]["precision"] == pytest.approx(0.5)
    assert report["metrics"]["recall"] == pytest.approx(0.5)
    assert report["metrics"]["f1"] == pytest.approx(0.5)
    # chr1 IS simulated and no beds configured → FN attributed as 'unknown'
    assert report["fn_attribution"] == {"unknown": 1}


def test_runner_validates_optional_reference_when_provided(tmp_path):
    """A passed --reference that doesn't exist should fail fast via validate_input_path."""
    run_dir = _make_summary_dir(tmp_path)
    golden = _touch(tmp_path / "golden.vcf")
    called = _touch(tmp_path / "called.vcf")
    with pytest.raises(SystemExit) as excinfo:
        compare_vcfs_runner(
            golden_vcf=str(golden), called_vcf=str(called),
            neat_run_dir=str(run_dir),
            output_dir=str(tmp_path / "out"),
            reference=str(tmp_path / "missing_ref.fa"),
        )
    assert excinfo.value.code == 5


def test_runner_validates_optional_target_bed_when_provided(tmp_path):
    run_dir = _make_summary_dir(tmp_path)
    golden = _touch(tmp_path / "golden.vcf")
    called = _touch(tmp_path / "called.vcf")
    with pytest.raises(SystemExit) as excinfo:
        compare_vcfs_runner(
            golden_vcf=str(golden), called_vcf=str(called),
            neat_run_dir=str(run_dir),
            output_dir=str(tmp_path / "out"),
            target_bed=str(tmp_path / "missing_target.bed"),
        )
    assert excinfo.value.code == 5
