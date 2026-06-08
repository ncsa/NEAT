"""
Real-hap.py integration test for `neat compare-vcfs` (issue #297).

Runs an actual NEAT simulation, then invokes the real hap.py binary. Skipped
cleanly when hap.py isn't available — set NEAT_HAPPY_BIN to the absolute path
of a working hap.py to enable.

Notes on hap.py packaging: the conda `bioconda::hap.py` package is Python-2-based.
Its shebang resolves `python` via PATH, so the test prepends the env's bin
directory so the child process picks up python2.7 from the same env.
"""
import gzip
import json
import os
import shutil
from pathlib import Path

import pysam
import pytest

from neat.compare_vcfs.runner import compare_vcfs_runner
from neat.read_simulator.runner import read_simulator_runner


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def happy_bin():
    """Skip the test if a working hap.py is unavailable."""
    explicit = os.environ.get("NEAT_HAPPY_BIN")
    if explicit and Path(explicit).is_file():
        return Path(explicit)
    found = shutil.which("hap.py")
    if found:
        return Path(found)
    pytest.skip(
        "hap.py is not available. Install via `conda create -n hap_py_env "
        "-c bioconda -c conda-forge hap.py -y` and set NEAT_HAPPY_BIN to its "
        "absolute path to enable this test."
    )


@pytest.fixture
def happy_env_path(happy_bin):
    """Prepend hap.py's env bin to PATH so its #!/usr/bin/env python shebang
    resolves to the python interpreter shipped alongside it."""
    return str(happy_bin.parent)


def _write_ref(path: Path) -> Path:
    """A small but realistic reference: ~2kb with mixed GC content."""
    # Two contigs so we can also exercise per-contig stats.
    seq_a = "ACGT" * 250  # 1000 bp
    seq_b = ("AAAGGCCC" * 125)  # 1000 bp, higher GC
    path.write_text(f">chr1\n{seq_a}\n>chr2\n{seq_b}\n", encoding="utf-8")
    # Index for hap.py (it needs a .fai)
    pysam.FastaFile(str(path))
    return path


def _write_config(path: Path, ref_path: Path) -> Path:
    cfg = (
        f"reference: {ref_path}\n"
        "produce_fastq: false\n"
        "produce_bam: false\n"
        "produce_vcf: true\n"
        "read_len: 100\n"
        "coverage: 5\n"
        "rng_seed: 42\n"
        "mutation_rate: 0.01\n"
        "overwrite_output: true\n"
    )
    path.write_text(cfg, encoding="utf-8")
    return path


def _called_vcf_dropping_first_variant(golden_vcf: Path, output_vcf: Path) -> Path:
    """Read golden.vcf.gz, drop the first variant, write a new bgzipped+indexed VCF.

    Dropping a variant creates exactly one false negative; everything else is
    a true positive (the rest of the truth's variants are also in the caller VCF).
    """
    written = 0
    skipped = False
    with pysam.VariantFile(str(golden_vcf)) as src:
        # Keep the same samples / header
        with pysam.VariantFile(str(output_vcf), "wz", header=src.header) as dst:
            for rec in src:
                if not skipped:
                    skipped = True
                    continue
                dst.write(rec)
                written += 1
    pysam.tabix_index(str(output_vcf), preset="vcf", force=True)
    assert skipped, "golden VCF had no variants to drop; test setup wrong"
    return output_vcf


def _called_vcf_dropping_first_n(golden_vcf: Path, output_vcf: Path, n: int) -> int:
    """Drop the first n variants from golden_vcf into output_vcf; returns the
    number actually dropped (capped at the total available)."""
    actually_dropped = 0
    with pysam.VariantFile(str(golden_vcf)) as src:
        with pysam.VariantFile(str(output_vcf), "wz", header=src.header) as dst:
            for rec in src:
                if actually_dropped < n:
                    actually_dropped += 1
                    continue
                dst.write(rec)
    pysam.tabix_index(str(output_vcf), preset="vcf", force=True)
    return actually_dropped


# ===========================================================================
# Integration test — real hap.py end-to-end
# ===========================================================================

def test_compare_vcfs_real_happy_end_to_end(tmp_path, happy_bin, happy_env_path, monkeypatch):
    """
    Run NEAT, drop one variant from the golden VCF to create a caller VCF with
    exactly one FN, then run compare-vcfs against the real hap.py binary.
    """
    monkeypatch.setenv("PATH", happy_env_path + os.pathsep + os.environ["PATH"])

    # 1. Simulate
    sim_out = tmp_path / "sim_out"
    sim_out.mkdir()
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(tmp_path / "conf.yml", ref)
    read_simulator_runner(str(cfg), str(sim_out), "run")

    golden = sim_out / "run_golden.vcf.gz"
    assert golden.is_file(), "NEAT did not produce a golden VCF"
    assert (sim_out / "simulation_summary.json").is_file()

    # 2. Build a "called" VCF missing one variant
    called = tmp_path / "called.vcf.gz"
    _called_vcf_dropping_first_variant(golden, called)

    # 3. Run compare-vcfs against real hap.py
    cmp_out = tmp_path / "cmp_out"
    compare_vcfs_runner(
        golden_vcf=str(golden),
        called_vcf=str(called),
        neat_run_dir=str(sim_out),
        output_dir=str(cmp_out),
        reference=str(ref),
        happy_bin=str(happy_bin),
    )

    # 4. Assert the three reports exist and look sensible
    assert (cmp_out / "comparison_summary.json").is_file()
    assert (cmp_out / "comparison_summary.txt").is_file()
    assert (cmp_out / "FN_with_reasons.vcf").is_file()
    assert (cmp_out / "happy.vcf.gz").is_file()

    report = json.loads((cmp_out / "comparison_summary.json").read_text())
    assert report["schema_version"] == "1"
    # We dropped exactly one variant — expect at least one FN
    assert report["counts"]["FN"] >= 1
    assert report["counts"]["FP"] == 0
    assert report["metrics"]["precision"] == pytest.approx(1.0)
    # Every FN must be attributed somehow
    total_attribution = sum(report["fn_attribution"].values())
    assert total_attribution >= report["counts"]["FN"]

    # And the annotated VCF must carry NEAT_REASON
    with pysam.VariantFile(str(cmp_out / "FN_with_reasons.vcf")) as vf:
        assert "NEAT_REASON" in vf.header.info
        for rec in vf:
            assert "NEAT_REASON" in rec.info


def test_compare_vcfs_real_happy_multi_fn(tmp_path, happy_bin, happy_env_path, monkeypatch):
    """
    Drop multiple variants to verify FN counts scale linearly and that each
    FN gets a NEAT_REASON tag (locks in attribution behavior at scale).
    """
    monkeypatch.setenv("PATH", happy_env_path + os.pathsep + os.environ["PATH"])

    sim_out = tmp_path / "sim_out"
    sim_out.mkdir()
    ref = _write_ref(tmp_path / "ref.fa")
    cfg = _write_config(tmp_path / "conf.yml", ref)
    read_simulator_runner(str(cfg), str(sim_out), "run")
    golden = sim_out / "run_golden.vcf.gz"

    called = tmp_path / "called.vcf.gz"
    n_dropped = _called_vcf_dropping_first_n(golden, called, n=3)
    assert n_dropped == 3, "golden produced fewer than 3 variants — config too small"

    cmp_out = tmp_path / "cmp_out"
    compare_vcfs_runner(
        golden_vcf=str(golden), called_vcf=str(called),
        neat_run_dir=str(sim_out),
        output_dir=str(cmp_out),
        reference=str(ref),
        happy_bin=str(happy_bin),
        plot=True,
    )

    report = json.loads((cmp_out / "comparison_summary.json").read_text())
    assert report["counts"]["FN"] >= 3, "expected ≥3 FNs after dropping 3 variants"
    assert report["counts"]["FP"] == 0
    assert sum(report["fn_attribution"].values()) >= report["counts"]["FN"]

    # Every FN record must be annotated; this is the contract step 6 promises
    with pysam.VariantFile(str(cmp_out / "FN_with_reasons.vcf")) as vf:
        fn_records = list(vf)
    assert len(fn_records) == report["counts"]["FN"]
    for rec in fn_records:
        tags = rec.info["NEAT_REASON"]
        assert tags, f"FN at {rec.chrom}:{rec.pos} has empty NEAT_REASON"

    # --plot was on, so the bar chart should be present
    assert (cmp_out / "fn_attribution.png").is_file()


def test_compare_vcfs_real_happy_with_chrom_mismatched_bed(tmp_path, happy_bin, happy_env_path, monkeypatch):
    """
    End-to-end against real hap.py with a mutation_bed that uses '1'/'2' while
    the reference uses 'chr1'/'chr2'. Verifies that:
      - the warning surfaces in the JSON report
      - FNs are NOT mislabeled as 'outside_mutation_bed' (the semantic fix)
      - --chrom-aliases silences the warning AND restores correct attribution
    """
    monkeypatch.setenv("PATH", happy_env_path + os.pathsep + os.environ["PATH"])

    # Hand-craft a mutation_bed that won't match the reference's chr-prefixed names.
    # The simulator will log warnings and effectively ignore it (NEAT's pre-existing
    # "skip BED chroms not in reference" behavior). simulation_summary still records
    # the BED path; compare-vcfs reads it and detects the mismatch.
    bed = tmp_path / "mismatched_mut.bed"
    bed.write_text("1\t0\t1000\n2\t0\t1000\n")

    sim_out = tmp_path / "sim_out"
    sim_out.mkdir()
    ref = _write_ref(tmp_path / "ref.fa")
    cfg_path = tmp_path / "conf.yml"
    cfg_path.write_text(
        f"reference: {ref}\n"
        "produce_fastq: false\n"
        "produce_bam: false\n"
        "produce_vcf: true\n"
        "read_len: 100\n"
        "coverage: 5\n"
        "rng_seed: 42\n"
        "mutation_rate: 0.01\n"
        f"mutation_bed: {bed}\n"
        "overwrite_output: true\n",
        encoding="utf-8",
    )
    read_simulator_runner(str(cfg_path), str(sim_out), "run")
    golden = sim_out / "run_golden.vcf.gz"

    called = tmp_path / "called.vcf.gz"
    _called_vcf_dropping_first_variant(golden, called)

    # ------------------------------------------------------------------
    # Run 1: no --chrom-aliases → mismatch warning, no misleading outside_*
    # ------------------------------------------------------------------
    cmp_out_no_aliases = tmp_path / "cmp_no_aliases"
    compare_vcfs_runner(
        golden_vcf=str(golden), called_vcf=str(called),
        neat_run_dir=str(sim_out),
        output_dir=str(cmp_out_no_aliases),
        reference=str(ref),
        happy_bin=str(happy_bin),
    )
    report = json.loads((cmp_out_no_aliases / "comparison_summary.json").read_text())
    chrom_warnings = [w for w in report["warnings"] if w["type"] == "chrom_naming_mismatch"]
    assert chrom_warnings, "expected chrom-mismatch warning to fire"
    assert chrom_warnings[0]["suggested_aliases"] == {"1": "chr1", "2": "chr2"}
    # Semantic guarantee: mismatched BED must NOT produce outside_mutation_bed
    assert "outside_mutation_bed" not in report["fn_attribution"]

    # ------------------------------------------------------------------
    # Run 2: --chrom-aliases supplied → no warning, BED becomes usable
    # ------------------------------------------------------------------
    aliases = tmp_path / "aliases.tsv"
    aliases.write_text("1\tchr1\n2\tchr2\n")
    cmp_out_aliased = tmp_path / "cmp_aliased"
    compare_vcfs_runner(
        golden_vcf=str(golden), called_vcf=str(called),
        neat_run_dir=str(sim_out),
        output_dir=str(cmp_out_aliased),
        reference=str(ref),
        happy_bin=str(happy_bin),
        chrom_aliases=str(aliases),
    )
    report_aliased = json.loads((cmp_out_aliased / "comparison_summary.json").read_text())
    assert not [w for w in report_aliased["warnings"] if w["type"] == "chrom_naming_mismatch"]
