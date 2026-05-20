"""
Tests for neat/compare_vcfs/happy.py — hap.py subprocess invocation and
output-VCF parsing into TP/FN/FP buckets.
"""
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from neat.compare_vcfs.happy import (
    HappyExecutionError,
    HappyParseError,
    _classify,
    parse_happy_output,
    run_happy,
)


# ---------------------------------------------------------------------------
# run_happy
# ---------------------------------------------------------------------------

def test_run_happy_builds_expected_command(tmp_path, monkeypatch):
    """The subprocess invocation must include all forwarded options in order."""
    captured = {}

    def fake_run(cmd, capture_output, text):
        captured["cmd"] = cmd
        # Pretend hap.py wrote its output
        Path(str(tmp_path / "happy") + ".vcf.gz").touch()
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    run_happy(
        happy_bin=Path("/bin/hap.py"),
        golden_vcf=Path("/tmp/g.vcf"),
        called_vcf=Path("/tmp/c.vcf"),
        output_prefix=tmp_path / "happy",
        reference=Path("/tmp/ref.fa"),
        target_bed=Path("/tmp/t.bed"),
    )

    assert captured["cmd"][:5] == [
        "/bin/hap.py", "/tmp/g.vcf", "/tmp/c.vcf", "-o", str(tmp_path / "happy")
    ]
    assert "-r" in captured["cmd"] and "/tmp/ref.fa" in captured["cmd"]
    assert "-T" in captured["cmd"] and "/tmp/t.bed" in captured["cmd"]


def test_run_happy_omits_optional_flags_when_not_given(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        Path(str(tmp_path / "happy") + ".vcf.gz").touch()
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    run_happy(
        happy_bin=Path("/bin/hap.py"),
        golden_vcf=Path("/tmp/g.vcf"),
        called_vcf=Path("/tmp/c.vcf"),
        output_prefix=tmp_path / "happy",
    )
    assert "-r" not in captured["cmd"]
    assert "-T" not in captured["cmd"]


def test_run_happy_returns_output_vcf_path(tmp_path, monkeypatch):
    expected = Path(str(tmp_path / "happy") + ".vcf.gz")

    def fake_run(cmd, **kw):
        expected.touch()
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    got = run_happy(
        happy_bin=Path("/bin/hap.py"),
        golden_vcf=Path("/tmp/g.vcf"),
        called_vcf=Path("/tmp/c.vcf"),
        output_prefix=tmp_path / "happy",
    )
    assert got == expected


def test_run_happy_raises_when_returncode_nonzero(tmp_path, monkeypatch):
    def fake_run(cmd, **kw):
        return SimpleNamespace(returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(HappyExecutionError, match="exited 1"):
        run_happy(
            happy_bin=Path("/bin/hap.py"),
            golden_vcf=Path("/tmp/g.vcf"),
            called_vcf=Path("/tmp/c.vcf"),
            output_prefix=tmp_path / "happy",
        )


def test_run_happy_raises_when_output_vcf_missing(tmp_path, monkeypatch):
    """Even with returncode=0, missing the expected output VCF is a fatal contract violation."""
    def fake_run(cmd, **kw):
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(HappyExecutionError, match="output VCF is missing"):
        run_happy(
            happy_bin=Path("/bin/hap.py"),
            golden_vcf=Path("/tmp/g.vcf"),
            called_vcf=Path("/tmp/c.vcf"),
            output_prefix=tmp_path / "happy",
        )


def test_run_happy_passes_extra_args(tmp_path, monkeypatch):
    captured = {}

    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        Path(str(tmp_path / "happy") + ".vcf.gz").touch()
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    run_happy(
        happy_bin=Path("/bin/hap.py"),
        golden_vcf=Path("/tmp/g.vcf"),
        called_vcf=Path("/tmp/c.vcf"),
        output_prefix=tmp_path / "happy",
        extra_args=["--engine", "vcfeval"],
    )
    assert "--engine" in captured["cmd"] and "vcfeval" in captured["cmd"]


# ---------------------------------------------------------------------------
# _classify
# ---------------------------------------------------------------------------

def test_classify_tp_when_truth_tp():
    assert _classify("TP", ".") == "TP"


def test_classify_tp_when_query_tp():
    assert _classify(".", "TP") == "TP"


def test_classify_tp_when_both_tp():
    assert _classify("TP", "TP") == "TP"


def test_classify_fn():
    assert _classify("FN", ".") == "FN"


def test_classify_fp():
    assert _classify(".", "FP") == "FP"


def test_classify_returns_none_for_nocall():
    assert _classify(".", ".") is None


def test_classify_returns_none_for_unknown_codes():
    assert _classify("N", "N") is None


def test_classify_tp_beats_other_signals():
    """If a record has TP on one side and FN on the other, hap.py considers it matched."""
    assert _classify("TP", "FN") == "TP"
    assert _classify("FN", "TP") == "TP"


# ---------------------------------------------------------------------------
# parse_happy_output — synthetic VCFs
# ---------------------------------------------------------------------------

_HAPPY_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
##FORMAT=<ID=BVT,Number=1,Type=String,Description="Variant Type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
"""


def _make_happy_vcf(tmp_path: Path, lines: list[str]) -> Path:
    """Write a tiny hap.py-shaped VCF and bgzip+index it for pysam."""
    raw = tmp_path / "happy.vcf"
    raw.write_text(_HAPPY_HEADER + "\n".join(lines) + ("\n" if lines else ""))
    bgz_path = pysam.tabix_index(str(raw), preset="vcf", force=True)
    return Path(bgz_path)


def test_parse_happy_output_buckets_tp_fn_fp(tmp_path):
    vcf = _make_happy_vcf(tmp_path, [
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t1|0:TP:SNP\t1|0:TP:SNP",
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t1|0:FN:SNP\t.:.:.",
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t.:.:.\t1|0:FP:SNP",
        "chr1\t400\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t.:.:.\t.:.:.",  # no-call, dropped
    ])
    buckets = parse_happy_output(vcf)
    assert len(buckets["TP"]) == 1
    assert len(buckets["FN"]) == 1
    assert len(buckets["FP"]) == 1
    assert buckets["TP"][0].pos == 100
    assert buckets["FN"][0].pos == 200
    assert buckets["FP"][0].pos == 300


def test_parse_happy_output_empty_file_returns_empty_buckets(tmp_path):
    vcf = _make_happy_vcf(tmp_path, [])
    buckets = parse_happy_output(vcf)
    assert buckets == {"TP": [], "FN": [], "FP": []}


def test_parse_happy_output_raises_on_single_sample(tmp_path):
    """A VCF with only TRUTH is not a hap.py output."""
    header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH
chr1\t100\t.\tA\tT\t.\t.\t.\tGT:BD\t1|0:TP
"""
    raw = tmp_path / "single.vcf"
    raw.write_text(header)
    bgz = Path(pysam.tabix_index(str(raw), preset="vcf", force=True))
    with pytest.raises(HappyParseError, match="TRUTH and QUERY"):
        parse_happy_output(bgz)


def test_parse_happy_output_raises_when_BD_format_absent(tmp_path):
    """A VCF lacking the BD FORMAT is not a hap.py output."""
    header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t1|0\t1|0
"""
    raw = tmp_path / "noBD.vcf"
    raw.write_text(header)
    bgz = Path(pysam.tabix_index(str(raw), preset="vcf", force=True))
    with pytest.raises(HappyParseError, match="BD FORMAT"):
        parse_happy_output(bgz)


def test_parse_happy_output_preserves_record_chrom_pos(tmp_path):
    vcf = _make_happy_vcf(tmp_path, [
        "chr1\t150\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t1|0:FN:SNP\t.:.:.",
        "chr1\t250\t.\tA\tT\t.\t.\t.\tGT:BD:BVT\t1|0:FN:INDEL\t.:.:.",
    ])
    fns = parse_happy_output(vcf)["FN"]
    assert [(r.chrom, r.pos) for r in fns] == [("chr1", 150), ("chr1", 250)]
