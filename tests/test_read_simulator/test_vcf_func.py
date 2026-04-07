"""
Regression tests for neat/read_simulator/utils/vcf_func.py

Focus: parse_input_vcf, specifically the WP genotype condition fix.

Bug fixed on branch fix/vcf-wp-genotype-condition:
    BEFORE: "WP" in [x.split('=') for x in record[7].split(';')]
            — always False because "WP" (str) can never be a member of a list
              of lists.
    AFTER:  "WP" in [x.split('=')[0] for x in record[7].split(';')]
            — correctly extracts INFO keys and tests membership.

There are two call sites for the WP condition, exercising different code paths:
    Path A (line 164): has_format=True, FORMAT column lacks a GT field
    Path B (line 185): has_format=False (no FORMAT column at all)
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO

from neat.read_simulator.utils.vcf_func import parse_input_vcf
from neat.read_simulator.utils.options import Options
from neat.variants import ContigVariants

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_REF_SEQ = "ACGT" * 100   # 400 bp; pos n (0-based) = "ACGT"[n % 4]
# 0-based position 1 → 'C'; VCF POS (1-based) = 2
_SNV_POS_VCF = 2   # 1-based
_SNV_REF     = "C"
_SNV_ALT     = "T"


def _write_ref(tmp_path: Path) -> Path:
    ref = tmp_path / "ref.fa"
    ref.write_text(f">chr1\n{_REF_SEQ}\n", encoding="utf-8")
    return ref


def _make_ref_index(tmp_path: Path):
    ref = _write_ref(tmp_path)
    return SeqIO.index(str(ref), "fasta")


def _make_input_dict():
    return {"chr1": ContigVariants()}


def _make_options(seed: int = 0) -> Options:
    opts = Options(rng_seed=seed)
    opts.ploidy = 2
    return opts


def _write_vcf(path: Path, header_cols: str, info: str, format_col: str = None,
               sample_col: str = None) -> Path:
    """
    Write a minimal single-variant VCF.

    header_cols: the #CHROM header line columns after FILTER (e.g. '' or 'FORMAT\tSAMPLE')
    info:        the INFO field value
    format_col:  FORMAT field value (None → no FORMAT column)
    sample_col:  SAMPLE field value (None → no sample column)
    """
    extra_header_cols = f"\t{header_cols}" if header_cols else ""
    extra_data_cols = ""
    if format_col is not None:
        extra_data_cols += f"\t{format_col}"
    if sample_col is not None:
        extra_data_cols += f"\t{sample_col}"

    path.write_text(
        "##fileformat=VCFv4.1\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{extra_header_cols}\n"
        f"chr1\t{_SNV_POS_VCF}\t.\t{_SNV_REF}\t{_SNV_ALT}\t42\tPASS\t{info}{extra_data_cols}\n",
        encoding="utf-8",
    )
    return path


def _get_variant(input_dict):
    """Return the first (and expected only) variant from the chr1 ContigVariants."""
    cv = input_dict["chr1"]
    assert cv.variant_locations, "No variants were added to ContigVariants"
    loc = cv.variant_locations[0]
    return cv.contig_variants[loc][0]


# ===========================================================================
# Path B — no FORMAT column, WP in INFO (the cleanest WP path)
# ===========================================================================

def test_no_format_wp_in_info_genotype_matches_wp_value(tmp_path):
    """
    Path B (line 185): no FORMAT column, WP=0|1 in INFO.

    Before the fix the WP condition always evaluated to False and genotype was
    generated randomly. After the fix the genotype must equal [0, 1].
    """
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info="WP=0|1")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


def test_no_format_wp_slash_notation_converted(tmp_path):
    """WP genotype written with '/' separator is correctly converted to '|'."""
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info="WP=0/1")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


def test_no_format_wp_homozygous_alt(tmp_path):
    """WP=1|1 (homozygous alt) is parsed correctly."""
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info="WP=1|1")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [1, 1])


def test_no_format_wp_among_other_info_fields(tmp_path):
    """WP is found even when other INFO fields precede it."""
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info="DP=42;WP=0|1;DB")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


def test_no_format_no_wp_genotype_is_randomly_generated(tmp_path):
    """Without WP and without FORMAT, genotype is assigned randomly (not from WP)."""
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info="DP=42")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options(seed=0))

    # We can't predict the random genotype, but the variant must exist and
    # genotype must be a valid numpy array of length == ploidy (2).
    variant = _get_variant(input_dict)
    assert variant.genotype is not None
    assert len(variant.genotype) == 2


def test_no_format_dot_info_genotype_randomly_generated(tmp_path):
    """INFO field '.' (missing) does not trigger WP path; genotype is random."""
    vcf = _write_vcf(tmp_path / "in.vcf", header_cols="", info=".")
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options(seed=0))

    variant = _get_variant(input_dict)
    assert len(variant.genotype) == 2


# ===========================================================================
# Path A — has FORMAT column, FORMAT lacks GT, WP in INFO
# ===========================================================================

def test_has_format_no_gt_wp_in_info_genotype_from_wp(tmp_path):
    """
    Path A (line 164): FORMAT column present but without GT, WP=0|1 in INFO.

    Before the fix the WP condition was never reached. After the fix it is
    reached and genotype is set from the WP field.

    Note: there is a secondary variable-shadowing bug on line 176 that causes
    normal_sample_field to be constructed incorrectly (record[9] indexes into
    the string segment, not the original record). The genotype itself is still
    set correctly — this test validates that part of the fix only.
    """
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info="WP=0|1",
        format_col="DP",       # FORMAT has DP, not GT
        sample_col="40",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


def test_has_format_no_gt_no_wp_genotype_randomly_generated(tmp_path):
    """
    Path A else branch (line 178): FORMAT present but no GT and no WP →
    genotype is generated randomly. Validates that the WP fix does not break
    the fallback random-genotype path.
    """
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info="DP=50",
        format_col="DP",
        sample_col="50",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options(seed=0))

    variant = _get_variant(input_dict)
    assert variant.genotype is not None
    assert len(variant.genotype) == 2


# ===========================================================================
# Standard GT path — unaffected by WP fix
# ===========================================================================

def test_has_format_with_gt_genotype_from_sample_column(tmp_path):
    """
    FORMAT column contains GT — genotype is read from the SAMPLE column.
    This is the most common VCF path; the WP fix must not affect it.
    """
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info=".",
        format_col="GT",
        sample_col="0|1",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


def test_has_format_gt_homozygous_ref(tmp_path):
    """GT=0|0 (homozygous ref) is parsed correctly via the standard GT path."""
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info=".",
        format_col="GT",
        sample_col="0|0",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 0])


def test_has_format_gt_slash_notation(tmp_path):
    """GT=0/1 (unphased, slash separator) is handled correctly."""
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info=".",
        format_col="GT",
        sample_col="0/1",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    np.testing.assert_array_equal(variant.genotype, [0, 1])


# ===========================================================================
# General parse_input_vcf correctness
# ===========================================================================

def test_variant_not_in_reference_skipped(tmp_path):
    """A variant whose chromosome is absent from the reference is silently skipped."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chrX\t2\t.\tC\tT\t42\tPASS\t.\n",
        encoding="utf-8",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    assert input_dict["chr1"].variant_locations == []


def test_ref_mismatch_skipped(tmp_path):
    """A variant whose REF field doesn't match the reference is skipped."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{_SNV_POS_VCF}\t.\tG\tT\t42\tPASS\t.\n",  # REF='G' but actual is 'C'
        encoding="utf-8",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    assert input_dict["chr1"].variant_locations == []


def test_missing_qual_gets_default(tmp_path):
    """A '.' QUAL is replaced by the default '42'."""
    vcf = _write_vcf(
        tmp_path / "in.vcf",
        header_cols="FORMAT\tSAMPLE",
        info=".",
        format_col="GT",
        sample_col="0|1",
    )
    # Rewrite with QUAL='.'
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"chr1\t{_SNV_POS_VCF}\t.\t{_SNV_REF}\t{_SNV_ALT}\t.\tPASS\t.\tGT\t0|1\n",
        encoding="utf-8",
    )
    ref = _make_ref_index(tmp_path)
    input_dict = _make_input_dict()

    parse_input_vcf(input_dict, vcf, 2, ref, _make_options())

    variant = _get_variant(input_dict)
    assert variant.qual_score == "42"