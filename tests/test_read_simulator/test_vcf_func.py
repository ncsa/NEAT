"""
Tests for neat/read_simulator/utils/vcf_func.py
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO

from neat.read_simulator.utils.options import Options
from neat.read_simulator.utils.vcf_func import (
    parse_input_vcf,
    retrieve_genotype,
    variant_genotype,
)
from neat.variants import SingleNucleotideVariant
from neat.variants.contig_variants import ContigVariants
from neat.variants.deletion import Deletion
from neat.variants.insertion import Insertion
from neat.variants.unknown_variant import UnknownVariant


# Shared fixtures and helpers

# Reference sequence: chr1 = ACGTACGTACGTACGTACGT (20 bp)
#                     chr2 = TTGGTTGGTTGG          (12 bp)
_REF_TEXT = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTGGTTGGTTGG\n"

# VCF column separator
_TAB = "\t"


@pytest.fixture()
def ref_fasta(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(_REF_TEXT, encoding="utf-8")
    return SeqIO.index(str(fa), "fasta")


@pytest.fixture()
def empty_input_dict():
    return {"chr1": ContigVariants(), "chr2": ContigVariants()}


@pytest.fixture()
def opts():
    o = Options(rng_seed=42)
    o.ploidy = 2
    return o


def _write_vcf(tmp_path: Path, name: str, lines: list[str]) -> Path:
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _vcf_header_no_format():
    return [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]


def _vcf_header_with_format(sample="SAMPLE1"):
    return [
        "##fileformat=VCFv4.2",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]


# retrieve_genotype

def _make_vcf_record(format_field, sample_field, info="."):
    """Build a minimal 10-column VCF record list."""
    return ["chr1", 0, ".", "A", "G", "30", "PASS", info, format_field, sample_field]


def test_retrieve_genotype_phased():
    record = _make_vcf_record("GT", "0|1")
    gt = retrieve_genotype(record)
    np.testing.assert_array_equal(gt, [0, 1])


def test_retrieve_genotype_unphased_converted():
    """/ separator is normalised to | before splitting."""
    record = _make_vcf_record("GT", "0/1")
    gt = retrieve_genotype(record)
    np.testing.assert_array_equal(gt, [0, 1])


def test_retrieve_genotype_homozygous_ref():
    record = _make_vcf_record("GT", "0|0")
    np.testing.assert_array_equal(retrieve_genotype(record), [0, 0])


def test_retrieve_genotype_homozygous_alt():
    record = _make_vcf_record("GT", "1|1")
    np.testing.assert_array_equal(retrieve_genotype(record), [1, 1])


def test_retrieve_genotype_gt_among_other_fields():
    """GT may appear after other FORMAT fields."""
    record = _make_vcf_record("DP:GT:GQ", "30:0|1:99")
    np.testing.assert_array_equal(retrieve_genotype(record), [0, 1])


def test_retrieve_genotype_cancer_uses_column_10():
    """is_cancer=True should read from column index 10, not 9."""
    record = _make_vcf_record("GT", "0|0") + ["1|1"]  # 11 columns; col 10 = "1|1"
    gt = retrieve_genotype(record, is_cancer=True)
    np.testing.assert_array_equal(gt, [1, 1])


# variant_genotype

def test_variant_genotype_no_match():
    gt = variant_genotype(2, np.array([0, 0]), 1)
    np.testing.assert_array_equal(gt, [0, 0])


def test_variant_genotype_het():
    gt = variant_genotype(2, np.array([0, 1]), 1)
    np.testing.assert_array_equal(gt, [0, 1])


def test_variant_genotype_homozygous_alt():
    gt = variant_genotype(2, np.array([1, 1]), 1)
    np.testing.assert_array_equal(gt, [1, 1])


def test_variant_genotype_second_alt():
    """which_alt=2 matches ploids where full_genotype==2."""
    gt = variant_genotype(3, np.array([0, 1, 2]), 2)
    np.testing.assert_array_equal(gt, [0, 0, 1])


def test_variant_genotype_returns_correct_ploidy_length():
    gt = variant_genotype(4, np.array([1, 0, 1, 0]), 1)
    assert len(gt) == 4
    # Ploids where full_genotype == which_alt (1) get 1; others get 0
    np.testing.assert_array_equal(gt, [1, 0, 1, 0])


# parse_input_vcf — variant type classification

def test_parse_snv(tmp_path, ref_fasta, empty_input_dict, opts):
    """REF and ALT both length 1 → SingleNucleotideVariant."""
    vcf = _write_vcf(tmp_path, "snv.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",   # pos 1 (VCF) = pos 0 (0-based): ref=A
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants.get(0, [])
    assert len(variants) == 1
    assert isinstance(variants[0], SingleNucleotideVariant)
    assert variants[0].alt == "G"


def test_parse_deletion(tmp_path, ref_fasta, empty_input_dict, opts):
    """len(ref) > len(alt) and ref.startswith(alt) → Deletion."""
    # pos 5 (VCF) = pos 4 (0-based): reference is ACGT
    vcf = _write_vcf(tmp_path, "del.vcf", _vcf_header_no_format() + [
        "chr1\t5\t.\tACGT\tA\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants.get(4, [])
    assert len(variants) == 1
    assert isinstance(variants[0], Deletion)


def test_parse_insertion(tmp_path, ref_fasta, empty_input_dict, opts):
    """len(alt) > len(ref) and alt.startswith(ref) → Insertion."""
    # pos 9 (VCF) = pos 8 (0-based): reference is A
    vcf = _write_vcf(tmp_path, "ins.vcf", _vcf_header_no_format() + [
        "chr1\t9\t.\tA\tACGTACGT\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants.get(8, [])
    assert len(variants) == 1
    assert isinstance(variants[0], Insertion)


def test_parse_unknown_variant(tmp_path, ref_fasta, empty_input_dict, opts):
    """MNV (multi-nucleotide, same length > 1) falls through to UnknownVariant."""
    # pos 1 (VCF) = pos 0 (0-based): ref=AC, alt=TG (len==2, not 1, so not SNV)
    vcf = _write_vcf(tmp_path, "unk.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tAC\tTG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants.get(0, [])
    assert len(variants) == 1
    assert isinstance(variants[0], UnknownVariant)


# parse_input_vcf — filtering / skipping

def test_chrom_not_in_reference_skipped(tmp_path, ref_fasta, empty_input_dict, opts):
    """Variants on chromosomes absent from the reference are silently skipped."""
    vcf = _write_vcf(tmp_path, "chrom.vcf", _vcf_header_no_format() + [
        "chrX\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    total = sum(len(cv.variant_locations) for cv in empty_input_dict.values())
    assert total == 0


def test_ref_mismatch_skipped(tmp_path, ref_fasta, empty_input_dict, opts):
    """Variant whose REF doesn't match the reference sequence is skipped."""
    # chr1 pos 0 is 'A', but we claim it's 'T'
    vcf = _write_vcf(tmp_path, "mismatch.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tT\tG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert len(empty_input_dict["chr1"].variant_locations) == 0


def test_duplicate_position_same_alt_skipped(tmp_path, ref_fasta, empty_input_dict, opts):
    """A second variant at the same position with the same ALT is a duplicate and skipped,
    regardless of genotype. add_variant deduplicates by (position, type, ALT)."""
    vcf = _write_vcf(tmp_path, "dup_same_alt.vcf", _vcf_header_with_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT\t0|1",
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT\t1|0",  # same ALT 'G', different genotype
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert len(empty_input_dict["chr1"].contig_variants[0]) == 1


def test_duplicate_position_different_alt_both_accepted(tmp_path, ref_fasta, empty_input_dict, opts):
    """Two variants at the same position with different ALTs are not duplicates — both accepted."""
    vcf = _write_vcf(tmp_path, "dup_diff_alt.vcf", _vcf_header_with_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT\t0|1",
        "chr1\t1\t.\tA\tC\t30\tPASS\t.\tGT\t0|1",  # different ALT 'C'
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert len(empty_input_dict["chr1"].contig_variants[0]) == 2


def test_comment_and_header_lines_not_parsed_as_variants(tmp_path, ref_fasta, empty_input_dict, opts):
    """## header lines and #CHROM line are never treated as variant records."""
    vcf = _write_vcf(tmp_path, "headers.vcf", [
        "##fileformat=VCFv4.2",
        "##source=test",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert len(empty_input_dict["chr1"].variant_locations) == 1


# parse_input_vcf — QUAL handling

def test_missing_qual_replaced_with_42(tmp_path, ref_fasta, empty_input_dict, opts):
    """QUAL field '.' is replaced with the default value '42'."""
    vcf = _write_vcf(tmp_path, "qual.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t.\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants[0]
    assert variants[0].qual_score == "42"


# parse_input_vcf — FORMAT / genotype handling

def test_with_format_gt_uses_sample_genotype(tmp_path, ref_fasta, empty_input_dict, opts):
    """FORMAT column with GT field reads genotype from the sample column."""
    vcf = _write_vcf(tmp_path, "gt.vcf", _vcf_header_with_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT\t0|1",
    ])
    sample_cols = parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert sample_cols == {"SAMPLE1": 7}
    variants = empty_input_dict["chr1"].contig_variants[0]
    # genotype should reflect 0|1: only the second ploid carries the variant
    np.testing.assert_array_equal(variants[0].genotype, [0, 1])


def test_without_format_genotype_is_generated(tmp_path, ref_fasta, empty_input_dict, opts):
    """No FORMAT column → genotype is randomly generated (non-None)."""
    vcf = _write_vcf(tmp_path, "nogt.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants[0]
    assert variants[0].genotype is not None
    assert len(variants[0].genotype) == 2


def test_format_without_gt_generates_random_genotype(tmp_path, ref_fasta, empty_input_dict, opts):
    """FORMAT column present but no GT field → random genotype generated."""
    vcf = _write_vcf(tmp_path, "fmtngt.vcf", _vcf_header_with_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tDP\t42",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    variants = empty_input_dict["chr1"].contig_variants[0]
    assert variants[0].genotype is not None
    assert len(variants[0].genotype) == 2


def test_no_format_returns_empty_sample_columns(tmp_path, ref_fasta, empty_input_dict, opts):
    """No FORMAT column → returned sample_columns is empty."""
    vcf = _write_vcf(tmp_path, "nosample.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    result = parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert not result


def test_format_exits_if_no_sample_column(tmp_path, ref_fasta, empty_input_dict, opts):
    """FORMAT present but no sample column after it → sys.exit."""
    vcf = _write_vcf(tmp_path, "fmtnosample.vcf", [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT",
    ])
    with pytest.raises(SystemExit):
        parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)


# parse_input_vcf — multiple ALTs

def test_multiple_alts_each_gets_variant(tmp_path, ref_fasta, empty_input_dict, opts):
    """A comma-separated ALT field produces one variant object per alt allele."""
    vcf = _write_vcf(tmp_path, "multialt.vcf", _vcf_header_with_format() + [
        "chr1\t1\t.\tA\tG,C\t30\tPASS\t.\tGT\t1|2",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    # Two SNVs at position 0
    variants = empty_input_dict["chr1"].contig_variants.get(0, [])
    assert len(variants) == 2
    alts = {v.alt for v in variants}
    assert alts == {"G", "C"}


# parse_input_vcf — multiple contigs and is_input flag

def test_variants_routed_to_correct_contig(tmp_path, ref_fasta, empty_input_dict, opts):
    """Variants on different chromosomes end up in the correct ContigVariants."""
    # chr2 starts with TTGG... so pos 1 (0-based 0) = T
    vcf = _write_vcf(tmp_path, "multi.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
        "chr2\t1\t.\tT\tC\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    assert len(empty_input_dict["chr1"].variant_locations) == 1
    assert len(empty_input_dict["chr2"].variant_locations) == 1


def test_parsed_variants_marked_as_input(tmp_path, ref_fasta, empty_input_dict, opts):
    """All variants from an input VCF have is_input=True."""
    vcf = _write_vcf(tmp_path, "isinput.vcf", _vcf_header_no_format() + [
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    parse_input_vcf(empty_input_dict, vcf, 2, ref_fasta, opts)
    for v in empty_input_dict["chr1"].contig_variants[0]:
        assert v.is_input is True

