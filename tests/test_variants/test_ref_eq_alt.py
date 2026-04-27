"""
Regression tests for Issue #266 — REF == ALT in output VCF.

xfail tests indicate each bug exists before fixes are applied.
Once a fix is in place, remove the corresponding xfail marker.
"""

import io
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy.random import default_rng

from neat.models.mutation_model import MutationModel
from neat.variants import SingleNucleotideVariant
from neat.variants.contig_variants import ContigVariants


# Helpers

_SEQ = "ACGTACGTACGTACGT"
_REC = SeqRecord(Seq(_SEQ), id="chr1", name="chr1", description="")


# Issue #266 — generate_snv can return a ref base as ALT

def _all_diagonal_model():
    """MutationModel whose trinuc_trans_matrices are identity (100% on diagonal)."""
    diagonal_matrix = np.eye(4)
    all_diagonal = np.stack([diagonal_matrix] * 16)
    return MutationModel(trinuc_trans_matrices=all_diagonal)


def test_generate_snv_diagonal_model_cannot_produce_ref_eq_alt_regression():
    """Asserts fixed behavior (alt != central base); xfails because bug produces alt == central."""
    model = _all_diagonal_model()
    # ACA has a central base of 'C'. With the identity matrix, rng.choice always picks 'C'.
    snv = model.generate_snv(Seq("ACA"), reference_location=5, rng=default_rng(0))
    # Fixed behavior: alt must not equal the ref base. Currently fails.
    assert snv.alt != "C", "generate_snv returned ref base as ALT (REF==ALT bug)"


def test_generate_snv_default_model_never_produces_ref_eq_alt():
    """Default model: alt must never equal the central (reference) base."""
    model = MutationModel()
    rng = default_rng(42)
    for trinuc in ["ACA", "GCG", "TAT", "CGC", "AGA", "TGT", "ACG", "GCA"]:
        central = trinuc[1]
        snv = model.generate_snv(Seq(trinuc), reference_location=10, rng=rng)
        assert snv.alt != central, (
            f"generate_snv produced REF==ALT ({central!r}) for trinuc {trinuc!r}"
        )


def test_generate_snv_diagonal_model_avoids_ref_base_after_fix():
    """After fix: diagonal custom model must not return the reference base as ALT."""
    model = _all_diagonal_model()
    rng = default_rng(0)
    for trinuc in ["ACA", "GCG", "TAT", "CGC"]:
        central = trinuc[1]
        snv = model.generate_snv(Seq(trinuc), reference_location=5, rng=rng)
        assert snv.alt != central, (
            f"generate_snv produced REF==ALT ({central!r}) for trinuc {trinuc!r}"
        )


def test_generate_snv_near_diagonal_model_avoids_ref_base_after_fix():
    """After fix: even 99%-diagonal custom model must not produce REF==ALT."""
    near_diag = np.full((4, 4), 0.01 / 3)
    np.fill_diagonal(near_diag, 0.99)
    model = MutationModel(trinuc_trans_matrices=np.stack([near_diag] * 16))
    rng = default_rng(7)
    for _ in range(200):
        snv = model.generate_snv(Seq("ACA"), reference_location=5, rng=rng)
        assert snv.alt != "C", "generate_snv produced REF==ALT with near-diagonal model"


# Issue #266 — parse_input_vcf accepts REF==ALT variants from user VCFs

def _write_vcf(tmp_path: Path, name: str, lines: list) -> Path:
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n")
    return p


def _make_opts(tmp_path):
    from neat.read_simulator.utils.options import Options
    opts = Options(rng_seed=42)
    opts.ploidy = 2
    opts.produce_vcf = True
    opts.vcf = tmp_path / "out.vcf.gz"
    return opts


def test_parse_input_vcf_rejects_ref_eq_alt_regression(tmp_path):
    """Asserts fixed behavior (variant skipped) with xfails because bug accepts REF==ALT."""
    from Bio import SeqIO
    from neat.read_simulator.utils.vcf_func import parse_input_vcf

    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{_SEQ}\n")
    ref_fasta = SeqIO.index(str(fa), "fasta")

    vcf = _write_vcf(tmp_path, "refalt.vcf", [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1\t.\tA\tA\t30\tPASS\t.",
    ])
    input_dict = {"chr1": ContigVariants()}
    parse_input_vcf(input_dict, vcf, 2, ref_fasta, _make_opts(tmp_path))
    # Fixed behavior: variant should be skipped (locations == 0). Currently fails (bug accepts it).
    assert len(input_dict["chr1"].variant_locations) == 0, \
        "REF==ALT variant was accepted — bug still present"


def test_parse_input_vcf_skips_ref_eq_alt(tmp_path):
    """After fix: parse_input_vcf skips REF==ALT variants with a warning."""
    from Bio import SeqIO
    from neat.read_simulator.utils.vcf_func import parse_input_vcf

    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{_SEQ}\n")
    ref_fasta = SeqIO.index(str(fa), "fasta")

    vcf = _write_vcf(tmp_path, "refalt.vcf", [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1\t.\tA\tA\t30\tPASS\t.",
    ])
    input_dict = {"chr1": ContigVariants()}
    parse_input_vcf(input_dict, vcf, 2, ref_fasta, _make_opts(tmp_path))
    assert len(input_dict["chr1"].variant_locations) == 0


def test_parse_input_vcf_accepts_valid_snv(tmp_path):
    """After fix: valid SNVs (REF != ALT) are still accepted normally."""
    from Bio import SeqIO
    from neat.read_simulator.utils.vcf_func import parse_input_vcf

    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{_SEQ}\n")
    ref_fasta = SeqIO.index(str(fa), "fasta")

    vcf = _write_vcf(tmp_path, "valid.vcf", [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1\t.\tA\tG\t30\tPASS\t.",
    ])
    input_dict = {"chr1": ContigVariants()}
    parse_input_vcf(input_dict, vcf, 2, ref_fasta, _make_opts(tmp_path))
    assert len(input_dict["chr1"].variant_locations) == 1


# Issue #266 — get_ref_alt / write path returns REF==ALT without any guard

def test_get_ref_alt_snv_ref_eq_alt_is_possible():
    """Demonstrate that get_ref_alt() returns REF==ALT for a badly constructed SNV.

    This is not xfail — it documents that the data model allows the condition.
    The guard must exist in write_block_vcf, not get_ref_alt itself.
    """
    # _SEQ[0] == 'A'; get_ref_alt returns ('A', 'A')
    snv = SingleNucleotideVariant(0, "A", np.array([0, 1]), 40)
    ref, alt = ContigVariants.get_ref_alt(snv, _REC, 0)
    assert str(ref) == str(alt) == "A"
