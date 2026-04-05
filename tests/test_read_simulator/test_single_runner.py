"""
Unit and integration tests for neat/read_simulator/single_runner.py

Covers:
  - initialize_all_models
  - write_block_vcf
  - read_simulator_single (integration)
"""

import gzip
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.models import (
    FragmentLengthModel,
    MutationModel,
    SequencingErrorModel,
    TraditionalQualityModel,
)
from neat.read_simulator.single_runner import (
    initialize_all_models,
    read_simulator_single,
    write_block_vcf,
)
from neat.read_simulator.utils.options import Options
from neat.read_simulator.utils.output_file_writer import OutputFileWriter
from neat.variants import ContigVariants, SingleNucleotideVariant


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_opts(tmp_path: Path, *, rng_seed: int = 0) -> Options:
    """Return a minimal Options object suitable for testing."""
    opts = Options(rng_seed=rng_seed)
    opts.paired_ended = False
    opts.produce_fastq = False
    opts.produce_vcf = False
    opts.produce_bam = False
    opts.temp_dir_path = tmp_path
    opts.reference = str(tmp_path / "ref.fa")
    opts.rng_seed = rng_seed
    opts.fq1 = None
    opts.fq2 = None
    opts.vcf = None
    opts.bam = None
    opts.mutation_model = None
    opts.error_model = None
    opts.fragment_model = None
    opts.fragment_mean = None
    opts.fragment_st_dev = None
    opts.mutation_rate = None
    opts.ploidy = 2
    opts.min_mutations = 0
    opts.mutation_bed = None
    opts.read_len = 50
    opts.coverage = 2
    return opts


def _write_ref(tmp_path: Path, seq: str = "ACGT" * 100, name: str = "chr1") -> Path:
    ref = tmp_path / "ref.fa"
    ref.write_text(f">{name}\n{seq}\n")
    return ref


def _make_vcf_ofw(tmp_path: Path) -> OutputFileWriter:
    """Return an OutputFileWriter configured for gzip VCF output."""
    opts = _make_opts(tmp_path)
    opts.produce_vcf = True
    opts.vcf = tmp_path / "out.vcf.gz"
    return OutputFileWriter(options=opts, vcf_format="gzip")


def _make_contig_variants(positions_alts: list, genotype_ploidy: int = 2) -> ContigVariants:
    """Build a ContigVariants object with SNVs at the given (pos, alt) pairs."""
    cv = ContigVariants()
    genotype = np.array([1] + [0] * (genotype_ploidy - 1))
    for pos, alt in positions_alts:
        snv = SingleNucleotideVariant(position1=pos, alt=alt, genotype=genotype, qual_score=40)
        cv.add_variant(snv)
    return cv


# ===========================================================================
# initialize_all_models — default paths (no custom model files)
# ===========================================================================

class TestInitializeAllModels:

    def test_returns_four_objects(self, tmp_path):
        opts = _make_opts(tmp_path)
        result = initialize_all_models(opts)
        assert len(result) == 4

    def test_default_mut_model_type(self, tmp_path):
        opts = _make_opts(tmp_path)
        mut_model, *_ = initialize_all_models(opts)
        assert isinstance(mut_model, MutationModel)

    def test_default_seq_error_model_type(self, tmp_path):
        opts = _make_opts(tmp_path)
        _, seq_error_model, _, _ = initialize_all_models(opts)
        assert isinstance(seq_error_model, SequencingErrorModel)

    def test_default_qual_score_model_type(self, tmp_path):
        opts = _make_opts(tmp_path)
        _, _, qual_score_model, _ = initialize_all_models(opts)
        assert isinstance(qual_score_model, TraditionalQualityModel)

    def test_default_fraglen_model_type(self, tmp_path):
        opts = _make_opts(tmp_path)
        _, _, _, fraglen_model = initialize_all_models(opts)
        assert isinstance(fraglen_model, FragmentLengthModel)

    def test_mut_model_rng_is_set(self, tmp_path):
        opts = _make_opts(tmp_path)
        mut_model, *_ = initialize_all_models(opts)
        assert mut_model.rng is opts.rng

    def test_custom_mutation_rate_overrides_model(self, tmp_path):
        opts = _make_opts(tmp_path)
        opts.mutation_rate = 0.005
        mut_model, *_ = initialize_all_models(opts)
        assert mut_model.avg_mut_rate == pytest.approx(0.005)

    def test_none_mutation_rate_does_not_override(self, tmp_path):
        opts = _make_opts(tmp_path)
        opts.mutation_rate = None
        mut_model_default, *_ = initialize_all_models(opts)
        default_rate = MutationModel().avg_mut_rate
        assert mut_model_default.avg_mut_rate == pytest.approx(default_rate)

    def test_fragment_mean_creates_fraglen_model_with_mean(self, tmp_path):
        opts = _make_opts(tmp_path)
        opts.fragment_mean = 200.0
        opts.fragment_st_dev = 40.0
        _, _, _, fraglen_model = initialize_all_models(opts)
        assert fraglen_model.fragment_mean == pytest.approx(200.0)

    def test_no_fragment_mean_uses_read_len_times_two(self, tmp_path):
        opts = _make_opts(tmp_path)
        opts.fragment_mean = None
        opts.read_len = 75
        _, _, _, fraglen_model = initialize_all_models(opts)
        assert fraglen_model.fragment_mean == pytest.approx(75 * 2.0)

    def test_fraglen_model_std_dev_set_from_fragment_mean(self, tmp_path):
        opts = _make_opts(tmp_path)
        opts.fragment_mean = None
        opts.read_len = 100
        _, _, _, fraglen_model = initialize_all_models(opts)
        expected_std = 100 * 2.0 * 0.2
        assert fraglen_model.fragment_st_dev == pytest.approx(expected_std)


# ===========================================================================
# write_block_vcf
# ===========================================================================

class TestWriteBlockVcf:

    def _ref_index(self, tmp_path: Path, seq: str = "ACGT" * 50) -> dict:
        """Write a FASTA and return a SeqIO index dict."""
        ref_path = _write_ref(tmp_path, seq=seq)
        return SeqIO.index(str(ref_path), "fasta")

    def test_single_snv_written(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = _make_contig_variants([(5, "T")])
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            content = fh.read()
        assert "chr1" in content
        assert "\t6\t" in content  # 1-based position for pos1=5

    def test_snv_ref_base_correct(self, tmp_path):
        # ref is ACGT... so position 4 is 'A' (0-based)
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = _make_contig_variants([(4, "G")])
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            content = fh.read()
        # position 4 in ACGT*50 is 'A'
        assert "\tA\t" in content

    def test_snv_alt_base_correct(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = _make_contig_variants([(4, "G")])
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            content = fh.read()
        assert "\tG\t" in content

    def test_qual_score_in_output(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = ContigVariants()
        genotype = np.array([1, 0])
        snv = SingleNucleotideVariant(position1=10, alt="C", genotype=genotype, qual_score=55)
        cv.add_variant(snv)
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            content = fh.read()
        assert "55" in content

    def test_empty_contig_variants_writes_nothing(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = ContigVariants()
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            content = fh.read()
        assert content == ""

    def test_multiple_snvs_all_written(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = _make_contig_variants([(2, "T"), (8, "C"), (15, "G")])
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            lines = [l for l in fh.readlines() if l.strip()]
        assert len(lines) == 3

    def test_variants_written_in_sorted_order(self, tmp_path):
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        # Add in reverse order; output should still be sorted
        cv = ContigVariants()
        genotype = np.array([1, 0])
        for pos, alt in [(20, "T"), (5, "G"), (12, "C")]:
            cv.add_variant(SingleNucleotideVariant(position1=pos, alt=alt, genotype=genotype, qual_score=30))
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            lines = [l for l in fh.readlines() if l.strip()]
        positions = [int(l.split("\t")[1]) for l in lines]
        assert positions == sorted(positions)

    def test_vcf_line_has_nine_tabs(self, tmp_path):
        """A VCF data line must have exactly 9 tab-separated fields (10 columns)."""
        ref_index = self._ref_index(tmp_path)
        ofw = _make_vcf_ofw(tmp_path)
        cv = _make_contig_variants([(3, "A")])
        write_block_vcf(cv, "chr1", 0, ref_index, ofw)
        ofw.flush_and_close_files(False)
        with gzip.open(tmp_path / "out.vcf.gz", "rt") as fh:
            line = fh.readline().strip()
        assert line.count("\t") == 9


# ===========================================================================
# read_simulator_single — integration tests
# ===========================================================================

class TestReadSimulatorSingle:
    """Integration tests for the full simulation loop."""

    def _run(self, tmp_path: Path, *, coverage: int = 2, read_len: int = 50,
             produce_fastq: bool = True, produce_vcf: bool = False,
             rng_seed: int = 42):
        ref_seq = "ACGT" * 100  # 400 bp
        ref_path = _write_ref(tmp_path, seq=ref_seq)

        opts = _make_opts(tmp_path, rng_seed=rng_seed)
        opts.read_len = read_len
        opts.coverage = coverage
        opts.produce_fastq = produce_fastq
        opts.produce_vcf = produce_vcf
        opts.reference = str(ref_path)

        if produce_fastq:
            opts.fq1 = tmp_path / "out.fq1.gz"
        if produce_vcf:
            opts.vcf = tmp_path / "out.vcf.gz"

        target_regions = [(0, 400, True)]
        discard_regions = [(0, 400, False)]
        mutation_regions = [(0, 400, 0.01)]

        return read_simulator_single(
            1,
            0,
            opts,
            None,
            "chr1",
            0,
            ContigVariants(),
            target_regions,
            discard_regions,
            mutation_regions,
        )

    def test_returns_four_element_tuple(self, tmp_path):
        result = self._run(tmp_path)
        assert len(result) == 4

    def test_thread_idx_preserved(self, tmp_path):
        thread_idx, *_ = self._run(tmp_path)
        assert thread_idx == 1

    def test_contig_name_preserved(self, tmp_path):
        _, contig_name, *_ = self._run(tmp_path)
        assert contig_name == "chr1"

    def test_local_variants_is_contig_variants(self, tmp_path):
        _, _, local_variants, _ = self._run(tmp_path)
        assert isinstance(local_variants, ContigVariants)

    def test_file_dict_has_fq1_key(self, tmp_path):
        _, _, _, file_dict = self._run(tmp_path)
        assert "fq1" in file_dict

    def test_file_dict_has_fq2_key(self, tmp_path):
        _, _, _, file_dict = self._run(tmp_path)
        assert "fq2" in file_dict

    def test_file_dict_has_vcf_key(self, tmp_path):
        _, _, _, file_dict = self._run(tmp_path)
        assert "vcf" in file_dict

    def test_file_dict_has_bam_key(self, tmp_path):
        _, _, _, file_dict = self._run(tmp_path)
        assert "bam" in file_dict

    def test_fq1_file_exists_after_run(self, tmp_path):
        self._run(tmp_path)
        assert (tmp_path / "out.fq1.gz").exists()

    def test_fq1_file_is_valid_gzip(self, tmp_path):
        self._run(tmp_path)
        with gzip.open(tmp_path / "out.fq1.gz", "rt") as fh:
            content = fh.read()
        assert len(content) > 0

    def test_fq1_content_contains_fastq_records(self, tmp_path):
        self._run(tmp_path)
        with gzip.open(tmp_path / "out.fq1.gz", "rt") as fh:
            lines = fh.readlines()
        # FASTQ records: 4 lines each; must have at least one record
        assert len(lines) >= 4
        # First line of first record should start with '@'
        assert lines[0].startswith("@")

    def test_reproducibility_with_same_seed(self, tmp_path):
        """Two runs with identical seeds must produce identical FASTQ output."""
        tmp_a = tmp_path / "a"
        tmp_b = tmp_path / "b"
        tmp_a.mkdir()
        tmp_b.mkdir()
        _write_ref(tmp_a)
        _write_ref(tmp_b)

        def _run_in(p):
            opts = _make_opts(p, rng_seed=99)
            opts.read_len = 50
            opts.coverage = 2
            opts.produce_fastq = True
            opts.fq1 = p / "out.fq1.gz"
            opts.reference = str(p / "ref.fa")
            return read_simulator_single(
                1, 0, opts, None, "chr1", 0, ContigVariants(),
                [(0, 400, True)], [(0, 400, False)], [(0, 400, 0.01)],
            )

        _run_in(tmp_a)
        _run_in(tmp_b)

        with gzip.open(tmp_a / "out.fq1.gz", "rb") as fa, \
             gzip.open(tmp_b / "out.fq1.gz", "rb") as fb:
            assert fa.read() == fb.read()

    def test_vcf_not_produced_when_produce_vcf_false(self, tmp_path):
        _, _, _, file_dict = self._run(tmp_path, produce_fastq=True, produce_vcf=False)
        assert file_dict["vcf"] is None