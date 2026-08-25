from neat.read_simulator.utils.options import Options

from pathlib import Path as _PathAlias
import numpy as _np
import textwrap as _textwrap
import pytest as _pytest


def _project_root() -> _PathAlias:
    return _PathAlias(__file__).resolve().parents[2]


def test_basic_options():
    reference = _project_root() / "data" / "H1N1.fa"
    base_options = Options(reference)
    assert base_options.reference == reference


def test_output_prefix_and_paths_single_end(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="neat_unit", overwrite_output=True)
    assert opts.output_dir == tmp_path
    assert opts.output_prefix == "neat_unit"
    opts.paired_ended = False
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    opts.log_configuration()
    assert opts.fq1 == tmp_path / "neat_unit.fastq.gz"
    assert opts.fq2 is None
    assert opts.bam is None
    assert opts.vcf is None


def test_output_paths_paired_end(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="pe", overwrite_output=True,
                   paired_ended=True, fragment_mean=200, fragment_st_dev=50)
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    opts.log_configuration()
    assert opts.fq1 == tmp_path / "pe_r1.fastq.gz"
    assert opts.fq2 == tmp_path / "pe_r2.fastq.gz"


def test_rng_seed_reproducible():
    opts1 = Options(rng_seed=123)
    opts2 = Options(rng_seed=123)
    a1 = opts1.rng.integers(0, 1000000, size=10)
    a2 = opts2.rng.integers(0, 1000000, size=10)
    assert (a1 == a2).all()
    opts3 = Options()
    assert isinstance(opts3.rng_seed, (int, _np.integer))


def test_from_cli_single_end_with_threads_and_splits(tmp_path: _PathAlias):
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 75
        coverage: 5
        ploidy: 2
        paired_ended: false

        produce_bam: false
        produce_vcf: false
        produce_fastq: true

        avg_seq_error: 0.01
        rescale_qualities: true
        quality_offset: 33
        rng_seed: 42
        overwrite_output: true

        parallel_block_size: 500000
        threads: 2
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli.yml"
    yml_path.write_text(cfg, encoding="utf-8")

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "fromcli", yml_path)

    assert opts.reference == _project_root() / "data" / "H1N1.fa"
    assert opts.read_len == 75
    assert opts.coverage == 5
    assert opts.ploidy == 2
    assert opts.rng_seed == 42

    assert opts.output_dir == outdir
    assert opts.output_prefix == "fromcli"
    assert opts.fq1 == outdir / "fromcli.fastq.gz"
    assert opts.fq2 is None
    assert opts.bam is None
    assert opts.vcf is None

    assert opts.threads == 2
    assert opts.splits_dir.is_dir()
    assert opts.splits_dir.name == "splits"


def test_from_cli_paired_end_fragments(tmp_path: _PathAlias):
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 101
        coverage: 2
        ploidy: 2
        paired_ended: true
        fragment_mean: 200
        fragment_st_dev: 30

        produce_bam: false
        produce_vcf: false
        produce_fastq: true

        rng_seed: 7
        overwrite_output: true

        threads: 1
        """
    ).strip() + "\n"

    yml_path = tmp_path / "neat_from_cli_pe.yml"
    yml_path.write_text(cfg, encoding="utf-8")
    outdir = tmp_path / "peout"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "peprefix", yml_path)

    assert opts.paired_ended is True
    assert opts.fragment_mean == 200
    assert opts.fragment_st_dev == 30
    assert opts.fq1 == outdir / "peprefix_r1.fastq.gz"
    assert opts.fq2 == outdir / "peprefix_r2.fastq.gz"


def test_default_values():
    opts = Options()
    assert opts.read_len == 101
    assert opts.coverage == 10
    assert opts.ploidy == 2
    assert opts.paired_ended is False
    assert opts.produce_fastq is True
    assert opts.produce_bam is False
    assert opts.produce_vcf is False
    assert opts.quality_offset == 33
    assert opts.threads == 1
    # parallel_block_size default is 0 (sentinel for auto-tune from genome length and
    # thread count). An explicit positive int in YAML overrides; see runner for the
    # auto-tune logic. The splitting strategy itself is no longer a user-facing option —
    # it's derived from `threads` at runtime.
    assert opts.parallel_block_size == 0
    assert opts.overwrite_output is False
    assert opts.rescale_qualities is False
    assert opts.min_mutations == 0
    assert opts.output_prefix == "neat_sim"
    assert opts.output_files == []
    # N handling defaults to the realistic exclude policy with a 50%-N drop threshold.
    assert opts.n_handling == "exclude"
    assert opts.n_max_fraction == 0.5


def test_rng_seed_zero():
    """Seed value 0 is valid and should not auto-generate a seed."""
    opts = Options(rng_seed=0)
    assert opts.rng_seed == 0
    # Should produce deterministic output
    a = opts.rng.integers(0, 1_000_000, size=5)
    opts2 = Options(rng_seed=0)
    b = opts2.rng.integers(0, 1_000_000, size=5)
    assert (a == b).all()


def test_copy_with_changes(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, rng_seed=1)
    new_ref = tmp_path / "other.fa"
    new_fq1 = tmp_path / "r1.fastq.gz"

    copy = opts.copy_with_changes(reference=new_ref, fq1=new_fq1)

    assert copy.reference == new_ref
    assert copy.fq1 == new_fq1
    # Unchanged fields should carry over
    assert copy.rng_seed == opts.rng_seed
    assert copy.read_len == opts.read_len
    # Original should be unmodified
    assert opts.reference == ref
    assert opts.fq1 is None


def test_copy_with_changes_no_args():
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, rng_seed=2)
    copy = opts.copy_with_changes()
    assert copy.reference == ref
    assert copy.coverage == opts.coverage


def test_check_and_log_error_none_passthrough():
    """None value should not raise or exit."""
    Options.check_and_log_error("any_key", None, 0, 100)  # no exception


def test_check_and_log_error_numeric_in_range():
    Options.check_and_log_error("coverage", 10, 1, 1000000)  # no exception


def test_check_and_log_error_numeric_out_of_range(capsys):
    with _pytest.raises(SystemExit):
        Options.check_and_log_error("coverage", 0, 1, 1000000)


def test_check_and_log_error_fractional_coverage_ok():
    """Fractional coverage is accepted under the float schema (issue #242)."""
    Options.check_and_log_error("coverage", 0.5, 1e-6, 1000000)  # no exception


def test_check_and_log_error_zero_coverage_exits():
    """Zero (or negative) coverage is rejected by the positive lower bound."""
    with _pytest.raises(SystemExit):
        Options.check_and_log_error("coverage", 0, 1e-6, 1000000)


def test_from_cli_fractional_coverage(tmp_path: _PathAlias):
    """A YAML config with fractional coverage parses and round-trips the float value."""
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 75
        coverage: 0.5
        ploidy: 2
        paired_ended: false
        produce_fastq: true
        rng_seed: 42
        overwrite_output: true
        """
    ).strip() + "\n"

    yml_path = tmp_path / "frac.yml"
    yml_path.write_text(cfg, encoding="utf-8")
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "frac", yml_path)
    assert opts.coverage == 0.5


def test_check_and_log_error_choice_valid():
    # The `choice`-type validator backs the `n_handling` schema field (exclude/telomere).
    # A synthetic key keeps this unit focused on the validator itself.
    Options.check_and_log_error("_test_choice_key", "alpha", "choice", ["alpha", "beta"])


def test_check_and_log_error_choice_invalid():
    with _pytest.raises(SystemExit):
        Options.check_and_log_error("_test_choice_key", "gamma", "choice", ["alpha", "beta"])


def test_n_handling_choice_validation():
    """The n_handling schema choices are accepted; anything else exits."""
    choices = ("exclude", "telomere")
    Options.check_and_log_error("n_handling", "exclude", "choice", choices)
    Options.check_and_log_error("n_handling", "telomere", "choice", choices)
    with _pytest.raises(SystemExit):
        Options.check_and_log_error("n_handling", "ttaggg", "choice", choices)


def test_n_max_fraction_range_validation():
    """n_max_fraction is bounded to [0.0, 1.0]."""
    Options.check_and_log_error("n_max_fraction", 0.0, 0.0, 1.0)
    Options.check_and_log_error("n_max_fraction", 1.0, 0.0, 1.0)
    with _pytest.raises(SystemExit):
        Options.check_and_log_error("n_max_fraction", 1.5, 0.0, 1.0)


def test_check_options_no_output_files_exits():
    opts = Options(rng_seed=0)
    opts.produce_fastq = False
    opts.produce_bam = False
    opts.produce_vcf = False
    with _pytest.raises(SystemExit):
        opts.check_options()


def test_check_options_paired_with_fragment_model_clears_mean_stdev():
    opts = Options(rng_seed=0, paired_ended=True,
                   fragment_model="some_model.pkl",
                   fragment_mean=300.0, fragment_st_dev=30.0)
    opts.check_options()
    assert opts.fragment_mean is None
    assert opts.fragment_st_dev is None


def test_log_configuration_produces_bam_and_vcf(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="out",
                   overwrite_output=True, produce_fastq=True,
                   produce_bam=True, produce_vcf=True)
    opts.log_configuration()
    assert opts.bam == tmp_path / "out_golden.bam"
    assert opts.vcf == tmp_path / "out_golden.vcf.gz"
    assert opts.bam in opts.output_files
    assert opts.vcf in opts.output_files


def test_read_yaml_deprecated_parallel_mode_warns(tmp_path: _PathAlias, caplog):
    """An old YAML config with `parallel_mode: ...` parses without error and emits a
    deprecation warning. The key is silently ignored (the splitting strategy is now
    derived from `threads`), so existing user configs keep working across the
    deprecation."""
    ref = _project_root() / "data" / "H1N1.fa"
    cfg = _textwrap.dedent(
        f"""
        reference: {ref}
        parallel_mode: size
        threads: 1
        """
    ).strip() + "\n"
    yml_path = tmp_path / "old_config.yml"
    yml_path.write_text(cfg, encoding="utf-8")
    import logging as _logging
    with caplog.at_level(_logging.WARNING):
        Options.from_cli(tmp_path, "out", yml_path)
    assert any(
        "parallel_mode" in rec.message and "deprecated" in rec.message
        for rec in caplog.records
    )


def test_deprecated_cleanup_splits_warns(tmp_path: _PathAlias, caplog):
    """cleanup_splits and reuse_splits must fire a deprecation warning and not crash."""
    ref = _project_root() / "data" / "H1N1.fa"
    cfg = _textwrap.dedent(
        f"""
        reference: {ref}
        cleanup_splits: true
        reuse_splits: false
        threads: 1
        """
    ).strip() + "\n"
    yml_path = tmp_path / "old_config.yml"
    yml_path.write_text(cfg, encoding="utf-8")
    import logging as _logging
    with caplog.at_level(_logging.WARNING):
        Options.from_cli(tmp_path, "out", yml_path)
    keys_warned = {rec.message.split("`")[1] for rec in caplog.records if "deprecated" in rec.message}
    assert "cleanup_splits" in keys_warned
    assert "reuse_splits" in keys_warned


def test_log_configuration_fragment_mean_less_than_read_len_exits(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="out",
                   overwrite_output=True, read_len=150,
                   fragment_mean=100.0, fragment_st_dev=10.0)
    with _pytest.raises(SystemExit):
        opts.log_configuration()


def test_log_configuration_fragment_mean_without_stdev_exits(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="out",
                   overwrite_output=True, read_len=100,
                   fragment_mean=300.0, fragment_st_dev=None)
    with _pytest.raises(SystemExit):
        opts.log_configuration()


def test_log_configuration_paired_without_model_or_mean_exits(tmp_path: _PathAlias):
    ref = _project_root() / "data" / "H1N1.fa"
    opts = Options(reference=ref, output_dir=tmp_path, output_prefix="out",
                   overwrite_output=True, paired_ended=True,
                   fragment_model=None, fragment_mean=None)
    with _pytest.raises(SystemExit):
        opts.log_configuration()




# ===========================================================================
# 3' adapter readthrough options
# ===========================================================================

def _adapter_options(**overrides):
    opts = Options(
        reference=_project_root() / "data" / "H1N1.fa",
        paired_ended=True,
        fragment_mean=80,
        fragment_st_dev=20,
    )
    for key, value in overrides.items():
        setattr(opts, key, value)
    return opts


def test_adapters_default_disabled():
    """Both short-insert features are off by default, so output is unchanged."""
    opts = Options()
    assert opts.adapters is False
    assert opts.keep_short_fragments is False
    assert opts.adapter_preset == "truseq"


def test_resolve_adapters_noop_when_disabled():
    """With neither feature on, nothing is resolved and no sequences appear."""
    opts = _adapter_options()
    opts.resolve_adapters()
    assert opts.adapter_r1 is None
    assert opts.adapter_r2 is None
    assert opts.keep_short_fragments is False


def test_resolve_adapters_truseq_preset():
    opts = _adapter_options(adapters=True, adapter_preset="truseq")
    opts.resolve_adapters()
    assert opts.adapter_r1 == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    assert opts.adapter_r2 == "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def test_resolve_adapters_nextera_preset():
    """Tagmentation puts the same Mosaic End on both reads."""
    opts = _adapter_options(adapters=True, adapter_preset="nextera")
    opts.resolve_adapters()
    assert opts.adapter_r1 == opts.adapter_r2 == "CTGTCTCTTATACACATCT"


def test_resolve_adapters_preset_overrides_user_sequences():
    """A named preset always wins, so the resolved pair can never mix the two sources."""
    opts = _adapter_options(adapters=True, adapter_preset="truseq",
                            adapter_r1="AAAA", adapter_r2="TTTT")
    opts.resolve_adapters()
    assert opts.adapter_r1 == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


def test_resolve_adapters_custom_preset_accepted():
    opts = _adapter_options(adapters=True, adapter_preset="custom",
                            adapter_r1="ACGTACGT", adapter_r2="TGCATGCA")
    opts.resolve_adapters()
    assert opts.adapter_r1 == "ACGTACGT"
    assert opts.adapter_r2 == "TGCATGCA"


@_pytest.mark.parametrize("r1,r2", [
    (None, "ACGT"),           # missing r1
    ("ACGT", None),           # missing r2
    ("acgt", "ACGT"),         # lowercase
    ("ACGTN", "ACGT"),        # N is not a definite base
    ("ACGT XYZ", "ACGT"),     # junk
])
def test_resolve_adapters_custom_preset_rejects_bad_sequences(r1, r2):
    opts = _adapter_options(adapters=True, adapter_preset="custom",
                            adapter_r1=r1, adapter_r2=r2)
    with _pytest.raises(SystemExit):
        opts.resolve_adapters()


def test_adapters_imply_keep_short_fragments():
    """Readthrough is only visible on short inserts, so enabling it must keep them."""
    opts = _adapter_options(adapters=True)
    opts.resolve_adapters()
    assert opts.keep_short_fragments is True


def test_keep_short_fragments_alone_needs_no_adapter():
    """The isolation control works on its own and resolves no adapter sequences."""
    opts = _adapter_options(keep_short_fragments=True)
    opts.resolve_adapters()
    assert opts.adapter_r1 is None
    assert opts.keep_short_fragments is True


def test_short_insert_features_require_a_fragment_distribution():
    """Short inserts have to come from somewhere; without a fragment source this is unusable."""
    opts = _adapter_options(adapters=True, fragment_mean=None, fragment_st_dev=None)
    opts.fragment_model = None
    with _pytest.raises(SystemExit):
        opts.resolve_adapters()


def test_short_fragment_mean_is_fatal_without_adapters(tmp_path: _PathAlias):
    """Below read_len every fragment would be resampled away, so the run is not what was asked for."""
    opts = Options(reference=_project_root() / "data" / "H1N1.fa", output_dir=tmp_path,
                   output_prefix="x", overwrite_output=True, paired_ended=True,
                   fragment_mean=80, fragment_st_dev=20, read_len=101)
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    with _pytest.raises(SystemExit):
        opts.log_configuration()


def test_short_fragment_mean_only_warns_with_adapters(tmp_path: _PathAlias, caplog):
    """Short-insert libraries are exactly the regime readthrough models, so this is allowed."""
    opts = Options(reference=_project_root() / "data" / "H1N1.fa", output_dir=tmp_path,
                   output_prefix="x", overwrite_output=True, paired_ended=True,
                   fragment_mean=80, fragment_st_dev=20, read_len=101, adapters=True)
    opts.produce_fastq = True
    opts.produce_bam = False
    opts.produce_vcf = False
    opts.log_configuration()
    assert any("below `read_len`" in r.message for r in caplog.records)


def _adapter_config(tmp_path: _PathAlias, extra: str) -> _PathAlias:
    cfg = _textwrap.dedent(
        f"""
        reference: {(_project_root() / 'data' / 'H1N1.fa').as_posix()}
        read_len: 101
        coverage: 2
        paired_ended: true
        fragment_mean: 80
        fragment_st_dev: 20
        produce_bam: false
        produce_vcf: false
        produce_fastq: true
        rng_seed: 42
        overwrite_output: true
        """
    ).strip() + "\n" + _textwrap.dedent(extra).strip() + "\n"
    path = tmp_path / "adapters.yml"
    path.write_text(cfg, encoding="utf-8")
    return path


def test_adapter_keys_reach_options_from_config(tmp_path: _PathAlias):
    """
    The config keys must be wired to real attributes. An option present in the `defs` schema but
    absent from Options.__init__ silently does nothing (see `no_coverage_bias`), so read the
    resolved values back off a full from_cli round trip rather than trusting the schema alone.
    """
    yml = _adapter_config(tmp_path, "adapters: true\nadapter_preset: nextera\n")
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "ad", yml)

    assert opts.adapters is True
    assert opts.adapter_preset == "nextera"
    assert opts.adapter_r1 == "CTGTCTCTTATACACATCT"
    assert opts.keep_short_fragments is True


def test_keep_short_fragments_key_reaches_options_from_config(tmp_path: _PathAlias):
    yml = _adapter_config(tmp_path, "keep_short_fragments: true\n")
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "ks", yml)

    assert opts.keep_short_fragments is True
    assert opts.adapters is False
    assert opts.adapter_r1 is None


def test_custom_adapter_sequences_from_config(tmp_path: _PathAlias):
    """Quoted DNA strings survive YAML parsing and validation intact."""
    yml = _adapter_config(
        tmp_path,
        'adapters: true\nadapter_preset: custom\nadapter_r1: "ACGTACGT"\nadapter_r2: "TGCATGCA"\n',
    )
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    opts = Options.from_cli(outdir, "cu", yml)

    assert opts.adapter_r1 == "ACGTACGT"
    assert opts.adapter_r2 == "TGCATGCA"


def test_invalid_adapter_preset_rejected(tmp_path: _PathAlias):
    yml = _adapter_config(tmp_path, "adapters: true\nadapter_preset: nonsense\n")
    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    with _pytest.raises(SystemExit):
        Options.from_cli(outdir, "bad", yml)
